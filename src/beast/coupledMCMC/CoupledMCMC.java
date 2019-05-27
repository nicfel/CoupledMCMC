package beast.coupledMCMC;


import beast.core.*;
import beast.core.util.Log;
import beast.util.Randomizer;
import beast.util.XMLParser;
import beast.util.XMLProducer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

@Citation(value= "Altekar G, Dwarkadas S, Huelsenbeck J and Ronquist F (2004). \n" +
				"  Parallel Metropolis Coupled Markov Chain Monte Carlo For Bayesian Phylogenetic Inference.\n" +
				"  Bioinformatics, 20(3), 407-415."
		, year = 2004, firstAuthorSurname = "Altekar",
		DOI="10.1093/bioinformatics/btg427")
@Description("")
public class CoupledMCMC extends MCMC {
	public Input<Integer> nrOfChainsInput = new Input<Integer>("chains", " number of chains to run in parallel (default 2)", 2);
	public Input<Integer> resampleEveryInput = new Input<Integer>("resampleEvery", "number of samples in between resampling (and possibly swappping) states", 10000);
	public Input<String> tempDirInput = new Input<>("tempDir","directory where temporary files are written","");
	
	// input of the difference between temperature scalers
	public Input<Double> deltaTemperatureInput = new Input<>("deltaTemperature","temperature difference between the i-th and the i-th+1 chain", 0.1);
	public Input<Double> maxTemperatureInput = new Input<>("maxTemperature","temperature scaler, the higher this value, the hotter the chains");	
	public Input<Boolean> logHeatedChainsInput = new Input<>("logHeatedChains","if true, log files for heated chains are also printed", true);
	
//	public Input<Boolean> optimiseTemperatureInput = new Input<>("optimiseTemperature","if specified, teh temperature is optimzed after n swaps", false);
//	public Input<Integer> nrExchangesInput = new Input<>("nrExchanges","if specified, the temperature is optimzed after n swaps", 1);
//	public Input<Integer> optimizeDelayInput = new Input<>("optimizeDelay","after this many iterations, the temperature will be optimized (if optimising is set to true)", 100000);
//	public Input<Integer> optimizeEveryInput = new Input<>("optimizeEvery","only optimizes the temperature every n-th potential step", 1);

	public Input<Boolean> preScheduleInput = new Input<>("preSchedule","if true, how long chains are run for is scheduled at the beginning", true);
	
	// nr of samples between re-arranging states
	int resampleEvery;	
	double maxTemperature;
	int startOptimising;	
	
	/** plugins representing MCMC with model, loggers, etc **/
	HeatedChain [] chains;
	/** threads for running MCMC chains **/
	Thread [] threads;
	/** keep track of time taken between logs to estimate speed **/
    long startLogTime;

	// keep track of when threads finish in order to optimise thread usage
	long [] finishTimes;
		
	private ArrayList<Long>[] runTillIteration;


	@Override
	public void initAndValidate() {
		if (nrOfChainsInput.get() < 1) {
			throw new RuntimeException("chains must be at least 1");
		}
		if (nrOfChainsInput.get() == 1) {
			Log.warning.println("Warning: coupled MCMC needs at least 2 chains, but the number of chains is 1. Running plain MCMC.");
		}
		// initialize the differently heated chains
		chains = new HeatedChain[nrOfChainsInput.get()];
		
		resampleEvery = resampleEveryInput.get();				
		
		if (maxTemperatureInput.get()!=null){
			System.err.println("calculating delta T from max Temperature and number of chains");
			maxTemperature = maxTemperatureInput.get();
		}else{
			maxTemperature = deltaTemperatureInput.get()*(nrOfChainsInput.get()-1);
		}
		
	} // initAndValidate
	
	private void initRun(){
		String sXML = new XMLProducer().toXML(this);
		
		// removes coupled MCMC parts of the xml
		
		
		sXML = sXML.replaceAll("chains=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("resampleEvery=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("tempDir=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("deltaTemperature=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("maxTemperature=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("optimiseTemperature=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("logHeatedChains=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("optimizeDelay=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("optimizeEvery=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("nrExchanges=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("preSchedule=['\"][^ ]*['\"]", "");
		
		sXML = sXML.replaceAll("spec=\"Logger\"", "");
		sXML = sXML.replaceAll("<logger", "<coupledLogger spec=\"beast.coupledMCMC.CoupledLogger\"");
		sXML = sXML.replaceAll("</logger", "</coupledLogger");
		
		// check if the loggers have a same issue
        String sMCMCMC = this.getClass().getName();
		while (sMCMCMC.length() > 0) {
			sXML = sXML.replaceAll("\\b" + CoupledMCMC.class.getName() + "\\b", HeatedChain.class.getName());
			if (sMCMCMC.indexOf('.') >= 0) {
				sMCMCMC = sMCMCMC.substring(sMCMCMC.indexOf('.') + 1);
			} else {
				sMCMCMC = "";
			}
		}
//		long nSeed = Randomizer.getSeed();
			
		// create new chains		
		for (int i = 0; i < chains.length; i++) {
			XMLParser parser = new XMLParser();
			String sXML2 = sXML;
			if (i>0){
				sXML2 = sXML2.replaceAll("fileName=\"", "fileName=\"chain" + i+ "");
			}
			
			try {
				chains[i] = (HeatedChain) parser.parseFragment(sXML2, true);
	
				// remove all screen loggers
				for (int j = chains[i].coupledLoggersInput.get().size()-1; j >=0 ; j--){
					if (chains[i].coupledLoggersInput.get().get(j).getID().contentEquals("screenlog")){
						chains[i].coupledLoggersInput.get().remove(j);
					}
				}
				
				// remove all loggers of heated chains if they are not logged
				if (!logHeatedChainsInput.get() && i != 0){
					for (int j = 0; j < chains[i].coupledLoggersInput.get().size(); j++)
						chains[i].coupledLoggersInput.get().get(j).setSuppressLogging(true);					
				}
				// remove only the screen logger
				
				// initialize each chain individually
				chains[i].setResampleEvery(resampleEvery);
				chains[i].setTemperature(i, getTemperature(i));
				
				// needed to avoid error of putting the working dir twice
				String[] splittedFileName = stateFileName.split("/");
				
				
				chains[i].setStateFile(
						splittedFileName[splittedFileName.length-1].replace(".state", "." + i + "state"), restoreFromFile);				
				chains[i].setChainNr(i);
				chains[i].run();
				
//				System.out.println(chains[0].coupledLoggersInput.get().get(0).fileName);
//				System.out.println(chains[0].getStateFileName());
//				System.exit(0);


			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}		
		// ensure that each chain has the same starting point
		threads = new Thread[chains.length];
		finishTimes = new long[chains.length];
		
//		for (int k = 0; k < chains.length; k++) {			
//			threads[k] = new HeatedChainThread(k, maxIteration);
//			threads[k].start();
//		}
//		
//		System.out.println(maxIteration);
//		System.exit(0);

		chainLength = chainLengthInput.get();
		
//		startOptimising = (int) optimizeDelayInput.get()/resampleEvery * nrExchangesInput.get();
		
		// pre schedule which chains to swap when
		if (preScheduleInput.get()){
			buildSchedule();
		}
		
		if (restoreFromFile){
			System.out.println("restoring from file, printing to screen but not to loggers will start again from 0");
			System.out.println("we further assume that all chains ended in the same state, if logging heated chains" +
			" the different heated chains can have different amount of interations");
		}
		
		
		System.out.println("sample\tswapsColdCain\tswapProbability");

		
	}
	
	private double getTemperature(int i){
		return i*maxTemperature/(nrOfChainsInput.get()-1);
	}
	
	class HeatedChainThread extends Thread {
		final int chainNr;
		long runUntil;

		HeatedChainThread(int chainNr, long runUntil) {
			this.chainNr = chainNr;
			this.runUntil = runUntil;
		}		
		
		public void run() {
			try {
				finishTimes[chainNr] = chains[chainNr].runTillResample(runUntil);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	@Override 
	public void run() throws IOException {
		
		initRun();
		
		if (preScheduleInput.get()){
			runPrescheduled();
			return;
		}
	} // run
	
	private void runPrescheduled(){
		
		int totalSwaps = 0;
		int successfullSwaps = 0, successfullSwaps0 = 0;
		
		

		
		// run each thread until it's next swapping time
		// start threads with individual chains here.
		threads = new Thread[chains.length];
		for (int k = 0; k < chains.length; k++) {
			
			long runk = chainLength;
			
			if (runTillIteration[k].size()>0)  runk = runTillIteration[k].get(0);

			
			threads[k] = new HeatedChainThread(k, runk);
			threads[k].start();
		}
		
		startLogTime = -1;
		long startSample = 0;
		
		for (long sampleNr = resampleEvery; sampleNr < chainLength; sampleNr += resampleEvery) {
			
				
			
			// get the chains to swap
			int i=-1, j=-1;
			for (int k = 0; k < threads.length; k++){
				if (runTillIteration[k].size()>0){
					if (runTillIteration[k].get(0) == sampleNr){
						if (i>=0)
							j = k;
						else
							i = k;
					}
				}
			}
			
			try {
				threads[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			
			try {
				threads[j].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			

					
					
			if (chains[i].getBeta() <  chains[j].getBeta()) {
				int tmp = i; i = j; j = tmp;
			}
		
			
			double p1before = chains[i].getCurrentLogLikelihood();
			double p2before = chains[j].getCurrentLogLikelihood();

			// robust calculations can be extremly expensive, just calculate the new probs instead 
			double p1after = chains[i].getUnscaledCurrentLogLikelihood() * chains[j].getBeta();
			double p2after = chains[j].getUnscaledCurrentLogLikelihood() * chains[i].getBeta();
			
			double logAlpha = (p1after + p2after) - (p1before  + p2before);
			if (Math.exp(logAlpha) > Randomizer.nextDouble()) {

				successfullSwaps++;
				if (i == 0) {
					successfullSwaps0++;
				}

				// swap temperatures    
				double beta = chains[i].getBeta();
				chains[i].setBeta(chains[j].getBeta());						
				chains[j].setBeta(beta);
				
				// swap loggers and the state file names
				swapLoggers(chains[i], chains[j]);				
				
				// swap Operator tuning
				swapOperatorTuning(chains[i], chains[j]);
			}
			totalSwaps++;

			runTillIteration[i].remove(0);
			runTillIteration[j].remove(0);
			
			long runi = chainLength, runj = chainLength;
			
			if (runTillIteration[i].size()>0)  runi = runTillIteration[i].get(0);
			if (runTillIteration[j].size()>0)  runj = runTillIteration[j].get(0);
			
			threads[i] = new HeatedChainThread(i, runi);
			threads[i].start();
			
			threads[j] = new HeatedChainThread(j, runj);
			threads[j].start();
			
			System.out.print("\t" + sampleNr + "\t" + successfullSwaps0 + "\t" + ((double) successfullSwaps/totalSwaps) + " ");
			if (startLogTime>0){			
	            final long logTime = System.currentTimeMillis();
	            final int secondsPerMSamples = (int) ((logTime - startLogTime) * 1000.0 / (sampleNr - startSample + 1.0));
	            final String timePerMSamples =
	                    (secondsPerMSamples >= 3600 ? secondsPerMSamples / 3600 + "h" : "") +
	                            (secondsPerMSamples >= 60 ? (secondsPerMSamples % 3600) / 60 + "m" : "") +
	                            (secondsPerMSamples % 60 + "s");

            
	            System.out.print(timePerMSamples + "/Msamples\n");
			}else{
	            System.out.print("--\n");
			}
			if (sampleNr>=10000 && startLogTime<0){
				startSample = sampleNr;
				startLogTime = System.currentTimeMillis();
			}

				
//			System.err.println(sampleNr + " " + i + " " + j + " succesfull swap fraction: " + (double) successfullSwaps/totalSwaps + " maximal Temperature: " + maxTemperature + " ");
			
		}
		// ensure that every chains ran to the end even if it's not participating in the last swap
		for (int i = 0; i < threads.length; i++){
			try {
				threads[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		

		System.err.println("#Successfull swaps = " + successfullSwaps);
		System.err.println("#Successfull swaps with cold chain = " + successfullSwaps0);
		// wait 5 seconds for the log to complete
		try {
			Thread.sleep(5000);
		} catch (InterruptedException e) {
			// ignore
		}

	}
		
    public double updateTemperature(int successfullSwaps, int totalSwaps) {
        double successfullFraction = (double) successfullSwaps/totalSwaps;
        final double target = 0.35;

        double count = (totalSwaps + 1.0);
//        switch (transform) {
//            case log:
//                double count = Math.log(totalSwaps + 1.0);
//                break;
//            case sqrt:
//                double count = Math.sqrt(totalSwaps);
//                break;
//            case none:
//            	break;
//            default:
//            	break;
//        }

        double deltaP = successfullFraction - target;
    	System.out.println(successfullFraction + " " + target + " " + Math.exp(1/count * deltaP));

        
  
        maxTemperature *= Math.exp(1/count * deltaP);
        
//        maxTemperature*=0.999999;
        maxTemperature = Math.min(maxTemperature, 1.0);
//        deltaP += Math.log(1.0 / temperatureScaler - 1.0);        
        
        if (deltaP > -Double.MAX_VALUE && deltaP < Double.MAX_VALUE) {
            return deltaP;
        }
        return 0;
    }
	
	/* swaps the states of mcmc1 and mcmc2 */
	void swapLoggers(HeatedChain mcmc1, HeatedChain mcmc2) {
		int mcmc2size = mcmc2.coupledLoggersInput.get().size();
		int mcmc1size = mcmc1.coupledLoggersInput.get().size();
		
		for (int i = 0; i < mcmc2size; i++){
			for (int j = 0; j < mcmc1size; j++){
				if (mcmc2.coupledLoggersInput.get().get(i).getID().contentEquals(mcmc1.coupledLoggersInput.get().get(j).getID())){
					if (!logHeatedChainsInput.get()){
						boolean suppressLogging1 = mcmc1.coupledLoggersInput.get().get(i).getSuppressLogging();
						boolean suppressLogging2 = mcmc2.coupledLoggersInput.get().get(j).getSuppressLogging();
						
						mcmc1.coupledLoggersInput.get().get(i).setSuppressLogging(suppressLogging2);
						mcmc2.coupledLoggersInput.get().get(j).setSuppressLogging(suppressLogging1);

					}

					PrintStream printstream2 = mcmc2.coupledLoggersInput.get().get(i).getM_out();
					PrintStream printstream1 = mcmc1.coupledLoggersInput.get().get(j).getM_out();
										
					mcmc2.coupledLoggersInput.get().get(i).setM_out(printstream1);
					mcmc1.coupledLoggersInput.get().get(j).setM_out(printstream2);
//					printstream1.close();
//					printstream2.close();
				}					
			}			
		}		
		
		// swap the state file Names as well
		String stateFileName = mcmc1.getStateFileName();
		mcmc1.setStateFileName(mcmc2.getStateFileName());
		mcmc2.setStateFileName(stateFileName);
	}


	void swapOperatorTuning(HeatedChain mcmc1, HeatedChain mcmc2) {
		List<Operator> operatorList1 = new ArrayList<Operator>();
		List<Operator> operatorList2 = new ArrayList<Operator>();		
		
		operatorList1 = mcmc1.getOperatorSchedule().operatorsInput.get();
		operatorList2 = mcmc2.getOperatorSchedule().operatorsInput.get();
		
		
		for (int i = 0; i < operatorList1.size(); i++){
			Operator operator1 = operatorList1.get(i);
			Operator operator2 = operatorList2.get(i);

		    int m_nNrRejected = operator1.get_m_nNrRejected();
		    int m_nNrAccepted = operator1.get_m_nNrAccepted();
		    int m_nNrRejectedForCorrection = operator1.get_m_nNrRejectedForCorrection();
		    int m_nNrAcceptedForCorrection = operator1.get_m_nNrAcceptedForCorrection();
		    
		    operator1.setAcceptedRejected(operator2.get_m_nNrRejected(), operator2.get_m_nNrAccepted(), operator2.get_m_nNrRejectedForCorrection(), operator2.get_m_nNrAcceptedForCorrection());
		    operator2.setAcceptedRejected(m_nNrAccepted, m_nNrRejected, m_nNrAcceptedForCorrection, m_nNrRejectedForCorrection);
		    
//		    System.out.println(operator1.getCoercableParameterValue() + " " + operator2.getCoercableParameterValue());
		    
		    double coercableParameterValue = operator1.getCoercableParameterValue();
		    operator1.setCoercableParameterValue(operator2.getCoercableParameterValue());
		    operator2.setCoercableParameterValue(coercableParameterValue);
		    
//		    System.out.println(operator1.getCoercableParameterValue() + " " + operator2.getCoercableParameterValue());
//		    System.out.println();


			
			
		}		
	}


	
	/* makes the schedule of when to swap which chains */
	@SuppressWarnings("unchecked")
	void buildSchedule(){
		// initialize Arrays
		runTillIteration = new ArrayList[chains.length];
		for (int i = 0; i < chains.length; i++){
			runTillIteration[i] = new ArrayList<Long>();
		}

		for (long sampleNr = resampleEvery; sampleNr < chainLength; sampleNr += resampleEvery) {
			int i,j;
			// resample state
			i = Randomizer.nextInt(chains.length);
			
			j = i;
			while (i == j) {
				j = Randomizer.nextInt(chains.length);
			}
			runTillIteration[i].add(sampleNr);
			runTillIteration[j].add(sampleNr);
		}
	}
	
} // class MCMCMC



