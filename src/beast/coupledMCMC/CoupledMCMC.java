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
	public Input<Double> deltaTemperatureInput = new Input<>("deltaTemperature","temperature scaler, the higher this value, the hotter the chains", 0.1);
	public Input<Double> maxTemperatureInput = new Input<>("maxTemperature","temperature scaler, the higher this value, the hotter the chains");	
	public Input<Boolean> logHeatedChainsInput = new Input<>("logHeatedChains","if true, log files for heated chains are also printed", false);
	
	public Input<Boolean> optimiseTemperatureInput = new Input<>("optimiseTemperature","if specified, teh temperature is optimzed after n swaps", false);
	public Input<Integer> nrExchangesInput = new Input<>("nrExchanges","if specified, the temperature is optimzed after n swaps", 1);
	public Input<Integer> optimizeDelayInput = new Input<>("optimizeDelay","after this many iterations, the temperature will be optimized (if optimising is set to true)", 100000);
	public Input<Integer> optimizeEveryInput = new Input<>("optimizeEvery","only optimizes the temperature every n-th potential step", 1);

	public Input<Boolean> preScheduleInput = new Input<>("preSchedule","if true, how long chains are run for is scheduled at the beginning", false);
	
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

	List<StateNode> tmpStateNodes;

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
		
		// pre schedule which chains to swap when
		if (preScheduleInput.get()){
			buildSchedule();
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
//		        FileWriter outfile = new FileWriter(new File(tempDirInput.get() + stateFileName.replace("xml.state", "chain" + i + "_seed" + Randomizer.getSeed() + ".xml") ));
//		        outfile.write(sXML2);
//		        outfile.close();
				
				chains[i] = (HeatedChain) parser.parseFragment(sXML2, true);
	
				// remove all loggers
				if (!logHeatedChainsInput.get() && i != 0){
					chains[i].loggersInput.get().clear();					
				}// remove only the screen logger
				else if (i != 0 && logHeatedChainsInput.get()) {
					for (int j = 0; j < chains[i].loggersInput.get().size(); j++){
						if (chains[i].loggersInput.get().get(j).getID().contentEquals("screenlog")){
							chains[i].loggersInput.get().get(j).setIsLogging(false);;
						}
					}
				}	
				
				// initialize each chain individually
				chains[i].setResampleEvery(resampleEvery);
				chains[i].setTemperature(i, getTemperature(i));
				chains[i].setStateFile(stateFileName.replace(".state", "." + i + "state"), restoreFromFile);				
//				chains[i].setSeed(nSeed+i);
				chains[i].run();
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}		
		// get a copy of the list of state nodes to facilitate swapping states
		tmpStateNodes = startStateInput.get().stateNodeInput.get();

		chainLength = chainLengthInput.get();
		finishTimes = new long[chains.length];
		
		startOptimising = (int) optimizeDelayInput.get()/resampleEvery * nrExchangesInput.get();
	}
	
	private double getTemperature(int i){
		return i*maxTemperature/(nrOfChainsInput.get()-1);
	}
	
	class HeatedChainThread extends Thread {
		final int chainNr;

		HeatedChainThread(int chainNr) {
			this.chainNr = chainNr;
		}
		
		public void run() {
			try {
				finishTimes[chainNr] = chains[chainNr].runTillResample();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	@Override 
	public void run() throws IOException {
		
		initRun();
		
		int totalSwaps = 0;
		int successfullSwaps = 0, successfullSwaps0 = 0;
		int optimizationSteps = 0;

		for (int sampleNr = 0; sampleNr < chainLength; sampleNr += resampleEvery) {
			long startTime = System.currentTimeMillis();
			
			// start threads with individual chains here.
			threads = new Thread[chains.length];
			
			for (int k = 0; k < chains.length; k++) {
//				try {
//					chains[k].runTillResample();
//				} catch (Exception e) {
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//				}
				threads[k] = new HeatedChainThread(k);
				threads[k].start();
			}


			// wait for the chains to finish
	        startLogTime = System.currentTimeMillis();
			for (Thread thread : threads) {
				try {
					thread.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			
			if (chains.length > 1) {
				for (int ne = 0; ne < nrExchangesInput.get(); ne++){ 
					
					int i,j;
					// resample state
					i = Randomizer.nextInt(chains.length);
					
					j = i;
					while (i == j) {
						j = Randomizer.nextInt(chains.length);
					}
					
					if (i > j) {
						int tmp = i; i = j; j = tmp;
					}
										
					
					
					double p1before = chains[i].getCurrentLogLikelihood();
					double p2before = chains[j].getCurrentLogLikelihood();
	
					// robust calculations can be extremly expensive, just calculate the new probs instead 
					double p1after = chains[i].getCurrentLogLikelihood() / chains[i].getBeta() * chains[j].getBeta();
					double p2after = chains[j].getCurrentLogLikelihood() / chains[j].getBeta() * chains[i].getBeta();
					
					double logAlpha = (p1after + p2after) - (p1before  + p2before);

					if (Math.exp(logAlpha) > Randomizer.nextDouble()) {

						successfullSwaps++;
						if (i == 0) {
							successfullSwaps0++;
						}
//						System.err.println(chains[i].getCurrentLogLikelihood() / chains[i].getBeta() + " " + chains[j].getCurrentLogLikelihood() / chains[j].getBeta());

						// swap temperatures    
						double betai = chains[i].getBeta();
						double betaj = chains[j].getBeta();
						chains[j].setBeta(betai);
						chains[i].setBeta(betaj);
						// swap loggers
						swapLoggers(chains[i], chains[j]);
						
						
						chains[i].calcCurrentLogLikelihoodRobustly();
						chains[j].calcCurrentLogLikelihoodRobustly();
						

						
//						System.err.println(chains[j].getCurrentLogLikelihood() / chains[j].getBeta() + " " + chains[i].getCurrentLogLikelihood() / chains[i].getBeta());

//						System.err.println("swap " +i+ " and " + j);
//						System.exit(0);
					}
					totalSwaps++;

				}
				
				optimizationSteps++;
//				if (optimiseTemperatureInput.get() && 
//						totalSwaps > startOptimising*2 &&
//						optimizationSteps % optimizeEveryInput.get()==0){
//					// adapt the temperature scaler
//					updateTemperature(successfullSwaps, totalSwaps-startOptimising);
//					for (int k = 1; k < chains.length; k++) {
//						chains[k].setBeta(k, getTemperature(k));
//					}
//				}	
				System.err.println("succesfull swap fraction: " + (double) successfullSwaps/totalSwaps + " maximal Temperature: " + maxTemperature + " ");
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
	} // run
		
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
	void swapStates(MCMC mcmc1, MCMC mcmc2) {
		State state1 = mcmc1.startStateInput.get();
		State state2 = mcmc2.startStateInput.get();

	
		List<StateNode> stateNodes1 = state1.stateNodeInput.get();
		List<StateNode> stateNodes2 = state2.stateNodeInput.get();
		for (int i = 0; i < stateNodes1.size(); i++) {
			StateNode stateNode1 = stateNodes1.get(i);
			StateNode stateNode2 = stateNodes2.get(i);
			StateNode tmp = tmpStateNodes.get(i);
			tmp.assignFromWithoutID(stateNode1);
			stateNode1.assignFromWithoutID(stateNode2);
			stateNode2.assignFromWithoutID(tmp);
		}
	}
	
	/* swaps the states of mcmc1 and mcmc2 */
	void swapLoggers(MCMC mcmc1, MCMC mcmc2) {
		int mcmc2size = mcmc2.loggersInput.get().size();
		int mcmc1size = mcmc2.loggersInput.get().size();
		
		
		for (int i = 0; i < mcmc2size; i++){
			for (int j = 0; j < mcmc1size; j++){
				if (mcmc2.loggersInput.get().get(i).getID().contentEquals(mcmc1.loggersInput.get().get(j).getID())){
					PrintStream printstream2 = new PrintStream(mcmc2.loggersInput.get().get(i).getPrintStream());
					PrintStream printstream1 = new PrintStream(mcmc1.loggersInput.get().get(j).getPrintStream());
					
					mcmc2.loggersInput.get().get(i).setPrintStream(printstream1);
					mcmc1.loggersInput.get().get(j).setPrintStream(printstream2);
					
					boolean isLogging2 = mcmc2.loggersInput.get().get(i).getIsLogging();
					boolean isLogging1 = mcmc1.loggersInput.get().get(j).getIsLogging();
					
					mcmc2.loggersInput.get().get(i).setIsLogging(isLogging1);
					mcmc1.loggersInput.get().get(j).setIsLogging(isLogging2);
				}
					
			}
			
		}		
	}

	
	/* makes the schedule of when to swap which chains */
	void buildSchedule(){
		
	}
	
} // class MCMCMC



