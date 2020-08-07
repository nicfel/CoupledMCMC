package beast.coupledMCMC;


import beast.app.beauti.BeautiDoc;
import beast.core.*;
import beast.core.util.Log;
import beast.util.HeapSort;
import beast.util.Randomizer;
import beast.util.XMLParser;
import beast.util.XMLProducer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BetaDistributionImpl;

		
@Description("Serial temperaing aka umbrella sampling")
public class SerialMCMC extends MCMC {
	public Input<Integer> nrOfChainsInput = new Input<Integer>("chains", " number of chains to run in parallel (default 2)", 2);
	public Input<Integer> resampleEveryInput = new Input<Integer>("resampleEvery", "number of samples in between resampling (and possibly swappping) states", 100);
	public Input<String> tempDirInput = new Input<>("tempDir","directory where temporary files are written","");
	
	// input of the difference between temperature scalers
	public Input<Double> deltaTemperatureInput = new Input<>("deltaTemperature","temperature difference between the i-th and the i-th+1 chain", 0.1);
	public Input<Double> maxTemperatureInput = new Input<>("maxTemperature","maximal temperature of the hottest chain");	
	public Input<Boolean> logHeatedChainsInput = new Input<>("logHeatedChains","if true, log files for heated chains are also printed", false);
	
	// Input on whether the temperature between chains should be optimized
	public Input<Boolean> optimiseInput = new Input<>("optimise","if true, the temperature is automatically optimised to reach a target acceptance probability", true);
	
	public Input<Integer> optimiseDelayInput = new Input<>("optimiseDelay","after this many epochs/swaps, the temperature will be optimized (if optimising is set to true)", 100);
	public Input<Double> targetAcceptanceProbabilityInput = new Input<>("target", "target acceptance probability of swaps", 0.234);
	
	public Input<Boolean> useBetaDistributionInput = new Input<>("useBetaDistribution","if true, the spacing between chains is assumed to be for the quantiles of a beta distribution with alpha=1 and beta=tuneable", false);

	// nr of samples between re-arranging states
	int resampleEvery;	
	
	double deltaTemperature;
	
	private boolean optimise;
	
	/** plugins representing MCMC with model, loggers, etc **/
	HeatedChain [] chains;
	
	/** threads for running MCMC chains **/
	Thread [] threads;
	
	/** beta distribution for spacing*/
	org.apache.commons.math.distribution.BetaDistributionImpl m_dist;
	
	/** keep track of time taken between logs to estimate speed **/
    long startLogTime;
    
    /** keeps track of the last n swaps and if they were accepted or not **/
    List<Boolean> acceptedSwaps;

	// keep track of when threads finish in order to optimise thread usage
	long [] finishTimes;
		
	private ArrayList<Long>[] runTillIteration;
	
	int totalSwaps = 0;
	int successfullSwaps = 0, successfullSwaps0 = 0;
	
	long sampleOffset = 0;
	
	// defines which scheme for spacing between the heated chains to use
	private enum Spacing{Geometric, Beta };
	private Spacing spacing;

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
		
		deltaTemperature = deltaTemperatureInput.get();
		
		acceptedSwaps = new ArrayList<>();				
		
		optimise = optimiseInput.get();
		
		if (useBetaDistributionInput.get()) {
			m_dist = new BetaDistributionImpl(1, deltaTemperature);
			spacing = Spacing.Beta;
		}else {
			spacing = Spacing.Geometric;
		}
				
	} // initAndValidate
	
	private void initRun(){
		XMLProducer p = new XMLProducer();
		String sXML = p.toXML(this, new ArrayList<>());
		
		// removes coupled MCMC parts of the xml		
		sXML = sXML.replaceAll("chains=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("resampleEvery=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("tempDir=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("deltaTemperature=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("maxTemperature=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("optimise=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("logHeatedChains=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("optimizeDelay=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("optimizeEvery=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("nrExchanges=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("preSchedule=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("target=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("neighbourSwapping=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("useBetaDistribution=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("spec=\"Logger\"", "");
		sXML = sXML.replaceAll("<logger", "<coupledLogger spec=\"beast.coupledMCMC.CoupledLogger\"");
		sXML = sXML.replaceAll("</logger", "</coupledLogger");
		
		// check if the loggers have a same issue
        String sMCMCMC = this.getClass().getName();
		while (sMCMCMC.length() > 0) {
			sXML = sXML.replaceAll("\\b" + SerialMCMC.class.getName() + "\\b", HeatedChain.class.getName());
			if (sMCMCMC.indexOf('.') >= 0) {
				sMCMCMC = sMCMCMC.substring(sMCMCMC.indexOf('.') + 1);
			} else {
				sMCMCMC = "";
			}
		}
			
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
				// initialize each chain individually
				chains[i].setResampleEvery(resampleEvery);
				chains[i].setTemperature(i, getTemperature(i));
				
				// needed to avoid error of putting the working dir twice
				String[] splittedFileName = stateFileName.split("/");
				
				chains[i].setStateFile(
						splittedFileName[splittedFileName.length-1].replace(".state", "." + i + "state"), restoreFromFile);
				chains[i].setChainNr(i);
				chains[i].run();

			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}		
		// ensure that each chain has the same starting point
		threads = new Thread[chains.length];
		finishTimes = new long[chains.length];
		

		chainLength = chainLengthInput.get();
		
		
		if (restoreFromFile){
			System.out.println("restoring from file, printing to screen but not to loggers will start again from 0");
			System.out.println("we further assume that all chains ended in the same state, if logging heated chains" +
			" the different heated chains can have different amount of interations");
		}		
		
	}
	
	private double getTemperature(int i){
		if (spacing==Spacing.Geometric) {
			return i*deltaTemperature;
		}else {
			double beta=-1;
			try {
				beta = 1-m_dist.cumulativeProbability(i/(double) chains.length);
			} catch (MathException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return 1/beta -1;
		}
		
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
		
		if (restoreFromFile) {
			try {
				String str = BeautiDoc.load(new File(stateFileName));
				String [] strs = str.split("\n");
				for (int i = 0; i < strs.length; i++) {
					String [] Spltstr = strs[i].replaceAll("\\s", "").split("=");
					if (Spltstr[0].contentEquals("sample"))
						sampleOffset = Long.parseLong(Spltstr[1]);
					else if (Spltstr[0].contentEquals("totalSwaps"))
						totalSwaps = Integer.parseInt(Spltstr[1]);
					else if (Spltstr[0].contentEquals("successfullSwaps"))
						successfullSwaps = Integer.parseInt(Spltstr[1]);
					else if (Spltstr[0].contentEquals("successfullSwaps0"))
						successfullSwaps0 = Integer.parseInt(Spltstr[1]);
					else if (Spltstr[0].contentEquals("deltaTemperature"))
						deltaTemperature = Double.parseDouble(Spltstr[1]);
					else if (Spltstr[0].contentEquals("lastSwaps")) {
						String[] tmp = Spltstr[1].replaceAll("\\[", "").replaceAll("\\]", "").replaceAll("\\s", "").split(",");
						acceptedSwaps = new ArrayList<>();
						for (int j = 0; j < tmp.length; j++)
							acceptedSwaps.add(Boolean.parseBoolean(tmp[j]));
					}					
				}			
				
				Log.warning("Restoring: totalSwaps=" + totalSwaps + 
						" successfullSwaps=" + successfullSwaps +
						" successfullSwaps0=" + successfullSwaps0 +
						" deltaTemperature=" + deltaTemperature);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		
		runNeigbours();
	} // run
	
	private void runNeigbours(){
		// run each thread until it's next swapping time
		// start threads with individual chains here.
		threads = new Thread[chains.length];
		for (int k = 0; k < chains.length; k++) {
			threads[k] = new HeatedChainThread(k, resampleEvery);
			threads[k].start();
		}
		
		startLogTime = -1;
		long startSample = 0;
		
		// print header for system output
		System.out.println("sample\tswapsColdCain\tswapProbability\tdeltaTemperature");		
		double currProb = ((double) successfullSwaps/totalSwaps);
		if (Double.isNaN(currProb))
			currProb = 0.0;
		System.out.print("\t" + (startSample + sampleOffset) + "\t" + successfullSwaps0 + "\t" + currProb + "\t" + deltaTemperature + "\n");

		
		for (long sampleNr = resampleEvery; sampleNr <= chainLength; sampleNr += resampleEvery) {	
			// get the chains to swap, conditioning on them being neighbours
			int chain_i = Randomizer.nextInt(chains.length-1);
			int i=-1,j=-1;
			
			for (int k = 0; k<chains.length; k++) {
				if (chains[k].getChainNr()==chain_i)
					i=k;
				if (chains[k].getChainNr()==chain_i+1)
					j=k;	
			}
						
			// look that every chain goes to the next place 
//			for (int k = 0; k < chains.length; k++) {
			{	int k = 0;
				try {
					threads[k].join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
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
				
				// swap chain numbers
				int chainNr = chains[i].getChainNr();
				chains[i].setChainNr(chains[j].getChainNr());
				chains[j].setChainNr(chainNr);
				
				// swap loggers and the state file names
				swapLoggers(chains[i], chains[j]);				
				
				// swap Operator tuning
				swapOperatorTuning(chains[i], chains[j]);		
				
				if (chains[j].getBeta() < chains[i].getBeta())
					throw new IllegalArgumentException("error in temperatures of chains");
				
				if (optimise) {
					acceptedSwaps.add(true);
				}
			} else {
				if (optimise) {
					acceptedSwaps.add(false);
				}
			}
			totalSwaps++;
			
			if (optimise && totalSwaps > optimiseDelayInput.get()) {
				double delta = getDelta();			
				
				deltaTemperature += delta;
	            
	            // boundary case checks
	    		if (maxTemperatureInput.get() != null){
	    			deltaTemperature = Math.max(deltaTemperature, maxTemperatureInput.get()/(chains.length-1)); 
	            } else if (deltaTemperature < 0) {
	            	deltaTemperature = 0;
	            }
	            
	            // figure out order of chains
	            int n = chains.length;
	            int [] order = new int[n];
	            for (int k = 0; k < n; k++) {
	            	order[k] = k;
	            }
	            double [] temp = new double[chains.length];
	            for (int k = 0; k < n; k++) {
	            	temp[k] = 1.0 - chains[k].getBeta();
	            }
	            HeapSort.sort(temp, order);
	            
	            // if the spacing is using a beta distribution, update the beta distribution
	            if (spacing==Spacing.Beta)
	            	m_dist = new BetaDistributionImpl(1, deltaTemperature);
	            
	            // set new temperatures asynchronously, in same order as before
	            for (int k = 0; k < n; k++) {
	            	chains[k].setTemperature(0, getTemperature(chains[k].getChainNr()));
	            }
	            
			}

			if (sampleNr < chainLength) {
				//for (int k = 0; k < chains.length; k++) {
				int k = 0;
					threads[k] = new HeatedChainThread(k, sampleNr+resampleEvery);
					threads[k].start();
				//}
			}

			
			System.out.print("\t" + (sampleNr + sampleOffset) + "\t" + successfullSwaps0 + "\t" + ((double) successfullSwaps/totalSwaps) + "\t" + deltaTemperature + " ");
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

			if (sampleNr % storeEveryInput.get() == 0) {
		        try {
		        	storeStateToFile(sampleNr);
		        } catch (Exception e) {
		            e.printStackTrace();
		        }
			}			
		}
		// ensure that every chains ran to the end even if it's not participating in the last swap
//		for (int i = 0; i < threads.length; i++){
//			try {
//				threads[i].join();
//			} catch (InterruptedException e) {
//				e.printStackTrace();
//			}
//		}
		
		System.err.println("#Swap attemps = " + totalSwaps);
		System.err.println("#Successfull swaps = " + successfullSwaps);
		System.err.println("#Successfull swaps with cold chain = " + successfullSwaps0);
		// wait 5 seconds for the log to complete
		try {
			Thread.sleep(5000);
		} catch (InterruptedException e) {
			// ignore
		}
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
				}					
			}			
		}		
				
		// swap the state file Names as well
		String stateFileName1 = mcmc1.getStateFileName();
		mcmc1.setStateFileName(mcmc2.getStateFileName());
		mcmc2.setStateFileName(stateFileName1);		
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
		    
		    
		    double coercableParameterValue = operator1.getCoercableParameterValue();
		    operator1.setCoercableParameterValue(operator2.getCoercableParameterValue());
		    operator2.setCoercableParameterValue(coercableParameterValue);
			
		}		
	}
	
	
	void storeStateToFile(long currentSample) throws FileNotFoundException {
        PrintStream out;
		try {
			out = new PrintStream(stateFileName + ".new");
			out.print("sample=" + (currentSample + sampleOffset) + "\n");
	        out.print("totalSwaps=" + totalSwaps + "\n");
	        out.print("successfullSwaps=" + successfullSwaps + "\n");
	        out.print("successfullSwaps0=" + successfullSwaps0 + "\n");
	        out.print("deltaTemperature=" + deltaTemperature + "\n");
        	out.print("lastSwaps=" + acceptedSwaps.toString() + "\n");
	
	        for (int k = 0; k < chains.length; k++) {
	        	out.print("beta_chain." + k + "=" + chains[k].getBeta() + "\n");
	        }
	        //out.print(new XMLProducer().toXML(this));
	        out.close();
	        File newStateFile = new File(stateFileName + ".new");
	        File oldStateFile = new File(stateFileName);
	        oldStateFile.delete();
	        // newStateFile.renameTo(oldStateFile); -- unstable under windows
	        Files.move(newStateFile.toPath(), oldStateFile.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	double getDelta() {
		double target = targetAcceptanceProbabilityInput.get();	

		double p = (double) successfullSwaps / totalSwaps;
		double p_local = 0.0;
		
		// compute the local acceptance probability
		acceptedSwaps.remove(0);
		// compute local acceptance probability
		for (int k = acceptedSwaps.size()/10; k < acceptedSwaps.size(); k++)
			if (acceptedSwaps.get(k))
				p_local += 1;
		
		p_local /= acceptedSwaps.size();
		
		
		
		// check that the local and global acceptance probability are at the same relative position
		// compared to the target acceptance probability
		if ((p<target) && (p_local>target))
			return 0.0;
		else if ((p>target) && (p_local<target))
			return 0.0;		
		

		double swapsTransformed;
		
//		if (optimise==Optimise.Log) {
//			swapsTransformed = Math.log(totalSwaps + 1.0);
//		}else if (optimise==Optimise.Sqrt) {
//			swapsTransformed = Math.sqrt(totalSwaps);			
//		}else {
			swapsTransformed = (double) totalSwaps;		
//		}
				
		// update maxTemperature		
		double delta = (p - target) / swapsTransformed;
		
		// prevent too large adaptions that lead to overshooting
		double maxstep = 1.0/(optimiseDelayInput.get()*10.0);
		
		if (delta>maxstep)
			delta = maxstep;
		if (delta<-maxstep)
			delta = -maxstep;
		
		
		return delta;
	}
	
}



