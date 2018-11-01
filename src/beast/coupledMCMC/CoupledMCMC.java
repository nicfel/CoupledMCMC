package beast.coupledMCMC;


import beast.core.*;
import beast.core.util.Log;
import beast.util.Randomizer;
import beast.util.XMLParser;
import beast.util.XMLProducer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

// Altekar G, Dwarkadas S, Huelsenbeck J and Ronquist F (2004). 
// Parallel Metropolis Coupled Markov Chain Monte Carlo For Bayesian Phylogenetic Inference.
// Bioinformatics, 20. ISSN 1367-4803, 
// http://dx.doi.org/10.1093/bioinformatics/btg427.

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
	
	public Input<Boolean> optimiseTemperatureInput = new Input<>("optimiseTemperature","if specified, teh temperature is optimzed after n swaps", true);
	public Input<Integer> nrExchangesInput = new Input<>("nrExchanges","if specified, the temperature is optimzed after n swaps", 1);
	public Input<Integer> acceptedWindowInput = new Input<>("acceptedWindow","window for which to calculate the acceptance ratio", 1000);
	public Input<Integer> stopUpdatingTemperatureInput = new Input<>("stopUpdatingTemperature","window for which to calculate the acceptance ratio", 100000);
	
	// nr of samples between re-arranging states
	int resampleEvery;	
	double maxTemperature;
	
	// keeps track of 
	boolean[] wasAccepted;
	int wacount;
	
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
			Log.warning.println("Warning: MCMCMC needs at least 2 chains to be effective, but chains=1. Running plain MCMC.");
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
		
		wasAccepted = new boolean[acceptedWindowInput.get()];
		wacount = 0;
				

		
	} // initAndValidate
	
	private void initRun(){
		String sXML = new XMLProducer().toXML(this);
		// removes coupled MCMC parts of the xml
		sXML = sXML.replaceAll("chains=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("resampleEvery=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("tempDir=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("deltaTemperature=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("maxTemperature=['\"][^ ]*['\"]", "");
		
		sXML = sXML.replaceAll("logHeatedChains=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("adaptTemperatureEvery=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("swapNeighboursOnly=['\"][^ ]*['\"]", "");
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
		        FileWriter outfile = new FileWriter(new File(tempDirInput.get() + stateFileName.replace("xml.state", "chain" + i + ".xml") ));
		        outfile.write(sXML2);
		        outfile.close();
				
				chains[i] = (HeatedChain) parser.parseFragment(sXML2, true);
	
				// remove all loggers
				if (!logHeatedChainsInput.get() && i != 0){
					chains[i].loggersInput.get().clear();					
				}// remove only the screen logger
				else if (i != 0 && logHeatedChainsInput.get()) {
					for (int j = 0; j < chains[i].loggersInput.get().size(); j++){
						if (chains[i].loggersInput.get().get(j).getID().contentEquals("screenlog")){
							chains[i].loggersInput.get().remove(j);
						}
					}
				}	
				
				// initialize each chain individually
				chains[i].setChainNr(i, resampleEvery, getTemperature(i));
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
	}
	
	private double getTemperature(int i){
//		return (1-Math.exp(-i))*maxTemperature;
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
//					System.err.println(successfullSwaps0 + " " + successfullSwaps + ": " + i + " <--> " + j + ": " + logAlpha +  ": " + ((double) successfullSwaps/totalSwaps) + ": " + deltaTemperature);
					if (Math.exp(logAlpha) > Randomizer.nextDouble()) {
						wasAccepted[wacount] = true;
						wacount++;
						if (wacount==wasAccepted.length)
							wacount = 0;
												
						if (i == 0) {
							successfullSwaps0++;
						}
//						System.err.println(i + " <--> " + j);
						swapStates(chains[i], chains[j]);
						chains[i].calcCurrentLogLikelihoodRobustly();
						chains[j].calcCurrentLogLikelihoodRobustly();
	
						//assignState(chains[j], chains[i]);
						//chains[j].calcCurrentLogLikelihoodRobustly();
					}else{
						wasAccepted[wacount] = false;
						wacount++;
						if (wacount==wasAccepted.length)
							wacount = 0;
						
					}
					totalSwaps++;
				}
				if (totalSwaps < stopUpdatingTemperatureInput.get()){
					if (optimiseTemperatureInput.get() && 
							totalSwaps > wasAccepted.length){
						// adapt the temperature scaler
						updateTemperature();
						for (int k = 1; k < chains.length; k++) {
							chains[k].setBeta(k, getTemperature(k));
						}
						 
						
					}
				}
				System.err.println("succesfull swap fraction: " + getSuccesfullFraction() + " maximal Temperature: " + maxTemperature);
//				for (int k = 0; k < chains.length; k++) {
//					System.err.print(chains[k].getBeta() + " ");
//				}
//				System.exit(0);
				
				
					
//				// tuning
//				for (int k = 1; k < chains.length; k++) {
//					chains[k].optimiseRunTime(startTime, finishTimes[k], finishTimes[0]);
//				}
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
	
	public double getSuccesfullFraction(){
		int successfullCount = 0;
    	for (int i = 0; i < wasAccepted.length; i++)
    		if (wasAccepted[i])
    			successfullCount++;
    	return (double) successfullCount /(double) wasAccepted.length;
		
	}
	
    public double updateTemperature() {
    	int totalSwaps = wasAccepted.length; 
    	double successfullFraction = getSuccesfullFraction();
        
        final double target = 0.35;

//        double count = (totalSwaps + 1.0);
//        switch (transform) {
//            case log:
                double count = Math.log(totalSwaps + 1.0);
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
  
//        maxTemperature *= Math.exp(1/count * deltaP);
        
        maxTemperature*=0.999999;
//        maxTemperature = Math.min(maxTemperature, 1.0);
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
	
	
} // class MCMCMC



