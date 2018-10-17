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
@Description("Metropolis-Coupled Markov Chain Monte Carlo" +
		"" +
		"Note that log file names should have $(seed) in their name so " +
		"that the first chain uses the actual seed in the file name and all subsequent chains add one to it." +
		"Furthermore, the log and tree log should have the same sample frequency.")
public class CoupledMCMC extends MCMC {
	public Input<Integer> nrOfChainsInput = new Input<Integer>("chains", " number of chains to run in parallel (default 2)", 2);
	public Input<Integer> resampleEveryInput = new Input<Integer>("resampleEvery", "number of samples in between resampling (and possibly swappping) states", 10000);
	public Input<String> tempDirInput = new Input<>("tempDir","directory where temporary files are written","");
	public Input<Double> temperatureScalerInput = new Input<>("temperatureScaler","temperature scaler, the higher this value, the hotter the chains", 0.01);
	public Input<Boolean> logHeatedChainsInput = new Input<>("logHeatedChains","if true, log files for heated chains are also printed", false);
	public Input<Integer> adaptTemperatureEveryInput = new Input<>("adaptTemperatureEvery","if specified, teh temperature is optimzed after n swaps", -1);
	public Input<Boolean> swapNeighboursOnlyInput = new Input<>("swapNeighboursOnly","if true, only neighbouring chains are swapped", false);

	
	
	// nr of samples between re-arranging states
	int resampleEvery;	
	double temperatureScaler;
	int adaptTemperatureEvery;
	
	
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
		temperatureScaler = temperatureScalerInput.get();
		adaptTemperatureEvery = adaptTemperatureEveryInput.get();
		
	} // initAndValidate
	
	private void initRun(){
		String sXML = new XMLProducer().toXML(this);
		sXML = sXML.replaceAll("chains=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("resampleEvery=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("tempDir=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("temperatureScaler=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("logHeatedChains=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("adaptTemperatureEvery=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("swapNeighboursOnly=['\"][^ ]*['\"]", "");

	
        String sMCMCMC = this.getClass().getName();
		while (sMCMCMC.length() > 0) {
			sXML = sXML.replaceAll("\\b" + CoupledMCMC.class.getName() + "\\b", HeatedChain.class.getName());
			if (sMCMCMC.indexOf('.') >= 0) {
				sMCMCMC = sMCMCMC.substring(sMCMCMC.indexOf('.') + 1);
			} else {
				sMCMCMC = "";
			}
		}
		long nSeed = Randomizer.getSeed();
			

		
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
				else if (i != 0 && !logHeatedChainsInput.get()) {
					for (int j = 0; j < chains[i].loggersInput.get().size(); j++){
						if (chains[i].loggersInput.get().get(j).getID().contentEquals("screenlog"))
							chains[i].loggersInput.get().remove(j);
					}
				}			
				
				// initialize each chain individually
				chains[i].setChainNr(i, resampleEvery, temperatureScaler);
				chains[i].setStateFile(stateFileName.replace(".state", "." + i + "state"), restoreFromFile);				
				
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
				int i,j;
				// resample state
				if (swapNeighboursOnlyInput.get()){
					i = Randomizer.nextInt(chains.length-1);
					j = i+1;					
				}else{
					i = Randomizer.nextInt(chains.length);
					j = i;
					while (i == j) {
						j = Randomizer.nextInt(chains.length);
					}
					if (i > j) {
						int tmp = i; i = j; j = tmp;
					}
				}
				
				
				double p1before = chains[i].getCurrentLogLikelihood();
				double p2before = chains[j].getCurrentLogLikelihood();

				// robust calculations can be extremly expensive, just calculate the new probs instead 
				double p1after = chains[i].getCurrentLogLikelihood() * chains[i].getTemperature() / chains[j].getTemperature();
				double p2after = chains[j].getCurrentLogLikelihood() * chains[j].getTemperature() / chains[i].getTemperature();
				
								
				double logAlpha = p1after - p1before + p2after - p2before;
				System.err.println(successfullSwaps0 + " " + successfullSwaps + ": " + i + " <--> " + j + ": " + logAlpha +  ": " + ((double) successfullSwaps/totalSwaps) + ": " + temperatureScaler);
				if (Math.exp(logAlpha) > Randomizer.nextDouble()) {
					if (totalSwaps>99)
						successfullSwaps++;
					
					if (i == 0) {
						successfullSwaps0++;
					}
					System.err.println(i + " <--> " + j);
					swapStates(chains[i], chains[j]);
					chains[i].calcCurrentLogLikelihoodRobustly();
					chains[j].calcCurrentLogLikelihoodRobustly();

					//assignState(chains[j], chains[i]);
					//chains[j].calcCurrentLogLikelihoodRobustly();
				}
				totalSwaps++;
				
				if (adaptTemperatureEvery!=-1 && 
						totalSwaps % adaptTemperatureEvery==0 &&
						totalSwaps > 110){
					// adapt the temperature scaler
					updateTemperature(totalSwaps, (double) successfullSwaps/(totalSwaps-100));
					for (int k = 1; k < chains.length; k++) {
						chains[k].setTemperature(k, temperatureScaler);
					}
					 
					
				}
					
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
	
    public double updateTemperature(int totalSwaps, double successfulFraction) {
    	
//        // do no optimisation for the first N optimisable operations
//        if (autoOptimizeDelayCount < autoOptimizeDelay || !autoOptimise) {
//            autoOptimizeDelayCount++;
//            return 0;
//        }
        
        final double target = 0.1;

//        double count = (totalSwaps + 1.0);
//        switch (transform) {
//            case log:
//                double count = Math.log(totalSwaps + 1.0);
//                break;
//            case sqrt:
                double count = Math.sqrt(totalSwaps);
//                break;
//            case none:
//            	break;
//            default:
//            	break;
//        }

        double deltaP = successfulFraction - target;
  
        temperatureScaler *= Math.exp(1/count * deltaP);
//        deltaP += Math.log(1.0 / temperatureScaler - 1.0);        
        
//        System.out.println(totalSwaps + " " + logAlpha + " " + target+ " " + deltaP);
//        
//        temperatureScaler = 1.0 / (Math.exp(deltaP) + 1.0);

        if (deltaP > -Double.MAX_VALUE && deltaP < Double.MAX_VALUE) {
            return deltaP;
        }
        return 0;
    }

	
	private void assignState(HeatedChain mcmc1, HeatedChain mcmc2) {
		State state1 = mcmc1.startStateInput.get();
		State state2 = mcmc2.startStateInput.get();
		List<StateNode> stateNodes1 = state1.stateNodeInput.get();
		List<StateNode> stateNodes2 = state2.stateNodeInput.get();
		for (int i = 0; i < stateNodes1.size(); i++) {
			StateNode stateNode1 = stateNodes1.get(i);
			StateNode stateNode2 = stateNodes2.get(i);
			stateNode1.assignFromWithoutID(stateNode2);
		}
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



