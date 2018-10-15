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
	public Input<Integer> resampleEveryInput = new Input<Integer>("resampleEvery", "number of samples in between resampling (and possibly swappping) states", 1000);
	public Input<String> heatedMCMCClassInput = new Input<String>("heatedMCMCClass", "Name of the class used for heated chains", HeatedChain.class.getName());
	public Input<String> tempDirInput = new Input<>("tempDir","directory where temporary files are written","/tmp/");
	public Input<Double> temperatureScalerInput = new Input<>("temperatureScaler","temperature scaler, the higher this value, the hotter the chains",0.1);
	public Input<String> stateFileNameInput = new Input<>("stateFileName","name of the state file", "state.backup.xml");

	
	
	// nr of samples between re-arranging states
	int resampleEvery = 10;
	
	double temperatureScaler = 0.01;
	
	
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
		chains = new HeatedChain[nrOfChainsInput.get()];
		
		resampleEvery = resampleEveryInput.get();
		
		temperatureScaler = temperatureScalerInput.get();

		// the difference between the various chains is
		// 1. it runs an MCMC, not a  MultiplMCMC
		// 2. remove chains attribute
		// 3. output logs change for every chain
		// 4. log to stdout is removed to prevent clutter on stdout
		
		String sXML = new XMLProducer().toXML(this);
		sXML = sXML.replaceAll("chains=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("heatedMCMCClass=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("resampleEvery=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("tempDir=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("temperatureScaler=['\"][^ ]*['\"]", "");
		sXML = sXML.replaceAll("stateFileName=['\"][^ ]*['\"]", "");

	
        String sMCMCMC = this.getClass().getName();
		while (sMCMCMC.length() > 0) {
			sXML = sXML.replaceAll("\\b"+ HeatedChain.class.getName()+ "\\b", heatedMCMCClassInput.get());
			if (sMCMCMC.indexOf('.') >= 0) {
				sMCMCMC = sMCMCMC.substring(sMCMCMC.indexOf('.')+1);
			} else {
				sMCMCMC = "";
			}
		}
		long nSeed = Randomizer.getSeed();
		
		

		
		
		// create new chains		
		for (int i = 0; i < chains.length; i++) {
			XMLParser parser = new XMLParser();
			String sXML2 = sXML;
			sXML2 = sXML2.replaceAll("\\$\\(seed\\)", nSeed+i+"");
			
			try {
		        FileWriter outfile = new FileWriter(new File(tempDirInput.get() + "MCMCMC.xml"));
		        outfile.write(sXML2);
		        outfile.close();
				
				chains[i] = (HeatedChain) parser.parseFragment(sXML2, true);
	
				// remove all loggers, except for main chain
				if (i != 0) {
					chains[i].loggersInput.get().clear();
				}
				String chainFileName = stateFileNameInput.get()  + "." + i + ".state";
				
				chains[i].setChainNr(i, resampleEvery, temperatureScalerInput.get());
				chains[i].setStateFile(chainFileName, restoreFromFile);
			
				chains[i].run();
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}
	
		// reopen log files for main chain, which were closed at the end of run(); 
		//Logger.FILE_MODE = LogFileMode.resume;
		//for (Logger logger : m_chains[0].loggersInput.get()) {
		//	logger.init();
		//}
		
		// get a copy of the list of state nodes to facilitate swapping states
		tmpStateNodes = startStateInput.get().stateNodeInput.get();

		chainLength = chainLengthInput.get();
		finishTimes = new long[chains.length];
	} // initAndValidate
	
	
	
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
		
		int totalSwaps = 0;
		int successfullSwaps = 0, successfullSwaps0 = 0;
//		for (HeatedMCMC chain:m_chains) {
//			chain.startStateInput.get().setEverythingDirty(true);
//		}
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
//			for (int dummy = 1; dummy < chains.length; dummy++){
				// resample state
				int i = Randomizer.nextInt(chains.length);
				int j = i;
				while (i == j) {
					j = Randomizer.nextInt(chains.length);
				}
				if (i > j) {
					int tmp = i; i = j; j = tmp;
				}
				
				
				double p1before = chains[i].getCurrentLogLikelihood();
				double p2before = chains[j].getCurrentLogLikelihood();
//				swapStates(chains[i], chains[j]);
//				double p1after = chains[i].calcCurrentLogLikelihoodRobustly();
//				double p2after = chains[j].calcCurrentLogLikelihoodRobustly();
				// robust calculations can be extremly expensive, just calculate the new probs instead 
				double p1after = chains[i].getCurrentLogLikelihood() * chains[i].getTemperature() / chains[j].getTemperature();
				double p2after = chains[j].getCurrentLogLikelihood() * chains[j].getTemperature() / chains[i].getTemperature();

				
//				System.out.println(p1after + " " + chains[i].getCurrentLogLikelihood() + " " + p2after + " " + chains[j].getCurrentLogLikelihood());
//				System.out.println(chains[i].getTemperature() + " " + chains[j].getTemperature()) ;
				
								
				double logAlpha = p1after - p1before + p2after - p2before;
				System.err.println(successfullSwaps0 + " " + successfullSwaps + ": " + i + " <--> " + j + ": " + logAlpha +  ": " + ((double) successfullSwaps/totalSwaps));
				if (Math.exp(logAlpha) < Randomizer.nextDouble()) {
					// swap fails
					//assignState(chains[i], chains[j]);
//					swapStates(chains[i], chains[j]);
					// not neeeded since the temperature scaler is only applied later
//					chains[i].calcCurrentLogLikelihoodRobustly();
//					chains[j].calcCurrentLogLikelihoodRobustly();

					
				} else {
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
				
				// tuning
				for (int k = 1; k < chains.length; k++) {
					chains[k].optimiseRunTime(startTime, finishTimes[k], finishTimes[0]);
				}
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



