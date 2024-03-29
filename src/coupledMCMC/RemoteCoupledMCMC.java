package coupledMCMC;


import beastfx.app.inputeditor.BeautiDoc;
import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Citation.Citations;
import beast.base.core.Input.Validate;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.core.Log;
import beast.base.util.HeapSort;
import beast.base.util.Randomizer;
import beast.base.parser.XMLProducer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

@Citations(
		{
		@Citation(value= "Müller, Nicola Felix, and Remco Bouckaert. Coupled MCMC in BEAST 2. bioRxiv (2019): 603514. ", 
				year = 2019, firstAuthorSurname = "Müller",
				DOI="10.1101/603514"),
		@Citation(value="Bouckaert, Remco, Timothy G. Vaughan, Joëlle Barido-Sottani, Sebastián Duchêne, \n"
				+ "  Mathieu Fourment, Alexandra Gavryushkina, Joseph Heled et al. \n"
				+ "  BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis. \n"
				+ "  PLoS computational biology 15, no. 4 (2019): e1006650.", 
		        year = 2019, firstAuthorSurname = "bouckaert",
				DOI="10.1371/journal.pcbi.1006650"),	
		@Citation(value= "Altekar G, Dwarkadas S, Huelsenbeck J and Ronquist F (2004). \n" +
				"  Parallel Metropolis Coupled Markov Chain Monte Carlo For Bayesian Phylogenetic Inference.\n" +
				"  Bioinformatics, 20(3), 407-415."
				, year = 2004, firstAuthorSurname = "Altekar",
				DOI="10.1093/bioinformatics/btg427")
		}
)		
@Description("Distributed Parallel Metropolis Coupled Markov Chain Monte Carlo")
public class RemoteCoupledMCMC extends MCMC {
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
	
	public Input<Boolean> preScheduleInput = new Input<>("preSchedule","if true, how long chains are run for is scheduled at the beginning", true);
	

	public Input<String> hostsInput = new Input<>("hosts","comma separated list of hosts -- should be equal to number of chains", Validate.REQUIRED);
	public Input<String> portsInput = new Input<>("ports","comma separated list of ports on hosts. If fewer ports than hosts are specified, it cycles through the list (so you can specify only 1 port, which is 5001 by default)", "5001");

	public Input<Integer> loggerportInput = new Input<>("loggerport", "port where the logger service listens -- should differ from ports", 5000);

	
	// nr of samples between re-arranging states
	int resampleEvery;	
	
	double deltaTemperature;
	
	private boolean optimise;
	
	/** plugins representing MCMC with model, loggers, etc **/
	RemoteHeatedChain [] chains;
	
	/** threads for running MCMC chains **/
	Thread [] threads;
	
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
	RemoteLoggerServer remoteLoggerServer;

	@Override
	public void initAndValidate() {
		if (nrOfChainsInput.get() < 1) {
			throw new RuntimeException("chains must be at least 1");
		}
		if (nrOfChainsInput.get() == 1) {
			Log.warning.println("Warning: coupled MCMC needs at least 2 chains, but the number of chains is 1. Running plain MCMC.");
		}
		// initialize the differently heated chains
		chains = new RemoteHeatedChain[nrOfChainsInput.get()];
		
		resampleEvery = resampleEveryInput.get();				
		
		deltaTemperature = deltaTemperatureInput.get();
		
		acceptedSwaps = new ArrayList<>();				
		
		optimise = optimiseInput.get();
		
		if (logHeatedChainsInput.get()) {
			throw new IllegalArgumentException("logHeatedChains should be set to false -- it only works with non-remote CoupledMCMC");
		}
		
		try {
			remoteLoggerServer = new RemoteLoggerServer(loggerportInput.get(), restoreFromFile);
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException("Cannot create logger service on port " + loggerportInput.get());
		}
		
//		// check how optimisation should be performed
//		if (optimiseInput.get()==null) {
//			optimise = Optimise.True;
//		}else {
//			String input = optimiseInput.get().replace("\\s+", "");
//			switch (input) {
//				case "false":
//					optimise = Optimise.False;	
//					break;
//				case "true":
//					optimise = Optimise.True;
//					break;
//				case "log":
//					optimise = Optimise.Log;
//					break;
//				case "sqrt":
//					optimise = Optimise.Sqrt;
//					break;
//				default:
//					throw new IllegalArgumentException("optimise input should either be \"false\", \"true\", \"log\", \"sqrt\" or not specified at all, which is qual to \"true\"");
//			}
//		}
		
	} // initAndValidate
	
	private void initRun(){
		XMLProducer p = new XMLProducer();
		//String sXML = p.toXML(this, new ArrayList<>());
		
		String loggerHost = "localhost";
		try {
			loggerHost = InetAddress.getLocalHost().getHostAddress();
		} catch (UnknownHostException e1) {
			e1.printStackTrace();
		}
		
		List<RemoteCoupledLogger> coupledLoggers = new ArrayList<>();
		for (Logger logger : loggersInput.get()) {
			if (logger.fileNameInput.get() != null) {
				RemoteCoupledLogger coupledLogger = new RemoteCoupledLogger(logger);
				coupledLoggers.add(coupledLogger);
			}			
		}
		
		HeatedChain heated = new HeatedChain();
		heated.initByName("distribution", posteriorInput.get(), 
				"operator", operatorsInput.get(),
				"state", startStateInput.get(),
				"init", initialisersInput.get(),
				"chainLength", chainLengthInput.get(),
				"storeEvery", storeEveryInput.get(),
				"numInitializationAttempts", numInitializationAttempts.get(),
				"coupledLogger", coupledLoggers
				);
		String sXML = p.toXML(heated, new ArrayList<>());

//		// removes coupled MCMC parts of the xml		
//		sXML = sXML.replaceFirst("hosts=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("ports=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("loggerport=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("chains=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("resampleEvery=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("tempDir=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("deltaTemperature=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("maxTemperature=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("optimise=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("logHeatedChains=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("optimizeDelay=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("optimizeEvery=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("nrExchanges=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("preSchedule=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceFirst("target=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("spec=\"Logger\"", "");
//		sXML = sXML.replaceAll("<logger", "<coupledLogger spec=\"beast.coupledMCMC.RemoteCoupledLogger\" host=\"" + loggerHost + "\" port=\"" + loggerportInput.get() + "\" ");
//		sXML = sXML.replaceAll("</logger", "</coupledLogger");
		
		// check if the loggers have a same issue
        String sMCMCMC = this.getClass().getName();
		while (sMCMCMC.length() > 0) {
			sXML = sXML.replaceAll("\\b" + RemoteCoupledMCMC.class.getName() + "\\b", HeatedChain.class.getName());
			if (sMCMCMC.indexOf('.') >= 0) {
				sMCMCMC = sMCMCMC.substring(sMCMCMC.indexOf('.') + 1);
			} else {
				sMCMCMC = "";
			}
		}
		
		String [] hosts = hostsInput.get().split(",");
		if (hosts.length != chains.length) {
			throw new IllegalArgumentException("Number of hosts should be equal to number of chains");
		}
		String [] ports = portsInput.get().split(",");
		
			
		// create new chains		
		for (int i = 0; i < chains.length; i++) {
//			XMLParser parser = new XMLParser();
			String sXML2 = sXML;
			if (i>0){
				sXML2 = sXML2.replaceAll("fileName=\"", "fileName=\"chain" + i+ "");
			}
			
			try {
				chains[i] = new RemoteHeatedChain(sXML2, hosts[i], Integer.parseInt(ports[i % ports.length]), restoreFromFile);
				
	
//				// remove all screen loggers
//				for (int j = chains[i].coupledLoggersInput.get().size()-1; j >=0 ; j--){
//					if (chains[i].coupledLoggersInput.get().get(j).getID().contentEquals("screenlog")){
//						chains[i].coupledLoggersInput.get().remove(j);
//					}
//				}
//				
//				// remove all loggers of heated chains if they are not logged
//				if (!logHeatedChainsInput.get() && i != 0){
//					for (int j = 0; j < chains[i].coupledLoggersInput.get().size(); j++)
//						chains[i].coupledLoggersInput.get().get(j).setSuppressLogging(true);					
//				}
				// initialize each chain individually
				chains[i].setTemperature(i, getTemperature(i));
				chains[i].setChainNr(i);
				// the following set resampleEvery but also triggers the remove MC3 to call run()
				chains[i].setResampleEvery(resampleEvery);
							
				// needed to avoid error of putting the working dir twice
				String[] splittedFileName = stateFileName.split("/");
				
				chains[i].setStateFile(splittedFileName[splittedFileName.length-1].replace(".state", "." + i + "state"), restoreFromFile);
//				chains[i].setChainNr(i);
//				chains[i].run();

			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}		
		// ensure that each chain has the same starting point
		threads = new Thread[chains.length];
		finishTimes = new long[chains.length];
		

		chainLength = chainLengthInput.get();
		
		// pre schedule which chains to swap when
		if (preScheduleInput.get()){
			buildSchedule();
		}
		
		if (restoreFromFile){
			System.out.println("restoring from file, printing to screen but not to loggers will start again from 0");
			System.out.println("we further assume that all chains ended in the same state, if logging heated chains" +
			" the different heated chains can have different amount of interations");
		}		
		
	}
	
	private double getTemperature(int i){
		return i*deltaTemperature;
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

		
		if (preScheduleInput.get()){
			runPrescheduled();
			return;
		}
	} // run
	
	private void runPrescheduled(){
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
		
		// print header for system output
		System.out.println("sample\tswapsColdCain\tswapProbability\tdeltaTemperature");		
		double currProb = ((double) successfullSwaps/totalSwaps);
		if (Double.isNaN(currProb))
			currProb = 0.0;
		System.out.print("\t" + (startSample + sampleOffset) + "\t" + successfullSwaps0 + "\t" + currProb + "\t" + deltaTemperature + "\n");

		
		for (long sampleNr = resampleEvery; sampleNr <= chainLength; sampleNr += resampleEvery) {	
			
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
			
			// ensure that i is the colder chain		
			if (chains[i].getChainNr() > chains[j].getChainNr()) {
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
				
				// swap chain numbers
				int chainNr = chains[i].getChainNr();
				chains[i].setChainNr(chains[j].getChainNr());
				chains[j].setChainNr(chainNr);
				
				// swap loggers and the state file names
				// swapLoggers(chains[i], chains[j]);				
				
				// swap Operator tuning
				swapOperatorTuning(chains[i], chains[j]);		
				
				if (chains[j].getBeta() < chains[i].getBeta())
					throw new IllegalArgumentException("error in temperatures of chains");
				
				if (optimise) {
					acceptedSwaps.add(true);
				}
			}else {
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
	            
	            
	            // set new temperatures asynchronously, in same order as before
	            for (int k = 0; k < n; k++) {
	            	chains[k].setTemperature(0, getTemperature(chains[k].getChainNr()));
	            }
	            
			}

			runTillIteration[i].remove(0);
			runTillIteration[j].remove(0);
			
			long runi = chainLength, runj = chainLength;
			
			if (runTillIteration[i].size()>0)  runi = runTillIteration[i].get(0);
			if (runTillIteration[j].size()>0)  runj = runTillIteration[j].get(0);
			
			threads[i] = new HeatedChainThread(i, runi);
			threads[i].start();
			
			threads[j] = new HeatedChainThread(j, runj);
			threads[j].start();
			
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
		for (int i = 0; i < threads.length; i++){
			try {
				threads[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
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
	
//	/* swaps the states of mcmc1 and mcmc2 */
//	void swapLoggers(RemoteHeatedChain mcmc1, RemoteHeatedChain mcmc2) {
//		int mcmc2size = mcmc2.coupledLoggersInput.get().size();
//		int mcmc1size = mcmc1.coupledLoggersInput.get().size();
//		
//		for (int i = 0; i < mcmc2size; i++){
//			for (int j = 0; j < mcmc1size; j++){
//				if (mcmc2.coupledLoggersInput.get().get(i).getID().contentEquals(mcmc1.coupledLoggersInput.get().get(j).getID())){
//					if (!logHeatedChainsInput.get()){
//						boolean suppressLogging1 = mcmc1.coupledLoggersInput.get().get(i).getSuppressLogging();
//						boolean suppressLogging2 = mcmc2.coupledLoggersInput.get().get(j).getSuppressLogging();
//						
//						mcmc1.coupledLoggersInput.get().get(i).setSuppressLogging(suppressLogging2);
//						mcmc2.coupledLoggersInput.get().get(j).setSuppressLogging(suppressLogging1);
//
//					}
//
//					PrintStream printstream2 = mcmc2.coupledLoggersInput.get().get(i).getM_out();
//					PrintStream printstream1 = mcmc1.coupledLoggersInput.get().get(j).getM_out();
//										
//					mcmc2.coupledLoggersInput.get().get(i).setM_out(printstream1);
//					mcmc1.coupledLoggersInput.get().get(j).setM_out(printstream2);
//				}					
//			}			
//		}		
//				
//		// swap the state file Names as well
//		String stateFileName1 = mcmc1.getStateFileName();
//		mcmc1.setStateFileName(mcmc2.getStateFileName());
//		mcmc2.setStateFileName(stateFileName1);		
//	}

	void swapOperatorTuning(RemoteHeatedChain mcmc1, RemoteHeatedChain mcmc2) {
		String tuning1 = mcmc1.getOperatorTuningInfo();
		String tuning2 = mcmc2.getOperatorTuningInfo();
		mcmc1.setOperatorTuningInfo(tuning2);
		mcmc2.setOperatorTuningInfo(tuning1);
	}
	
	/* makes the schedule of when to swap which chains */
	@SuppressWarnings("unchecked")
	void buildSchedule(){
		// initialize Arrays
		runTillIteration = new ArrayList[chains.length];
		for (int i = 0; i < chains.length; i++){
			runTillIteration[i] = new ArrayList<Long>();
		}

		for (long sampleNr = resampleEvery; sampleNr <= chainLength; sampleNr += resampleEvery) {
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



