package beast.coupledMCMC;



import beast.app.beauti.BeautiDoc;
import beast.core.*;
import beast.core.util.Log;
import beast.util.Randomizer;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.UnsupportedEncodingException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BetaDistributionImpl;
import org.json.JSONException;
import org.json.JSONObject;
import org.xml.sax.SAXException;

		
@Description("Serial temperaing aka umbrella sampling")
@Citation("Geyer, C.J., 2011. Importance sampling, simulated tempering, and umbrella sampling. In Handbook of Markov Chain Monte Carlo (pp. 295-311). Chapman & Hall/CRC, Boca Raton.")
public class SerialMCMC extends MCMC {
	public Input<Integer> nrOfChainsInput = new Input<Integer>("chains", "number of chains to run in parallel (default 2)", 2);
	public Input<Integer> resampleEveryInput = new Input<Integer>("resampleEvery", "number of samples in between resampling (and possibly swappping) states", 100);
	
	// input of the difference between temperature scalers
	public Input<Double> deltaTemperatureInput = new Input<>("deltaTemperature","temperature difference between the i-th and the i-th+1 chain", 0.01);
	
	// Input on whether the temperature between chains should be optimized
	public Input<Boolean> optimiseInput = new Input<>("optimise", "if true, the chain constants are automatically optimised to reach a target acceptance probability", true);
	
	public Input<Integer> optimiseDelayInput = new Input<>("optimiseDelay", "after this many epochs/swaps, the temperature will be optimized (if optimising is set to true)", 100);
	
	public Input<Boolean> useBetaDistributionInput = new Input<>("useBetaDistribution", "if true, the spacing between chains is assumed to be for the quantiles of a beta distribution with alpha=1 and beta=tuneable", false);

	public Input<Integer> logCounterInput = new Input<>("logCounter", "subsample frequency of screen log output: higher means less screen output", 100);

	
	// nr of samples between re-arranging states
	int resampleEvery;	
	
	double deltaTemperature;
	
	private boolean optimise;
	
	/** plugins representing MCMC with model, loggers, etc **/
	HeatedChain [] chains;
	
	/** threads for running MCMC chains **/
	// Thread [] threads;
	
	/** beta distribution for spacing*/
	org.apache.commons.math.distribution.BetaDistributionImpl m_dist;
	
	/** keep track of time taken between logs to estimate speed **/
    long startLogTime;
    
    /** keeps track of the last n swaps and if they were accepted or not **/
    // List<Boolean> acceptedSwaps;

	// keep track of when threads finish in order to optimise thread usage
	long [] finishTimes;
		
//	private ArrayList<Long>[] runTillIteration;
	
	int totalSwaps = 0;
	int successfullSwaps = 0, successfullSwaps0 = 0;
	
	long sampleOffset = 0;
	
	// defines which scheme for spacing between the heated chains to use
	private enum Spacing{Geometric, Beta };
	private Spacing spacing;
	
	/** normalisation constant, which helps compensate for difference in posteriors under different temperatures
	 * will be auto-optimised during run
	 **/ 
	private double [] c;
	
	private ScreenLogger screenlog;
	
	
	class ScreenLogger extends Logger {
		
		ScreenLogger(Logger other) {
			initByName("log", other.loggersInput.get(), "logEvery", other.everyInput.get());
		}
		
		@Override
		public void log(long sampleNr) {
	        if ((sampleNr < 0) || (sampleNr % every > 0)) {
	            return;
	        }
	        if (sampleOffset >= 0) {
	            if (sampleNr == 0) {
	                // don't need to duplicate the last line in the log
	                return;
	            }
	            sampleNr += sampleOffset;
	        }

	        ByteArrayOutputStream baos = new ByteArrayOutputStream();
	        PrintStream out = new PrintStream(baos);

	        for (final BEASTObject m_logger : loggersInput.get()) {
	            ((Loggable)m_logger).log(sampleNr, out);
	        }

	        // Acquire log string and trim excess tab
	        String logContent;
	        try {
	            logContent = baos.toString("ASCII").trim();
	        } catch (UnsupportedEncodingException e) {
	            throw new RuntimeException("ASCII string encoding not supported: required for logging!");
	        }

            // logContent = prettifyLogLine(logContent);
            m_out.print(logContent);
            m_out.print("\t");
	 	}
	}

	
	OperatorStats [] operatorStats;
	
	public class OperatorStats {
		String [] stats;
	    
	    OperatorStats(int nrOfOperators) {
	    	stats = new String[nrOfOperators];
	    }
	    
	    public void store(List<Operator> operators) {
	    	for (int i = 0; i < operators.size(); i++) {
	    		Operator operator = operators.get(i);
	    		
	    		StringWriter out    = new StringWriter();
	    	    PrintWriter  writer = new PrintWriter(out);
	    		operator.storeToFile(writer);
	        	stats[i] = out.toString();
	    	}
	    }
	    
	    public void restore(List<Operator> operators) {
	    	if (stats[0] == null) {
	    		return;
	    	}
			try {
	    	for (int i = 0; i < operators.size(); i++) {
		    		JSONObject o = new JSONObject(stats[i]);
		    		Operator operator = operators.get(i);
		    		operator.restoreFromFile(o);
	    	}
			} catch (JSONException e) {
				// should not get here if operators implement storeToFile() correctly
				e.printStackTrace();
			}
	    }
	    
	    public void storeToFile(PrintStream out) {
	    	for (String s : stats) {
	    		if (s != null) {
	    			s = s.replaceAll("\n", "\\n");
	    			out.println(s);
	    		}
	    	}
	    }
	    
	    public void restorFromFile(String [] strs, int start) {
	    	for (int i = 0; i < stats.length; i++) {
	    		stats[i] = strs[start + i].replaceAll("\\n", "\n");
	    	}
	    }
	} // class OperatorStats
	
	
	

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
		
		c = new double[chains.length];
		
		resampleEvery = resampleEveryInput.get();				
		
		deltaTemperature = deltaTemperatureInput.get();
		
		// acceptedSwaps = new ArrayList<>();				
		
		optimise = optimiseInput.get();
		
		if (useBetaDistributionInput.get()) {
			m_dist = new BetaDistributionImpl(1, deltaTemperature);
			spacing = Spacing.Beta;
		}else {
			spacing = Spacing.Geometric;
		}
				
				
	} // initAndValidate
	
	private void initRun(){
		operatorStats = new OperatorStats[chains.length];
		int nrOfOperators = operatorsInput.get().size();
		for (int i = 0; i < chains.length; i++) {
			operatorStats[i] = new OperatorStats(nrOfOperators);
		}
		
		
		
//		XMLProducer p = new XMLProducer();
//		String sXML = p.toXML(this, new ArrayList<>());
//		
//		// removes coupled MCMC parts of the xml		
//		sXML = sXML.replaceAll("chains=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("resampleEvery=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("tempDir=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("deltaTemperature=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("maxTemperature=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("optimise=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("logHeatedChains=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("optimizeDelay=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("optimizeEvery=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("nrExchanges=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("preSchedule=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("target=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("neighbourSwapping=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("useBetaDistribution=['\"][^ ]*['\"]", "");
//		sXML = sXML.replaceAll("spec=\"Logger\"", "");
//		sXML = sXML.replaceAll("<logger", "<coupledLogger spec=\"beast.coupledMCMC.CoupledLogger\"");
//		sXML = sXML.replaceAll("</logger", "</coupledLogger");
		
		// check if the loggers have a same issue
//        String sMCMCMC = this.getClass().getName();
//		while (sMCMCMC.length() > 0) {
//			sXML = sXML.replaceAll("\\b" + SerialMCMC.class.getName() + "\\b", HeatedChain.class.getName());
//			if (sMCMCMC.indexOf('.') >= 0) {
//				sMCMCMC = sMCMCMC.substring(sMCMCMC.indexOf('.') + 1);
//			} else {
//				sMCMCMC = "";
//			}
//		}
			
		chains[0] = new HeatedChain() {
			
			// adjusted so it won't restore from the state file
		    @Override
		    public void run() throws IOException, SAXException, ParserConfigurationException {
		        // set up state (again). Other plugins may have manipulated the
		        // StateNodes, e.g. set up bounds or dimensions
		        state.initAndValidate();
		        // also, initialise state with the file name to store and set-up whether to resume from file
		        state.setStateFileName(stateFileName);
		        operatorSchedule.setStateFileName(stateFileName);

		        burnIn = burnInInput.get();
		        chainLength = chainLengthInput.get();
		        int nInitialisationAttempts = 0;
		        state.setEverythingDirty(true);
		        posterior = posteriorInput.get();

		        if (restoreFromFile) {
//		            state.restoreFromFile();
//		            operatorSchedule.restoreFromFile();
		            burnIn = 0;
		            oldLogLikelihood = state.robustlyCalcPosterior(posterior);
		        } else {
		            do {
		                for (final StateNodeInitialiser initialiser : initialisersInput.get()) {
		                    initialiser.initStateNodes();
		                }
		                oldLogLikelihood = state.robustlyCalcPosterior(posterior);
		                nInitialisationAttempts += 1;
		            } while (Double.isInfinite(oldLogLikelihood) && nInitialisationAttempts < numInitializationAttempts.get());
		        }
		        // final long startTime = System.currentTimeMillis();
		        

		        // do the sampling
		        logAlpha = 0;
		        debugFlag = Boolean.valueOf(System.getProperty("beast.debug"));

		        System.err.println("Start likelihood: " + oldLogLikelihood + " " + (nInitialisationAttempts > 1 ? "after " + nInitialisationAttempts + " initialisation attempts" : ""));
		        if (Double.isInfinite(oldLogLikelihood) || Double.isNaN(oldLogLikelihood)) {
		            reportLogLikelihoods(posterior, "");
		            throw new RuntimeException("Could not find a proper state to initialise. Perhaps try another seed.");
		        }

		        loggers = coupledLoggersInput.get();

		        // put the loggers logging to stdout at the bottom of the logger list so that screen output is tidier.
		        Collections.sort(loggers, (o1, o2) -> {
		            if (o1.isLoggingToStdout()) {
		                return o2.isLoggingToStdout() ? 0 : 1;
		            } else {
		                return o2.isLoggingToStdout() ? -1 : 0;
		            }
		        });
		        // warn if none of the loggers is to stdout, so no feedback is given on screen
		        boolean hasStdOutLogger = false;
		        boolean hasScreenLog = false;
		        for (CoupledLogger l : loggers) {
		        	if (l.isLoggingToStdout()) {
		        		hasStdOutLogger = true;
		        	}
		        	if (l.getID() != null && l.getID().equals("screenlog")) {
		        		hasScreenLog = true;
		        	}
		        }
		        if (!hasStdOutLogger) {
		        	Log.warning.println("WARNING: If nothing seems to be happening on screen this is because none of the loggers give feedback to screen.");
		        	if (hasScreenLog) {
		        		Log.warning.println("WARNING: This happens when a filename  is specified for the 'screenlog' logger.");
		        		Log.warning.println("WARNING: To get feedback to screen, leave the filename for screenlog blank.");
		        		Log.warning.println("WARNING: Otherwise, the screenlog is saved into the specified file.");
		        	}
		        }

		        // initialises log such that log file headers are written, etc.
		        for (final Logger log : loggers) {
		            log.init();
		        }
		        
//		        Randomizer = new Randomizer();
		                
		    } // run;
		  
			
			// run MCMC inner loop for resampleEvery nr of samples
			// suppress storage of state
			protected long runTillResample(long runUntil) throws Exception {
				int corrections = 0;
				final boolean isStochastic = posterior.isStochastic();
				

				for (long sampleNr = currentSample; sampleNr <= runUntil; sampleNr++) {
		            final Operator operator = propagateState(sampleNr);
		            
		            if (debugFlag && sampleNr % 3 == 0 || sampleNr % 10000 == 0) {
		                // check that the posterior is correctly calculated at every third
		                // sample, as long as we are in debug mode
		            	final double originalLogP = isStochastic ? posterior.getNonStochasticLogP() : oldLogLikelihood;
		                final double logLikelihood = isStochastic ? state.robustlyCalcNonStochasticPosterior(posterior) : state.robustlyCalcPosterior(posterior);
		                if (isTooDifferent(logLikelihood, originalLogP)) {
		                    reportLogLikelihoods(posterior, "");
		                    Log.err.println("At sample " + sampleNr + "\nPosterior incorrectly calculated: " + originalLogP + " != " + logLikelihood
		                    		+ "(" + (originalLogP - logLikelihood) + ")"
		                            + " Operator: " + operator.getClass().getName());
//		                    System.exit(0);
		                }
		                if (sampleNr > NR_OF_DEBUG_SAMPLES * 3) {
		                    // switch off debug mode once a sufficient large sample is checked
		                    debugFlag = false;
		                    if (isTooDifferent(logLikelihood, originalLogP)) {
		                        // incorrect calculation outside debug period.
		                        // This happens infrequently enough that it should repair itself after a robust posterior calculation
		                        corrections++;
		                        if (corrections > 100) {
		                            // after 100 repairs, there must be something seriously wrong with the implementation
		                        	Log.err.println("Too many corrections. There is something seriously wrong that cannot be corrected");
		                            state.storeToFile(sampleNr);
		                            operatorSchedule.storeToFile();
		                            System.exit(1);
		                        }
		                        oldLogLikelihood = state.robustlyCalcPosterior(posterior);
		                    }
		                } else {
		                    if (isTooDifferent(logLikelihood, originalLogP)) {
		                        // halt due to incorrect posterior during intial debug period
		                        state.storeToFile(sampleNr);
		                        operatorSchedule.storeToFile();
		                        System.exit(1);
		                    }
		                }
		            } else {
		                if (sampleNr >= 0) {
		                	operator.optimize(logAlpha);
		                }
		            }
		            callUserFunction(sampleNr);
		            
		            if (posterior.getCurrentLogP() == Double.POSITIVE_INFINITY) {
		            	throw new RuntimeException("Encountered a positive infinite posterior. This is a sign there may be numeric instability in the model.");
		            }
		        }
		        if (corrections > 0) {
		        	Log.err.println("\n\nNB: " + corrections + " posterior calculation corrections were required. This analysis may not be valid!\n\n");
		        }

		        // the plus 1 is required such that the currentSample is not recomputed
		        currentSample = runUntil+1;
		        
		        return System.currentTimeMillis();
			}
		};
		
		screenlog = null;

		List<CoupledLogger> coupledLoggers = new ArrayList<>();
		for (Logger logger : loggersInput.get()) {
			if (logger.fileNameInput.get() != null) {
				CoupledLogger coupledLogger = new CoupledLogger();
				coupledLogger.initByName("logEvery", logger.everyInput.get(),
						"log", logger.loggersInput.get(),
						"fileName", logger.fileNameInput.get(),
						"mode", logger.modeInput.get()
						);
				coupledLogger.setID(logger.getID());
				coupledLoggers.add(coupledLogger);
			} else {
				screenlog = new ScreenLogger(logger);
			}
		}
		
		chains[0].initByName("operator", operatorsInput.get(), 
				"chainLength", chainLengthInput.get(),
				"state", startStateInput.get(),
				"init", initialisersInput.get(),
				"storeEvery", storeEveryInput.get(),
				"numInitializationAttempts", numInitializationAttempts.get(),
				"distribution", posteriorInput.get(),
				"coupledLogger", coupledLoggers,
				"sampleFromPrior", sampleFromPriorInput.get(),
				"operatorschedule", operatorScheduleInput.get()
				);
		
		int i = 0;
		// create new chains		
//		for (int i = 0; i < chains.length; i++) {
//			XMLParser parser = new XMLParser();
//			String sXML2 = sXML;
//			if (i>0){
//				sXML2 = sXML2.replaceAll("fileName=\"", "fileName=\"chain" + i+ "");
//			}
//			
			try {
//				chains[i] = (HeatedChain) parser.parseFragment(sXML2, true);
	
				// remove all screen loggers
//				for (int j = chains[i].coupledLoggersInput.get().size()-1; j >=0 ; j--){
//					if (chains[i].coupledLoggersInput.get().get(j).getID().contentEquals("screenlog")){
//						if (i == 0) {
//							screenlog = new ScreenLogger(chains[i].coupledLoggersInput.get().get(j)); 
//						}
//						chains[i].coupledLoggersInput.get().remove(j);
//					}
//				}
//				
//				// remove all loggers of heated chains if they are not logged
//				if (i != 0){
//					for (int j = 0; j < chains[i].coupledLoggersInput.get().size(); j++)
//						chains[i].coupledLoggersInput.get().get(j).setSuppressLogging(true);					
//				}
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
//		}
		for (i = 1; i < chains.length; i++) {
			chains[i] = new HeatedChain();
			
			chains[i].setResampleEvery(resampleEvery);
			chains[i].setTemperature(i, getTemperature(i));
			chains[i].setChainNr(i);
		}
		// ensure that each chain has the same starting point
		// threads = new Thread[chains.length];
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
	
	@Override 
	public void run() throws IOException {
		initRun();
		
		if (restoreFromFile) {
			try {
	            startStateInput.get().restoreFromFile();

				String str = BeautiDoc.load(new File(stateFileName));
				str = str.substring(str.indexOf("<!--") + 5);
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
					// else if (Spltstr[0].contentEquals("lastSwaps")) {
					//	String[] tmp = Spltstr[1].replaceAll("\\[", "").replaceAll("\\]", "").replaceAll("\\s", "").split(",");
						// acceptedSwaps = new ArrayList<>();
						//for (int j = 0; j < tmp.length; j++)
							// acceptedSwaps.add(Boolean.parseBoolean(tmp[j]));
					//}					
				}
				int delta = operatorsInput.get().size();
				for (int i = 0; i < chains.length; i++) {
					operatorStats[i].restorFromFile(strs, delta * i + 5 + chains.length);
				}
				Log.warning("Restoring: totalSwaps=" + totalSwaps + 
						" successfullSwaps=" + successfullSwaps +
						" successfullSwaps0=" + successfullSwaps0 +
						" deltaTemperature=" + deltaTemperature);
			} catch (SAXException|IOException|ParserConfigurationException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		
		runNeigbours();
	} // run
	
	private void runNeigbours(){
		for (int k = 0; k < 1 && k < chains.length; k++) {
			chains[k].currentSample = 0;
			try {
				chains[k].runTillResample(resampleEvery);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		startLogTime = -1;
		long startSample = 0;
		int [] chainCount = new int[chains.length];
		
		// print header for system output
		System.out.println("sample\tchain nr\tswapsColdCain\tswapProbability\tdeltaTemperature");		
		double currProb = ((double) successfullSwaps/totalSwaps);
		if (Double.isNaN(currProb))
			currProb = 0.0;
		try {
			screenlog.init();
			System.out.print("\t" + (startSample + sampleOffset) + "\t" + Arrays.toString(chainCount) + "\t" + successfullSwaps0 + "\t" + currProb + "\t" + deltaTemperature);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}


		int logCounter = 0;
		for (long sampleNr = resampleEvery; sampleNr <= chainLength; sampleNr += resampleEvery) {	
			// get the chains to swap, conditioning on them being neighbours
			int chain_i = chains[0].getChainNr(); //Randomizer.nextInt(chains.length-1);
			
			int neighbour = Randomizer.nextBoolean() ?
					(chain_i == 0 ? 1 : chain_i - 1) :
					(chain_i == chains.length - 1 ? chain_i - 1 : chain_i + 1);
			
			int i = 0, j = -1;
			
			for (int k = 0; k<chains.length; k++) {
				if (chains[k].getChainNr() == neighbour)
					j = k;	
			}
			
			// robust calculations can be extremely expensive, just calculate the new probs instead 
			double p1after = chains[i].getUnscaledCurrentLogLikelihood() * chains[j].getBeta();
			double p2after = chains[i].getUnscaledCurrentLogLikelihood() * chains[i].getBeta();
			
			double logAlpha = chain_i == 0 || chain_i == chains.length-1 ? -Math.log(0.5) : 0;
			logAlpha += neighbour == 0 || neighbour == chains.length-1 ? Math.log(0.5) : 0;
			logAlpha += (p1after - p2after); 
			logAlpha += c[chain_i] - c[neighbour];
						
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
				//swapLoggers(chains[i], chains[j]);				
				
				// swap Operator tuning
				//swapOperatorTuning(chains[i], chains[j]);
				operatorStats[chain_i].store(operatorsInput.get());
				operatorStats[neighbour].restore(operatorsInput.get());
								
			}
			totalSwaps++;
			
			if (optimise) {
				if (totalSwaps > optimiseDelayInput.get()) {
					// automatic weight estimation by using the mean of
					// observed logP for each temperature
					// Keeps on updating unlike the scheme in:
					// Park, S. and Pande, V.S., 2007. Choosing weights for simulated tempering. Physical Review E, 76(1), p.016703.
					// DOI: 10.1103/PhysRevE.76.016703
					double logP = chains[0].getUnscaledCurrentLogLikelihood();
					if (c[0] == 0) {
						for (int k = 0; k < chains.length; k++) {
							c[chains[k].getChainNr()] = logP * chains[k].getBeta(); 
						}
					} else {
						double delta = 1.0 / totalSwaps;
						for (int k = 0; k < chains.length; k++) {
							c[chains[k].getChainNr()] = c[chains[k].getChainNr()] * (1-delta) + logP * chains[k].getBeta() * delta; 
						}
					}
				}	            
			}

			if (sampleNr < chainLength) {
				chains[0].currentSample = sampleNr + 1;
				try {
					finishTimes[0] = chains[0].runTillResample(sampleNr+resampleEvery);
				} catch (Exception e) {
					e.printStackTrace();
				}
				if (chains[0].getChainNr() != 0) {
					sampleNr -= resampleEvery;
				}
			}

			
			chainCount[chains[0].getChainNr()]++;
			
			if (logCounter == logCounterInput.get()) {
				System.out.print("\t" + (sampleNr + sampleOffset) + "\t" + Arrays.toString(chainCount) + "\t" +
//						Arrays.toString(c) + "\t" + 
						successfullSwaps0 + "\t" + ((double) successfullSwaps/totalSwaps) + "\t" + deltaTemperature + " ");
				screenlog.log(0);
				logCounter = 0;
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
			} else {
				logCounter++;
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
		
		System.err.println("#Swap attemps = " + totalSwaps);
		System.err.println("#Successfull swaps = " + successfullSwaps);
		System.err.println("#Successfull swaps with cold chain = " + successfullSwaps0);
	}
	

//	/* swaps the states of mcmc1 and mcmc2 */
//	void swapLoggers(HeatedChain mcmc1, HeatedChain mcmc2) {
//		int mcmc2size = mcmc2.coupledLoggersInput.get().size();
//		int mcmc1size = mcmc1.coupledLoggersInput.get().size();
//		
//		for (int i = 0; i < mcmc2size; i++){
//			for (int j = 0; j < mcmc1size; j++){
//				if (mcmc2.coupledLoggersInput.get().get(i).getID().contentEquals(mcmc1.coupledLoggersInput.get().get(j).getID())){
//					boolean suppressLogging1 = mcmc1.coupledLoggersInput.get().get(i).getSuppressLogging();
//					boolean suppressLogging2 = mcmc2.coupledLoggersInput.get().get(j).getSuppressLogging();
//					
//					mcmc1.coupledLoggersInput.get().get(i).setSuppressLogging(suppressLogging2);
//					mcmc2.coupledLoggersInput.get().get(j).setSuppressLogging(suppressLogging1);
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
//
//	void swapOperatorTuning(HeatedChain mcmc1, HeatedChain mcmc2) {
//		List<Operator> operatorList1 = new ArrayList<Operator>();
//		List<Operator> operatorList2 = new ArrayList<Operator>();		
//		
//		operatorList1 = mcmc1.getOperatorSchedule().operatorsInput.get();
//		operatorList2 = mcmc2.getOperatorSchedule().operatorsInput.get();
//		
//		
//		for (int i = 0; i < operatorList1.size(); i++){
//			Operator operator1 = operatorList1.get(i);
//			Operator operator2 = operatorList2.get(i);
//
//		    int m_nNrRejected = operator1.get_m_nNrRejected();
//		    int m_nNrAccepted = operator1.get_m_nNrAccepted();
//		    int m_nNrRejectedForCorrection = operator1.get_m_nNrRejectedForCorrection();
//		    int m_nNrAcceptedForCorrection = operator1.get_m_nNrAcceptedForCorrection();
//		    
//		    operator1.setAcceptedRejected(operator2.get_m_nNrRejected(), operator2.get_m_nNrAccepted(), operator2.get_m_nNrRejectedForCorrection(), operator2.get_m_nNrAcceptedForCorrection());
//		    operator2.setAcceptedRejected(m_nNrAccepted, m_nNrRejected, m_nNrAcceptedForCorrection, m_nNrRejectedForCorrection);
//		    
//		    
//		    double coercableParameterValue = operator1.getCoercableParameterValue();
//		    operator1.setCoercableParameterValue(operator2.getCoercableParameterValue());
//		    operator2.setCoercableParameterValue(coercableParameterValue);
//			
//		}		
//	}
	
	
	void storeStateToFile(long currentSample) throws FileNotFoundException {
        PrintStream out;
		try {
			out = new PrintStream(stateFileName + ".new");
			String state = startStateInput.get().toXML(0);
			out.print(state);
			out.print("<!--\n");
			out.print("sample=" + (currentSample + sampleOffset) + "\n");
	        out.print("totalSwaps=" + totalSwaps + "\n");
	        out.print("successfullSwaps=" + successfullSwaps + "\n");
	        out.print("successfullSwaps0=" + successfullSwaps0 + "\n");
	        out.print("deltaTemperature=" + deltaTemperature + "\n");
	
	        for (int k = 0; k < chains.length; k++) {
	        	out.print("beta_chain." + k + "=" + chains[k].getBeta() + "\n");
	        }
	        for (int k = 0; k < chains.length; k++) {
	        	operatorStats[k].storeToFile(out);
	        }
			out.print("-->\n");
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
	
}



