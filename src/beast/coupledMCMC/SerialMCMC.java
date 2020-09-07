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

		
@Description("Serial MCMC aka serial tempering aka simulated tempering aka umbrella sampling")
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
    
	// keep track of when threads finish in order to optimise thread usage
//	long [] finishTimes;
	
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
	private int [] chainCount;
	private ScreenLogger screenlog;
	
	
	class ScreenLogger extends Logger {
		
		ScreenLogger(Logger other) {
			initByName("log", other.loggersInput.get(), "logEvery", other.everyInput.get());
		}
		
		@Override
		public void log(long sampleNr) {
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
	    		} else {
	    			out.println();
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
		} else {
			spacing = Spacing.Geometric;
		}
				
		chainCount = new int[chains.length];

	} // initAndValidate
	
	private void initRun(){
		operatorStats = new OperatorStats[chains.length];
		int nrOfOperators = operatorsInput.get().size();
		for (int i = 0; i < chains.length; i++) {
			operatorStats[i] = new OperatorStats(nrOfOperators);
		}
			
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
		            state.restoreFromFile();
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

		        // initialises log such that log file headers are written, etc.
		        for (final Logger log : loggers) {
		            log.init();
		        }
		        		                
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
						"mode", logger.modeInput.get(),
						"sanitiseHeaders", logger.sanitiseHeadersInput.get()
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
		try {
			// initialize each chain individually
			chains[i].setResampleEvery(resampleEvery);
			chains[i].setTemperature(i, getTemperature(i));
							
			chains[i].setStateFile(stateFileName, restoreFromFile);
			chains[i].setChainNr(i);
			chains[i].run();

		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		for (i = 1; i < chains.length; i++) {
			chains[i] = new HeatedChain();
			
			chains[i].setResampleEvery(resampleEvery);
			chains[i].setTemperature(i, getTemperature(i));
			chains[i].setChainNr(i);
		}
		

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
					else if (Spltstr[0].startsWith("beta_chain")) {
						int chainNr = Integer.parseInt(Spltstr[0].substring(10));
						chains[chainNr].setChainNr(chainNr);
						double beta = Double.parseDouble(Spltstr[1]);
						chains[chainNr].setBeta(beta);
						c[chainNr] = Double.parseDouble(Spltstr[2]);
						chainCount[chainNr] = Integer.parseInt(Spltstr[3]);
					}
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
		
		// print header for system output
		System.out.print("sample\tchain nr\tswapsColdCain\tswapProbability\tdeltaTemperature");		
		double currProb = ((double) successfullSwaps/totalSwaps);
		if (Double.isNaN(currProb))
			currProb = 0.0;
		try {
			screenlog.init();
//			System.out.print("\t" + (startSample + sampleOffset) + "\t" + Arrays.toString(chainCount) + "\t" + successfullSwaps0 + "\t" + currProb + "\t" + deltaTemperature);
		} catch (IOException e1) {
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
			
			double logAlpha = chain_i == 0 || chain_i == chains.length-1 ? -Math.log(2) : 0;
			logAlpha += neighbour == 0 || neighbour == chains.length-1 ? Math.log(2) : 0;
			logAlpha += (p1after - p2after); 
			logAlpha += c[chain_i] - c[neighbour];
			if (chainCount[chain_i] > 100 && chainCount[neighbour] > 100) {
				logAlpha += Math.log(chainCount[chain_i]) - Math.log(chainCount[neighbour]);
			}
						
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
				for (CoupledLogger logger : chains[0].loggers) {
					logger.setSuppressLogging(neighbour != 0);
				}
				
				// swap Operator tuning
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
					chains[0].runTillResample(sampleNr+resampleEvery);
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
	        	out.print("beta_chain " + chains[k].getChainNr() + "=" + chains[k].getBeta() + 
	        			"=" + c[chains[k].getChainNr()] + "=" + chainCount[chains[k].getChainNr()] + "\n");
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
			e.printStackTrace();
		}

	}
	
}



