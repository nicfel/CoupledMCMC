package beast.coupledMCMC;


import java.io.IOException;
import java.util.Collections;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.StateNodeInitialiser;
import beast.core.util.Evaluator;
import beast.core.util.Log;
import beast.util.Randomizer;

@Description("Base class for doing Metropolis coupled MCMC. Each instance represenst a chain at a different temperature.")
public class HeatedChain extends MCMC {
	
	// LAMBDA is temperature multiplier
	static double LAMBDA = 0.01;
	
	// temperature on which this chain runs
	protected double temperature = 1.0;
	
	// nr of samples between re-arranging states
	protected int resampleEvery = 1000;
	
	// keep track of total nr of states sampled, using currentSample
	protected int currentSample = 0;

	protected double getCurrentLogLikelihood() {
		return oldLogLikelihood / temperature;
	};

	// set chain number for a given lambda
	
	public void setChainNr(int i, int resampleEvery, double lambda) {
		LAMBDA = lambda; // unneeded step, but still here
		temperature = 1 + i  * LAMBDA;
		this.resampleEvery = resampleEvery;
	}
	
	public double getTemperature(){
		return temperature;
	}

	// sets all nodes dirty and recalculates everything
	protected double calcCurrentLogLikelihoodRobustly() {
		oldLogLikelihood = robustlyCalcPosterior(posterior);
		return getCurrentLogLikelihood();
	};
	
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
        //state.setEverythingDirty(true);
        posterior = posteriorInput.get();

        if (restoreFromFile) {
            state.restoreFromFile();
            operatorSchedule.restoreFromFile();
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
        final long startTime = System.currentTimeMillis();

        // do the sampling
        logAlpha = 0;
        debugFlag = Boolean.valueOf(System.getProperty("beast.debug"));

//        System.err.println("Start state:");
//        System.err.println(state.toString());

        System.err.println("Start likelihood: " + oldLogLikelihood + " " + (nInitialisationAttempts > 1 ? "after " + nInitialisationAttempts + " initialisation attempts" : ""));
        if (Double.isInfinite(oldLogLikelihood) || Double.isNaN(oldLogLikelihood)) {
            reportLogLikelihoods(posterior, "");
            throw new RuntimeException("Could not find a proper state to initialise. Perhaps try another seed.");
        }

        loggers = loggersInput.get();

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
        for (Logger l : loggers) {
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

        // initialises log so that log file headers are written, etc.
        for (final Logger log : loggers) {
            log.init();
        }

//        doLoop();
//
//        System.out.println();
//        operatorSchedule.showOperatorRates(System.out);
//
//        System.out.println();
//        final long endTime = System.currentTimeMillis();
//        System.out.println("Total calculation time: " + (endTime - startTime) / 1000.0 + " seconds");
//        close();
//
//        System.err.println("End likelihood: " + oldLogLikelihood);
////        System.err.println(state);
//        state.storeToFile(chainLength);
//        operatorSchedule.storeToFile();
        //Randomizer.storeToFile(stateFileName);
    } // run;
	
	
	// run MCMC inner loop for resampleEvery nr of samples
	protected long runTillResample() throws Exception {
	       int corrections = 0;
	       for (int sampleNr = currentSample; sampleNr < currentSample + resampleEvery; sampleNr++) {
	            final int currentState = sampleNr;

	            state.store(currentState);
//	            if (m_nStoreEvery > 0 && iSample % m_nStoreEvery == 0 && iSample > 0) {
//	                state.storeToFile(iSample);
//	            	operatorSchedule.storeToFile();
//	            }

	            final Operator operator = operatorSchedule.selectOperator();
	            //System.out.print("\n" + sampleNr + " " + operator.getName()+ ":");

	            final Distribution evaluatorDistribution = operator.getEvaluatorDistribution();
	            Evaluator evaluator = null;

	            if (evaluatorDistribution != null) {
	                evaluator = new Evaluator() {
	                    @Override
	                    public double evaluate() {
	                        double logP = 0.0;

	                        state.storeCalculationNodes();
	                        state.checkCalculationNodesDirtiness();

	                        try {
	                            logP = evaluatorDistribution.calculateLogP();
	                        } catch (Exception e) {
	                            e.printStackTrace();
	                            System.exit(1);
	                        }

	                        state.restore();
	                        state.store(currentState);

	                        return logP;
	                    }
	                };
	            }

	            final double logHastingsRatio = operator.proposal(evaluator);

	            if (logHastingsRatio != Double.NEGATIVE_INFINITY) {

	            	if (operator.requiresStateInitialisation()) {
	            		state.storeCalculationNodes();
	            		state.checkCalculationNodesDirtiness();
	            	}

	                newLogLikelihood = posterior.calculateLogP();

	                logAlpha = newLogLikelihood/temperature - oldLogLikelihood/temperature + logHastingsRatio; 
	                //System.out.println(logAlpha + " " + newLogLikelihood + " " + oldLogLikelihood);
	                if (logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
	                    // accept
	                    oldLogLikelihood = newLogLikelihood;
	                    state.acceptCalculationNodes();

	                    if (sampleNr >= 0) {
	                        operator.accept();
	                    }
	                    //System.out.print(" accept");
	                } else {
	                    // reject
	                    if (sampleNr >= 0) {
	                        operator.reject(newLogLikelihood == Double.NEGATIVE_INFINITY ? -1 : 0);
	                    }
	                    state.restore();
	                    state.restoreCalculationNodes();
	                    //System.out.print(" reject");
	                }
	                state.setEverythingDirty(false);
	            } else {
	                // operation failed
	                if (sampleNr >= 0) {
	                    operator.reject(-2);
	                }
	                state.restore();
					if (!operator.requiresStateInitialisation()) {
	                    state.setEverythingDirty(false);
	                    state.restoreCalculationNodes();
					}
	                //System.out.print(" direct reject");
	            }
	            log(sampleNr);

	            if (debugFlag && sampleNr % 3 == 0 || sampleNr % 10000 == 0) {
	                // check that the posterior is correctly calculated at every third
	                // sample, as long as we are in debug mode
	            	final double fNonStochasticLogP = posterior.getNonStochasticLogP();
	                final double fLogLikelihood = state.robustlyCalcNonStochasticPosterior(posterior);
	                if (Math.abs(fLogLikelihood - fNonStochasticLogP) > 1e-6) {
	                    reportLogLikelihoods(posterior, "");
	                    System.err.println("At sample " + sampleNr + "\nLikelihood incorrectly calculated: " + fNonStochasticLogP + " != " + fLogLikelihood
	                            + " Operator: " + operator.getClass().getName());
	                }
	                if (sampleNr > NR_OF_DEBUG_SAMPLES * 3) {
	                    // switch off debug mode once a sufficient large sample is checked
	                    debugFlag = false;
	                    if (Math.abs(fLogLikelihood - fNonStochasticLogP) > 1e-6) {
	                        // incorrect calculation outside debug period.
	                        // This happens infrequently enough that it should repair itself after a robust posterior calculation
	                        corrections++;
	                        if (corrections > 100) {
	                            // after 100 repairs, there must be something seriously wrong with the implementation
	                            System.err.println("Too many corrections. There is something seriously wrong that cannot be corrected");
	                            state.storeToFile(sampleNr);
	                            operatorSchedule.storeToFile();
	                            System.exit(0);
	                        }
	                        oldLogLikelihood = state.robustlyCalcPosterior(posterior);;
	                    }
	                } else {
	                    if (Math.abs(fLogLikelihood - fNonStochasticLogP) > 1e-6) {
	                        // halt due to incorrect posterior during intial debug period
	                        state.storeToFile(sampleNr);
	                        operatorSchedule.storeToFile();
	                        System.exit(0);
	                    }
	                }
	            } else {
	                if (sampleNr >= 0) {
	                	operator.optimize(logAlpha);
	                }
	            }
	            callUserFunction(sampleNr);

	            // make sure we always save just before exiting
	            if (storeEvery > 0 && (sampleNr + 1) % storeEvery == 0 || sampleNr == chainLength) {
	                /*final double fLogLikelihood = */
	                state.robustlyCalcNonStochasticPosterior(posterior);
	                state.storeToFile(sampleNr);
	                operatorSchedule.storeToFile();
	            }
	        }
	        if (corrections > 0) {
	            System.err.println("\n\nNB: " + corrections + " posterior calculation corrections were required. This analysis may not be valid!\n\n");
	        }
	        currentSample += resampleEvery;
	        return System.currentTimeMillis();
	}

	public void optimiseRunTime(long startTime, long endTime, long endTimeMainChain) {}
	
}
