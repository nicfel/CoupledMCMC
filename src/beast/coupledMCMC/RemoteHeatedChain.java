package beast.coupledMCMC;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.net.Socket;
import java.net.UnknownHostException;
import java.util.Collections;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Logger;
import beast.core.Operator;
import beast.core.OperatorSchedule;
import beast.core.StateNodeInitialiser;
import beast.core.util.Evaluator;
import beast.core.util.Log;
import beast.util.Randomizer;

@Description("Local representation of heated chain running on a remote machine")
public class RemoteHeatedChain {
	// initialize socket and input output streams 
	private Socket socket		 = null; 
	private DataInputStream input = null; 
	private DataInputStream sin = null; 
	private DataOutputStream out	 = null; 

	
	private double beta;
	private int chainNr = 0;
	private double oldLogLikelihood = 0;

    public RemoteHeatedChain(String xml, String host, int port){
    	// loggersInput.setRule(Input.Validate.FORBIDDEN);
    	

    		// establish a connection 
    		try
    		{ 
    			socket = new Socket(host, port); 
    			System.out.println("Connected"); 

    			// takes input from terminal 
    			input = new DataInputStream(System.in); 

    			// sends output to the socket 
    			out = new DataOutputStream(socket.getOutputStream()); 
    			sin = new DataInputStream(socket.getInputStream());
    			out.writeUTF(xml);
    			String response = sin.readUTF();
    			System.err.println(response);
    		} 
    		catch(UnknownHostException u) 
    		{ 
    			System.out.println(u); 
    		} 
    		catch(IOException i) 
    		{ 
    			System.out.println(i); 
    		} 

    	
    }
    


	protected double getCurrentLogLikelihood() {
		return oldLogLikelihood * beta;
	}
	
	protected double getUnscaledCurrentLogLikelihood() {
		try {
			out.writeUTF("getu");
			String response = sin.readUTF();
			System.err.println(response);
			return Double.parseDouble(response);
		} catch (IOException e) {
			e.printStackTrace();
		} 
		return Double.NaN;
	}

	
	protected void setResampleEvery(int resampleEvery) {
		try {
			out.writeUTF("setr," + resampleEvery);
			String response = sin.readUTF();
			System.err.println(response);
		} catch (IOException e) {
			e.printStackTrace();
		} 
//		this.resampleEvery = resampleEvery;
	}
		
//	// RRB: what is the purpose of the first argument?
	protected void setTemperature(int i, double temperature) {
		try {
			out.writeUTF("setb," + 1/(1 + temperature));
			String response = sin.readUTF();
			System.err.println(response);
		} catch (IOException e) {
			e.printStackTrace();
		} 
		this.beta = 1/(1 + temperature);
	}

	protected double getBeta(){
		return beta;
	}	
	
	protected void setBeta(double beta){
		try {
			out.writeUTF("setb," + beta);
			String response = sin.readUTF();
			System.err.println(response);
		} catch (IOException e) {
			e.printStackTrace();
		} 
		this.beta = beta;		
	}	
	
	protected int getChainNr(){
		return chainNr;
	}	
	
	protected void setChainNr(int chainNr){
		this.chainNr = chainNr;
		try {
			out.writeUTF("seti," + chainNr);
			String response = sin.readUTF();
			System.err.println(response);
		} catch (IOException e) {
			e.printStackTrace();
		} 
	}	

//	protected String getStateFileName(){
//		return stateFileName;
//	}
//	
//	protected void setStateFileName(String stateFileName){
//		this.stateFileName = stateFileName;
//		state.setStateFileName(stateFileName);
//	}
//
//	
//	protected double calcCurrentLogLikelihoodRobustly() {
//		oldLogLikelihood = robustlyCalcPosterior(posterior);
//		return getCurrentLogLikelihood();
//	};
//
//	
//    @Override
//    public void run() throws IOException, SAXException, ParserConfigurationException {
//        // set up state (again). Other plugins may have manipulated the
//        // StateNodes, e.g. set up bounds or dimensions
//        state.initAndValidate();
//        // also, initialise state with the file name to store and set-up whether to resume from file
//        state.setStateFileName(stateFileName);
//        operatorSchedule.setStateFileName(stateFileName);
//
//        burnIn = burnInInput.get();
//        chainLength = chainLengthInput.get();
//        int nInitialisationAttempts = 0;
//        state.setEverythingDirty(true);
//        posterior = posteriorInput.get();
//
//        if (restoreFromFile) {
//            state.restoreFromFile();
//            operatorSchedule.restoreFromFile();
//            burnIn = 0;
//            oldLogLikelihood = state.robustlyCalcPosterior(posterior);
//        } else {
//            do {
//                for (final StateNodeInitialiser initialiser : initialisersInput.get()) {
//                    initialiser.initStateNodes();
//                }
//                oldLogLikelihood = state.robustlyCalcPosterior(posterior);
//                nInitialisationAttempts += 1;
//            } while (Double.isInfinite(oldLogLikelihood) && nInitialisationAttempts < numInitializationAttempts.get());
//        }
//        final long startTime = System.currentTimeMillis();
//        
//
//        // do the sampling
//        logAlpha = 0;
//        debugFlag = Boolean.valueOf(System.getProperty("beast.debug"));
//
//        System.err.println("Start likelihood: " + oldLogLikelihood + " " + (nInitialisationAttempts > 1 ? "after " + nInitialisationAttempts + " initialisation attempts" : ""));
//        if (Double.isInfinite(oldLogLikelihood) || Double.isNaN(oldLogLikelihood)) {
//            reportLogLikelihoods(posterior, "");
//            throw new RuntimeException("Could not find a proper state to initialise. Perhaps try another seed.");
//        }
//
//        loggers = coupledLoggersInput.get();
//
//        // put the loggers logging to stdout at the bottom of the logger list so that screen output is tidier.
//        Collections.sort(loggers, (o1, o2) -> {
//            if (o1.isLoggingToStdout()) {
//                return o2.isLoggingToStdout() ? 0 : 1;
//            } else {
//                return o2.isLoggingToStdout() ? -1 : 0;
//            }
//        });
//        // warn if none of the loggers is to stdout, so no feedback is given on screen
//        boolean hasStdOutLogger = false;
//        boolean hasScreenLog = false;
//        for (CoupledLogger l : loggers) {
//        	if (l.isLoggingToStdout()) {
//        		hasStdOutLogger = true;
//        	}
//        	if (l.getID() != null && l.getID().equals("screenlog")) {
//        		hasScreenLog = true;
//        	}
//        }
//        if (!hasStdOutLogger) {
//        	Log.warning.println("WARNING: If nothing seems to be happening on screen this is because none of the loggers give feedback to screen.");
//        	if (hasScreenLog) {
//        		Log.warning.println("WARNING: This happens when a filename  is specified for the 'screenlog' logger.");
//        		Log.warning.println("WARNING: To get feedback to screen, leave the filename for screenlog blank.");
//        		Log.warning.println("WARNING: Otherwise, the screenlog is saved into the specified file.");
//        	}
//        }
//
//        // initialises log such that log file headers are written, etc.
//        for (final Logger log : loggers) {
//            log.init();
//        }
//        
////        Randomizer = new Randomizer();
//                
//    } // run;
//    
//    public OperatorSchedule getOperatorSchedule(){
//    	return operatorSchedule;
//    }    
//	
//	// run MCMC inner loop for resampleEvery nr of samples
	protected long runTillResample(long runUntil) throws Exception {
		try {
			out.writeUTF("runtill," + runUntil);
			String response = sin.readUTF();
			System.err.println(response);
			oldLogLikelihood = Double.parseDouble(response);
		} catch (IOException e) {
			e.printStackTrace();
		} 
		return System.currentTimeMillis();
	}
	
	public String getOperatorTuningInfo() {
		try {
			out.writeUTF("geto");
			String response = sin.readUTF();
			System.err.println(response);
			return response;
		} catch (IOException e) {
			e.printStackTrace();
		} 
		return null;
	}
	
	public void setOperatorTuningInfo(String tuning) {
		try {
			out.writeUTF("seto" + tuning);
			String response = sin.readUTF();
			System.err.println(response);
		} catch (IOException e) {
			e.printStackTrace();
		} 
	}

//		int corrections = 0;
//		final boolean isStochastic = posterior.isStochastic();
//		
//
//		for (long sampleNr = currentSample; sampleNr <= runUntil; sampleNr++) {
//            final long currentState = sampleNr;
//            final Operator operator = propagateState(sampleNr);
////            System.out.println(chainNr + " " + sampleNr);
//            
//            if (debugFlag && sampleNr % 3 == 0 || sampleNr % 10000 == 0) {
//                // check that the posterior is correctly calculated at every third
//                // sample, as long as we are in debug mode
//            	final double originalLogP = isStochastic ? posterior.getNonStochasticLogP() : oldLogLikelihood;
//                final double logLikelihood = isStochastic ? state.robustlyCalcNonStochasticPosterior(posterior) : state.robustlyCalcPosterior(posterior);
//                if (isTooDifferent(logLikelihood, originalLogP)) {
//                    reportLogLikelihoods(posterior, "");
//                    Log.err.println("At sample " + sampleNr + "\nPosterior incorrectly calculated: " + originalLogP + " != " + logLikelihood
//                    		+ "(" + (originalLogP - logLikelihood) + ")"
//                            + " Operator: " + operator.getClass().getName());
////                    System.exit(0);
//                }
//                if (sampleNr > NR_OF_DEBUG_SAMPLES * 3) {
//                    // switch off debug mode once a sufficient large sample is checked
//                    debugFlag = false;
//                    if (isTooDifferent(logLikelihood, originalLogP)) {
//                        // incorrect calculation outside debug period.
//                        // This happens infrequently enough that it should repair itself after a robust posterior calculation
//                        corrections++;
//                        if (corrections > 100) {
//                            // after 100 repairs, there must be something seriously wrong with the implementation
//                        	Log.err.println("Too many corrections. There is something seriously wrong that cannot be corrected");
//                            state.storeToFile(sampleNr);
//                            operatorSchedule.storeToFile();
//                            System.exit(1);
//                        }
//                        oldLogLikelihood = state.robustlyCalcPosterior(posterior);
//                    }
//                } else {
//                    if (isTooDifferent(logLikelihood, originalLogP)) {
//                        // halt due to incorrect posterior during intial debug period
//                        state.storeToFile(sampleNr);
//                        operatorSchedule.storeToFile();
//                        System.exit(1);
//                    }
//                }
//            } else {
//                if (sampleNr >= 0) {
//                	operator.optimize(logAlpha);
//                }
//            }
//            callUserFunction(sampleNr);
//            
//           
////            // make sure we always save just before exiting
//            if (storeEvery > 0 && (sampleNr+1) % storeEvery == 0 || sampleNr == chainLength) {
//                /*final double logLikelihood = */
//                state.robustlyCalcNonStochasticPosterior(posterior);
//                state.storeToFile(sampleNr);
//                operatorSchedule.storeToFile();
//            }
//            
//
//            
//            if (posterior.getCurrentLogP() == Double.POSITIVE_INFINITY) {
//            	throw new RuntimeException("Encountered a positive infinite posterior. This is a sign there may be numeric instability in the model.");
//            }
//        }
//        if (corrections > 0) {
//        	Log.err.println("\n\nNB: " + corrections + " posterior calculation corrections were required. This analysis may not be valid!\n\n");
//        }
//
//        // the plus 1 is required such that the currentSample is not recomputed
//        currentSample = runUntil+1;
//        
//        return System.currentTimeMillis();
//	}
//	
//    /**
//     * Perform a single MCMC propose+accept/reject step.
//     *
//     * @param sampleNr the index of the current MCMC step
//     * @return the selected {@link beast.core.Operator}
//     */
//	@Override
//    protected Operator propagateState(final long sampleNr) {
//        state.store(sampleNr);
////            if (m_nStoreEvery > 0 && sample % m_nStoreEvery == 0 && sample > 0) {
////                state.storeToFile(sample);
////            	operatorSchedule.storeToFile();
////            }
//
//        final Operator operator = operatorSchedule.selectOperator();
//
////        if (printDebugInfo) System.err.print("\n" + sampleNr + " " + operator.getName()+ ":");
//
//        final Distribution evaluatorDistribution = operator.getEvaluatorDistribution();
//        Evaluator evaluator = null;
//
//        if (evaluatorDistribution != null) {
//            evaluator = new Evaluator() {
//                @Override
//                public double evaluate() {
//                    double logP = 0.0;
//
//                    state.storeCalculationNodes();
//                    state.checkCalculationNodesDirtiness();
//
//                    try {
//                        logP = evaluatorDistribution.calculateLogP();
//                    } catch (Exception e) {
//                        e.printStackTrace();
//                        System.exit(1);
//                    }
//
//                    state.restore();
//                    state.store(sampleNr);
//
//                    return logP;
//                }
//            };
//        }
//        final double logHastingsRatio = operator.proposal(evaluator);
//
//        if (logHastingsRatio != Double.NEGATIVE_INFINITY) {
//
//            if (operator.requiresStateInitialisation()) {
//                state.storeCalculationNodes();
//                state.checkCalculationNodesDirtiness();
//            }
//
//            newLogLikelihood = posterior.calculateLogP();
//
//            logAlpha = newLogLikelihood*beta - oldLogLikelihood*beta + logHastingsRatio; 
//            
//            if (logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
//                // accept
//                oldLogLikelihood = newLogLikelihood;
//                state.acceptCalculationNodes();
//
//                if (sampleNr >= 0) {
//                    operator.accept();
//                }
//            } else {
//                // reject
//                if (sampleNr >= 0) {
//                    operator.reject(newLogLikelihood == Double.NEGATIVE_INFINITY ? -1 : 0);
//                }
//                state.restore();
//                state.restoreCalculationNodes();
//            }
//            state.setEverythingDirty(false);
//        } else {
//            // operation failed
//            if (sampleNr >= 0) {
//                operator.reject(-2);
//            }
//            state.restore();
//            if (!operator.requiresStateInitialisation()) {
//                state.setEverythingDirty(false);
//                state.restoreCalculationNodes();
//            }
//        }
//        log(sampleNr);
//        return operator;
//    }
//
//    private boolean isTooDifferent(double logLikelihood, double originalLogP) {
//    	//return Math.abs((logLikelihood - originalLogP)/originalLogP) > 1e-6;
//    	return Math.abs(logLikelihood - originalLogP) > 1e-6;
//	}
//
//	public void optimiseRunTime(long startTime, long endTime, long endTimeMainChain) {}
//	
//	protected long getCurrentSample(){
//		return currentSample;
//	}
//	
//	@Override
//    public void log(final long sampleNr) {
//        for (final CoupledLogger log : loggers) {
//            log.log(sampleNr);
//        }
//    } // log
//
//
//    @Override
//    /**
//     * Set up information related to the file for (re)storing the State.
//     * The Runnable implementation is responsible for making its
//     * State synchronising with the file *
//     * @param fileName
//     * @param isRestoreFromFile
//     */
//    public void setStateFile(final String fileName, final boolean isRestoreFromFile) {
//
//        stateFileName = fileName;
//        restoreFromFile = isRestoreFromFile;
//        
//        state.setStateFileName(stateFileName);
//    }
//	
//
}
