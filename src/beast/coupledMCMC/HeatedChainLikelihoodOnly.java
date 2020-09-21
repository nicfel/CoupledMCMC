package beast.coupledMCMC;

import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Operator;
import beast.core.util.Evaluator;
import beast.util.Randomizer;

@Description("Like heated chain, but samples log(prior) + beta * log(likelihood) "
		   + "instead of beta * (log(prior) + log(likelihood))")
public class HeatedChainLikelihoodOnly extends HeatedChain {

	
	Distribution likelihood;
	double oldLogP = Double.NaN, newLogP;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		BEASTInterface o = findLikelihood(posteriorInput.get());
		if (o == null) {
			throw new IllegalArgumentException("Expected an element with id='likelihood', but could not find one in the posterior");
		}
		if (!(o instanceof Distribution)) {
			throw new IllegalArgumentException("Expected the element with id='likelihood' to be a Distribution");
		}
		likelihood = (Distribution) o;
	}

	private BEASTInterface findLikelihood(BEASTInterface o) {
		for (BEASTInterface o2 : o.listActiveBEASTObjects()) {
			if (o2.getID() != null && o2.getID().equals("likelihood")) {
				return o2;
			}
			BEASTInterface likelihood = findLikelihood(o2);
			if (likelihood != null) {
				return likelihood;
			}
		}
		return null;
	}
	
	@Override
	protected double getCurrentLogLikelihood() {
		return oldLogLikelihood - oldLogP + this.beta * oldLogP;
	}
	
		
	@Override
	protected double getScaledLogLikelihood(double beta) {
		return oldLogLikelihood - oldLogP + beta * oldLogP;
	}

	
    /**
     * Perform a single MCMC propose+accept/reject step.
     *
     * @param sampleNr the index of the current MCMC step
     * @return the selected {@link beast.core.Operator}
     */
	@Override
    protected Operator propagateState(final long sampleNr) {
		if (Double.isNaN(oldLogP)) {
			oldLogP = likelihood.getCurrentLogP();
		}
		
        state.store(sampleNr);
//            if (m_nStoreEvery > 0 && sample % m_nStoreEvery == 0 && sample > 0) {
//                state.storeToFile(sample);
//            	operatorSchedule.storeToFile();
//            }

        final Operator operator = operatorSchedule.selectOperator();

//        if (printDebugInfo) System.err.print("\n" + sampleNr + " " + operator.getName()+ ":");

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
                    state.store(sampleNr);

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
    		newLogP = likelihood.getCurrentLogP();

            logAlpha = newLogLikelihood - newLogP + newLogP * beta 
            		- (oldLogLikelihood - oldLogP + oldLogP * beta)
            		+ logHastingsRatio; 
            
            if (logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
                // accept
                oldLogLikelihood = newLogLikelihood;
                oldLogP = newLogP;
                state.acceptCalculationNodes();

                if (sampleNr >= 0) {
                    operator.accept();
                }
            } else {
                // reject
                if (sampleNr >= 0) {
                    operator.reject(newLogLikelihood == Double.NEGATIVE_INFINITY ? -1 : 0);
                }
                state.restore();
                state.restoreCalculationNodes();
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
        }
        log(sampleNr);
        return operator;
    }
}
