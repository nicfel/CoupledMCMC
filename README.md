# Coupled MCMC for BEAST2
## How to set-up your BEAST2 analysis to run with coupled MCMC/parallel tempering  

In order to setupa pre-prepared xml to run with coupled MCMC, open the  `*.xml` and change the MCMC line in the xml.

To do so, go to the line with:

```
<run id="mcmc" spec="MCMC" chainLength="....." numInitializationAttempts="....">
```

To have a run with coupled MCMC, we have to replace that one line with:

```
<run id="mcmc" spec="beast.coupledMCMC.CoupledMCMC" chainLength="100000000" storeEvery="1000000" deltaTemperature="0.025" chains="2" resampleEvery="10000">
```
* `chainLength="100000000"` defines for how many iterations the chains is run
* `deltaTemperature="0.025"` defines the temperature difference between the chain *n* and chain *n-1*. This value should be changed such that the acceptance probability of a swap is between 0.25 and 0.6
* `chains="2"` defines the number of parallel chains that are run. The first chain is the one that explores the posterior just like a normal MCMC chain. All other chains are what's called *heated*. This means that MCMC moves of those chains have a higher probability of being accepted. While these heated chains don't explore the posterior properly, they can be used to propose new states to the one cold chain.   


## Citation

Parallel Metropolis coupled Markov chain Monte Carlo for Bayesian phylogenetic inference 
[https://academic.oup.com/bioinformatics/article/20/3/407/186341](https://academic.oup.com/bioinformatics/article/20/3/407/186341)
