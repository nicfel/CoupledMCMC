# Source code for coupled MCMC implementation

Implements coupled MCMC for Beast2

## Citation

Parallel Metropolis coupled Markov chain Monte Carlo for Bayesian phylogenetic inference 
[https://academic.oup.com/bioinformatics/article/20/3/407/186341](https://academic.oup.com/bioinformatics/article/20/3/407/186341)

### Set up the xml to run two chains

In order to setup the analysis to run with coupled MCMC, we have to open the  `*.xml` and change one line in the xml.
To do so, go to the line with:
```
<run id="mcmc" spec="MCMC" chainLength="....." numInitializationAttempts="....">
```
To have a run with coupled MCMC, we have to replace that one line with:
```
<run id="mcmc" spec="beast.coupledMCMC.CoupledMCMC" logHeatedChains="true" chainLength="100000000" storeEvery="1000000" deltaTemperature="0.025" chains="2" resampleEvery="10000">
```
* `logHeatedChains="true"` also logs the log files of the heated chains if true.
* `chainLength="100000000"` defines for how many iterations the chains is run
* `deltaTemperature="0.025"` defines the temperature difference between the chain *n* and chain *n-1*. This value should be changed such that the acceptance probability of a swap is between 0.25 and 0.6
* `chains="2"` defines the number of parallel chains that are run. The first chain is the one that explores the posterior just like a normal MCMC chain. All other chains are what's called *heated*. This means that MCMC moves of those chains have a higher probability of being accepted. While these heated chains don't explore the posterior properly, they can be used to propose new states to the one cold chain.   

