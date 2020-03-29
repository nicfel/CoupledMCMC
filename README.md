# Coupled MCMC for BEAST2


## Installing the CoupledMCMC package

CoupledMCMC is a [BEAST 2](http://beast2.org) package, which you can install through the 
[package manager](http://www.beast2.org/managing-packages/) that comes with BEAST.
Choose `CoupledMCMC` from the list of packages.

## How to set up your BEAST2 analysis to run with coupled MCMC/parallel tempering 

### By using the conversion BEAUti

If you want to run a Standard analysis (i.e. an analysis using the Standard template for MCMC analyses), you can choose the CoupledMCMC template instead under the `File => Templates => CoupledMCMC` menu. Then you can set up the analysis as you would a Standard analysis, except the MCMC tab is replaced by the CoupledMCMC tab.

If you want to run another analysis than a Standard one, you still can use CoupledMCMC per below.

### By using the conversion app

After you installed the `CoupledMCMC` package (version 0.1.5 or better), the `MCMC2CoupledMCMC` app becomes available in the app launcher.

>
> 1. Create MCMC analysis in BEAUti with any of the available templates, save as `mcmc.xml`
> 
> Now there are 2 ways to proceed:
> 
> 2a. from a terminal, run
>
>  /path/to/beast/bin/applauncher MCMC2CoupledMCMC -xml mcmc.xml -o mc3.xml
>
> This creates a file `mc3.xml` containing a CoupledMCMC analysis with the same model/operators/loggers etc as the `mcmc.xml` analysis.
>
> 2b. from BEAUti, use menu `File > Launch apps`, select `MCMC to Coupled MCMC converter` from the available apps, fill in form and click OK
>


### By editing an XML file in a text editor

In order to set up a pre-prepared xml to run with coupled MCMC, open the  `*.xml` and change the MCMC line in the xml.

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

Müller, Nicola Felix, and Remco Bouckaert. "Coupled MCMC in BEAST 2." bioRxiv (2019): 603514. ([abstract](https://www.biorxiv.org/content/10.1101/603514v1.abstract), [pdf](https://www.biorxiv.org/content/biorxiv/early/2019/04/09/603514.full.pdf));


Parallel Metropolis coupled Markov chain Monte Carlo for Bayesian phylogenetic inference 
[https://academic.oup.com/bioinformatics/article/20/3/407/186341](https://academic.oup.com/bioinformatics/article/20/3/407/186341)
