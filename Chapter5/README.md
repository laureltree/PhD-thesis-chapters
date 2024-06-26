
## Genomic signatures of inbreeding depression for a threatened Aotearoa New Zealand passerine

Here, we investigate inbreeding and inbreeding depression for hihi (*Notiomystis cincta*). Firstly, using a custom 45k SNP array, we explore genomic inbreeding patterns by inferring homozygous segments across the genome. Secondly, we measure genomic inbreeding effects on fitness, measured as lifetime reproductive success, and its three components, annual reproductive success, juvenile survival and adult annual survival, in 363 hihi. This manuscript has been submitted for publication and is currently being peer-reviewed.

#### ROH-based inbreeding in hihi
Inbreeding was measured from the SNP data using the hidden Markov model-based approach in the R package RZooRoH  that identifies homozygous-by-decent segments and allows for the estimation of a global inbreeding coefficient. For each individual, HBD probabilities were summed over the first ten HBD classes to give individual inbreeding coefficients for all birds. We further divided the whole-genome inbreeding level into very recent, middle and ancient inbreeding. We note that a probability-based RZooRoH approach will, on average, yield higher inbreeding values than binary estimates that are offered e.g. by PLINK, but our previous work indicates that these values are highly correlated for hihi.

#### ROH density across the population
We measured the average ROH density across the genome for all 363 hihi individuals by extracting all HBD segments per chromosome in RZooRoH and estimating mean HBD probabilities of all markers in non-overlapping 500kb windows using the R package windowscanr, following R code provided by [Stoffel et al.](https://github.com/mastoffel/sheep_ID). 

#### MCMCglmm modelling
The effects of whole-genome genomic inbreeding, sex and lifespan on lifetime reproductive success (LRS), the effects of inbreeding, age and sex on annual reproductive success (ARS) and the effects of genomic inbreeding and sex on juvenile and adult annual survival (JUS, ADS) were tested using MCMCglmm, an R package that fits generalised linear mixed models using Markov chain Monte Carlo techniques. We also fitted additional models with the separated inbreeding values, as very recent and middle inbreeding seemed to contribute most to inbreeding levels in highly inbred birds. Convergence was checked graphically and with the Heidelberger and Welch convergence test using the coda R package.

#### ROH genome wide association scan (GWAS)
We also used the list of identified HBD-segments larger than 300kb from the RZooRoH analysis to test for association between an allele of a SNP being in a ROH and hihi fitness. All R scripts and models regarding this mapping of inbreeding depression are modified from [Stoffel et al.](https://github.com/mastoffel/sheep_ID) unless otherwise indicated, and can be found in the data availability section of their publication together with a detailed description of their methods.
