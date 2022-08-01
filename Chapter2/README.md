## Who are you? A framework to identify and report genetic sample mix-ups

We surveyed the community and found that almost 80% of responding researchers have encountered sample mix-ups. However, many recent studies in the field of molecular ecology do not appear to systematically report individual assignment checks as part of their publications. Although checks may be done, lack of consistent reporting means that it is difficult to assess whether sample mix-ups have occurred or been detected. Here, we present an easy-to-follow sample verification framework that can utilise existing metadata, including species, population structure, sex and pedigree information (detailed methods in the [publication](https://doi.org/10.1111/1755-0998.13575)).

#### Detection of sample mix-ups
Initial comparison of heterozygosity values of genotyped Z-linked SNP array markers for the Tiritiri Matangi individuals. Later, sample mix-up was confirmed comparing pedigree (microsatellite-based) and genomic relatedness in R.

#### Blank, Control and Duplicate check
Duplicates were identified by the default criteria in the Axiom Analysis Suite software. Sample duplicate individuals can also be identified for example using the --genome function in the PLINK software, in COLONY during analysis set-up or from the sequoia function in the package sequoia in R.

#### Principal Component Analysis 
We used the PLINK function –make-rel to calculate relatedness and the function --pca to get principal components. We used R to plot the first two principal components and the third and fourth principal components.

#### Sex check
Z-linked SNPs were used to confirm the sex of the genotyped individuals by calculating homozygosity by locus for 2,694 sex markers in the 'GENHET' package in R. For comparison, we also made use of the –check-sex function in PLINK, with the default 0.2/0.8 F-statistic threshold.

#### Check for relatedness consistency 
As a next step, we compared the pedigree-based (kinship2) and genomic relatedness (GCTA) of first-degree relationships of the individual hihi in order to identify more sample errors. As PLINK only reports Mendelian errors per family, sequoia only per SNP, and COLONY only indirectly, we instead built a custom-matrix in R. 

#### Parentage assignment 
We used the R package sequoia and the stand-alone software COLONY to reconstruct pedigree relationships for all Tiritiri Matangi individuals. Sequoia makes use of a conservative hill-climbing algorithm to construct high likelihood pedigrees using relatively few markers, however, this heuristic, sequential approach also means that the inferred relationship between two individuals is not necessarily the true relationship. COLONY on the other hand can employ one likelihood or two pairwise likelihood methods for inferring parentage and sibship, and has the advantage of being rather conservative and hence producing trustworthy pedigrees where the assignment confidences are reported separately.

#### The use of whole-genome sequencing data
In order to check that these WGS individuals were correctly allocated on the SNP array, we reduced the WGS files to those sites that are also genotyped with the hihi SNP array. Those files from two different genotyping processes were merged, resulting in a file containing each of the ten individuals twice. PLINK’s --genome function then confirmed that all ten hihi pairs had a pairwise IBD estimate of one, meaning that those ten hihi were assigned the correct ID during the SNP array genotyping.

#### Identifying sources of error 
In order to address whether wet lab errors arose from systematic mis-plating of samples, we mapped individuals across all genotyping plates. For this, we made use of the R package platetools (https://cran.r-project.org/package=platetools) that allowed us to colour code the wells according to the sample classification.
