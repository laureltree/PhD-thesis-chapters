
### UoA ###

## RepeatAbel ####

## L Duntsch



############## TARSUS ################


### install repeatabel
install.packages("hglm")
# install.packages("RepeatABEL_1.1.tar.gz", repos = NULL)

## 
# The package performs GWAS for data where there are repeated observations on individuals
# that may be related. Various random effects can be fitted. Random polygenic
# effects are fitted by default, but also other random effects can be fitted including spatial effects. 
# A model without the fixed SNP effect is fitted using the
# hglm package for estimating variance components (see Ronnegard, Shen & Alam (2010) hglm: 
#                                 A package for fitting hierarchical generalized linear models. The R Journal 2:20-28). 
# These variance component estimates are subsequently used
# in the GWAS where each marker is fitted one at a time. Thus, this is similar
# to using the polygenic_hglm function and subsequently the mmscore function
# in GenABEL. The package consists of three main functions rGLS, preFitModel
# and simulate_PhenData


library(RepeatABEL)
library(hglm)
library(data.table)

# Our proposed method is a two‐stage method where the distribution of the residuals and random effects 
# are estimated in a preliminary model that does not include the SNP effects.
# In the second stage, a model including the SNP effect is fitted.

#  Although the focus of the package is on models for related individuals having repeated observations, 
# it can be used for linear models having arbitrary covariance structures.


#################################################################################


########## using it on 523 individuals

## 18 individuals removed because of missing data


#################################################################################

## done once:
## Phenotype data; more complete set

Pheno_fixed_random <- fread("data_subset_523.csv", header = T)
Pheno_fixed_random[1:10,]

## retain only the ones we need
Pheno_fixed_random1 <- Pheno_fixed_random[,c("FldgTarsus", "FldgMass" , "FldgHead.bill",  "ClutchSize", "id", 
                                             "Dam.ID", "Social.Sire.ID", "MonthYear", "Year")]
Pheno_fixed_random1[1:10,]
dim(Pheno_fixed_random1)

## get the sex as well from the tfam again
tfam <- fread("hihi_confirmed_31012020_noZ.tfam", header = F)
tfam1 <- tfam[,c(2,5)] 
colnames(tfam1)[1] <- "id"
colnames(tfam1)[2] <- "sex"
tfam1

## merge into full table
Pheno_fixed_random2 <- merge(tfam1, Pheno_fixed_random1, by = "id", all.y = TRUE)
dim(Pheno_fixed_random2)

## fix sex: change female = 0 and male = 1
Pheno_fixed_random2$sex[Pheno_fixed_random2$sex == 2] <- 0
Pheno_fixed_random2[1:10,]
dim(Pheno_fixed_random2)

write.table(Pheno_fixed_random2, file = "Pheno_fixed_random_523.phen", col.names = TRUE, row.names = FALSE, quote = FALSE)

#################################################################################


## and we need a new .out file with the reduced amount of individuals in PLINK
# make genotype file based on transposed plink format .tped and .tmap
convert.snp.tped('hihi_curated_523_noZ.tped', 
                 'hihi_curated_523_noZ.tfam', 'Hihi_6Feb2020.out', strand = "u", bcast = 10000)

# MAKE GWAA FILE (GENABEL DATA FILE)
#combine the genabel genotype file with the phenotypic data using the "load.gwaa.data" function - you need to specify : the name of the phenotype file and the name of the genotype file
phen_gwaa<-load.gwaa.data(phenofile = 'Pheno_fixed_random_523.phen', 
                          genofile = 'Hihi_6Feb2020.out', force = TRUE, makemap = FALSE, sort = TRUE)
phen_gwaa[1:10,1:10]
# The data frame PD includes the trait value y, covariates (e.g. sex) 
# and the ID of the individuals. 


# vec <- c("AX-172486735", "AX-171754927", "AX-171714808", "AX-171722226", "AX-171764222", "AX-171712697",
#                                        "AX-171757070", "AX-171730907", "AX-171755458", "AX-171736922")

# vec <- which(Top_10_SNPs_Tarsus$rn == "AX-172486735" | Top_10_SNPs_Tarsus$rn == "AX-171754927" |
#                Top_10_SNPs_Tarsus$rn == "AX-171714808" | Top_10_SNPs_Tarsus$rn == "AX-171722226" |
#                Top_10_SNPs_Tarsus$rn == "AX-171764222" | Top_10_SNPs_Tarsus$rn ==  "AX-171712697"|
#                Top_10_SNPs_Tarsus$rn == "AX-171757070" | Top_10_SNPs_Tarsus$rn == "AX-171730907" | 
#                Top_10_SNPs_Tarsus$rn == "AX-171755458"| Top_10_SNPs_Tarsus$rn == "AX-171736922")


#################################################################################


# SOME USEFUL FUNCTIONS
nids(phen_gwaa)  #523 - Number of individuals in the data set
snps<-nsnps(phen_gwaa) #42433 - Number of SNPs in the data set
descriptives.trait(phen_gwaa) # Summary of the phenotypic data
s <- summary(gtdata(phen_gwaa)) # Summary of the genotypic data (chromosomeID, call rates, HWE)
s[1:5, ]

PD <- fread("Pheno_fixed_random_523.phen", header = T)
PD[1:10,]
dim(PD)
# The data frame Phen.Data includes the trait value y, two covariates (age and sex) 
# and the ID of the individuals. 


###################################################
######### start analysis ########
###################################################


####### Using the rGLS function ###########

# estimate genomic heritability of the phenotypic trait using the function ‘rGLS’
# We wish to include age and sex as fixed effects so the function input is:
?rGLS
GWAS1 <- rGLS(FldgTarsus ~ sex + ClutchSize, genabel.data = phen_gwaa,
              phenotype.data = PD)
# [1] "GRM ready"
# [1] "Variance component estimation ready"
# [1] "Rotation matrix ready"
# [1] "Rotate LMM started"
# [1] "Rotate LMM ready"

summary(GWAS1)

## Extracting the genotypic variance, a 95% CI and the heritability

# From the GWAS1 object, the estimated variance components for the prefitted
# model, not including SNP effects, can be extracted as follows.

est.hglm <- GWAS1@call$hglm
cat("Genotypic and permanent env. variance components:","\n",
    est.hglm$varRanef,", resp.","\n",
    "The residual variance is", est.hglm$varFix,".","\n")
# Genotypic and permanent env. variance components:
#   0.2099345 0.1804772 , resp.
# The residual variance is  0.552378 .

##  computations for the genotypic variance
logVCE <- est.hglm$SummVC2[[1]]
cat("Estimate and SE for the genotypic variance on a natural log scale:",
    "\n", logVCE[1],"and",logVCE[2],", resp.", "\n \n",
    "Confidence interval: [", exp( logVCE[1] - 1.96*logVCE[2] ) , "," ,
    exp( logVCE[1] + 1.96*logVCE[2] ) , "]" , "\n")

## heritability:
cat("Heritability:",
    est.hglm$varRanef[1]/(est.hglm$varFix + sum(est.hglm$varRanef)),
    "\n")
## 0.223975



##### OR:

##### Using the preFitModel function ##################

## Fitting the same model as above
# The same results can be computed using the preFitModel as follows
fixed=FldgTarsus ~ sex + ClutchSize
?preFitModel
Mod1 <- preFitModel(fixed, random=~1|id,
                    genabel.data = phen_gwaa,
                    phenotype.data = PD,
                    corStruc=list( id=list("GRM","Ind") ))
#[1] "GRM ready"

GWAS1b <- rGLS(fixed, genabel.data = phen_gwaa,
               phenotype.data = PD, V = Mod1$V)
lambda(GWAS1b) #1.049009
# In the example we wish to fit polygenic effects and permanent environmental effects. The former
# requires a correlatioon structure given by the GRM whereas the permanent environmental effects are iid



## A model having several different random effects!!

# A model including polygenic effects,
# permanent environmental effects, and e.g. Year effect as random
Mod2 <- preFitModel(FldgTarsus ~ sex + ClutchSize, random=~1|id + 1|Year + 1|Dam.ID + 1|Social.Sire.ID + 1|MonthYear,
                    genabel.data = phen_gwaa, phenotype.data = PD,
                    corStruc=list( id=list("GRM") , Year=list("Ind") , 
                                   Dam.ID=list("Ind") , Social.Sire.ID=list("Ind") , MonthYear=list("Ind")) )
GWAS2 <- rGLS(FldgTarsus ~ sex + ClutchSize, genabel.data = phen_gwaa,
              phenotype.data = PD, V = Mod2$V)   #If V is not specified then a model including random polygenic effects 
# and permanent environmental effects is fitted (using the hglm package) to compute V

# GWAS1@call$hglm = Mod2$fitted.hglm
est.hglm <- Mod2$fitted.hglm
cat("Genotypic and permanent env. variance components:","\n",
    est.hglm$varRanef,", resp.","\n",
    "The residual variance is", est.hglm$varFix,".","\n")
##  computations for the genotypic variance
logVCE <- est.hglm$SummVC2[[1]]
cat("Estimate and SE for the genotypic variance on a natural log scale:",
    "\n", logVCE[1],"and",logVCE[2],", resp.", "\n \n",
    "Confidence interval: [", exp( logVCE[1] - 1.96*logVCE[2] ) , "," ,
    exp( logVCE[1] + 1.96*logVCE[2] ) , "]" , "\n")

## heritability:
cat("Heritability:",
    est.hglm$varRanef[1]/(est.hglm$varFix + sum(est.hglm$varRanef)),
    "\n")
## 0.08177321


lambda(GWAS2) #1.040771
estlambda(GWAS2[, "P1df"], plot=TRUE)


# plot(GWAS1, main="GWAS1 results")
# plot(GWAS2, main="GWAS2 results")

## done


#############################################################################################


### Body mass as fledgling

#Mod3 +
GWAS3 <- rGLS(FldgMass ~ sex + ClutchSize, genabel.data = phen_gwaa,
              phenotype.data = PD)

est.hglm <- GWAS3@call$hglm
## heritability:
cat("Heritability:",
    est.hglm$varRanef[1]/(est.hglm$varFix + sum(est.hglm$varRanef)),
    "\n")
## 0.1991406

Mod4 <- preFitModel(FldgMass ~ sex + ClutchSize, random=~1|id + 1|Year + 1|Dam.ID + 1|Social.Sire.ID + 1|MonthYear,
                    genabel.data = phen_gwaa, phenotype.data = PD,
                    corStruc=list( id=list("GRM") , Year=list("Ind") , 
                                   Dam.ID=list("Ind") , Social.Sire.ID=list("Ind") , MonthYear=list("Ind")) )
GWAS4 <- rGLS(FldgMass ~ sex + ClutchSize, genabel.data = phen_gwaa,
              phenotype.data = PD, V = Mod4$V)

# GWAS1@call$hglm = Mod2$fitted.hglm
est.hglm <- Mod4$fitted.hglm
cat("Genotypic and permanent env. variance components:","\n",
    est.hglm$varRanef,", resp.","\n",
    "The residual variance is", est.hglm$varFix,".","\n")
##  computations for the genotypic variance
logVCE <- est.hglm$SummVC2[[1]]
cat("Estimate and SE for the genotypic variance on a natural log scale:",
    "\n", logVCE[1],"and",logVCE[2],", resp.", "\n \n",
    "Confidence interval: [", exp( logVCE[1] - 1.96*logVCE[2] ) , "," ,
    exp( logVCE[1] + 1.96*logVCE[2] ) , "]" , "\n")

## heritability:
cat("Heritability:",
    est.hglm$varRanef[1]/(est.hglm$varFix + sum(est.hglm$varRanef)),
    "\n")
## 0.06540873


lambda(GWAS4) #0.9925516
estlambda(GWAS4[, "P1df"], plot=TRUE)

## plot it! and compare

# plot(GWAS3, main="GWAS3 results")
# plot(GWAS4, main="GWAS4 results")


#############################################################################################


### head.bill length as fledgling

#Mod5 +
GWAS5 <- rGLS(FldgHead.bill ~ sex + ClutchSize, genabel.data = phen_gwaa,
              phenotype.data = PD)

est.hglm <- GWAS5@call$hglm
## heritability:
cat("Heritability:",
    est.hglm$varRanef[1]/(est.hglm$varFix + sum(est.hglm$varRanef)),
    "\n")
## 0.1692283

Mod6 <- preFitModel(FldgHead.bill ~ sex + ClutchSize, random=~1|id + 1|Year + 1|Dam.ID + 1|Social.Sire.ID + 1|MonthYear,
                    genabel.data = phen_gwaa, phenotype.data = PD,
                    corStruc=list( id=list("GRM") , Year=list("Ind") , 
                                   Dam.ID=list("Ind") , Social.Sire.ID=list("Ind") , MonthYear=list("Ind")) )
GWAS6 <- rGLS(FldgHead.bill ~ sex + ClutchSize, genabel.data = phen_gwaa,
              phenotype.data = PD, V = Mod6$V)

# GWAS1@call$hglm = Mod2$fitted.hglm
est.hglm <- Mod6$fitted.hglm
cat("Genotypic and permanent env. variance components:","\n",
    est.hglm$varRanef,", resp.","\n",
    "The residual variance is", est.hglm$varFix,".","\n")
##  computations for the genotypic variance
logVCE <- est.hglm$SummVC2[[1]]
cat("Estimate and SE for the genotypic variance on a natural log scale:",
    "\n", logVCE[1],"and",logVCE[2],", resp.", "\n \n",
    "Confidence interval: [", exp( logVCE[1] - 1.96*logVCE[2] ) , "," ,
    exp( logVCE[1] + 1.96*logVCE[2] ) , "]" , "\n")

## heritability:
cat("Heritability:",
    est.hglm$varRanef[1]/(est.hglm$varFix + sum(est.hglm$varRanef)),
    "\n")
## 0.07069575

lambda(GWAS6) #0.9738626                    #values up to 1.10 are generally considered acceptable for GWAS
estlambda(GWAS6[, "P1df"], plot=TRUE)

## plot it! and compare

# plot(GWAS5, main="GWAS5 results")
# plot(GWAS6, main="GWAS6 results")





#### more suitable plotting

plot(GWAS2, ylim = c(1,6), main="GWAS in RepeatABEL: tarsus length", pch=20, col=c("#af8dc3","#7fbf7b"))
plot(GWAS4, ylim = c(1,6), main="GWAS in RepeatABEL: body mass", pch=20, col=c("#af8dc3","#7fbf7b"))
plot(GWAS6, ylim = c(1,6), main="GWAS in RepeatABEL: head-bill length", pch=20, col=c("#af8dc3","#7fbf7b"))

#####

## make pdf
pdf("GWAS_in_RepeatABEL_tarsus_length.pdf", useDingbats = FALSE, width = 12, height = 7)
{
  plot1<-plot(GWAS2, ylim = c(1,10), 
              main = "GWAS in RepeatABEL: tarsus length",
              pch=20,
              col=c("#af8dc3","#7fbf7b"))
  
  abline(h= c(-log10(0.05), -log10(0.05/snps)), lty= 2, col='black')
}
dev.off()

###  plot all in one go ##########

# Reported P-values are based on Wald tests and are corrected for
# population stratification (structure and relatedness) and the repeated sampling
# of the same individuals when using the repeated measures GWAS

plot1<-plot(GWAS2, ylim = c(1,8), 
            main = "GWAS in RepeatABEL: tarsus length",
            pch=20,
            col=c("#af8dc3","#7fbf7b"))
plot1
abline(h= c(-log10(0.05), -log10(0.05/snps)), lty= 2, col='black')  #Bonferroni correction 1.178328e-06

# The summary for the SNPs, which show most significant association can be produced with
results.GWAS_tl<-descriptives.scan(GWAS2, top=10)
results.GWAS_tl


plot2<-plot(GWAS4, ylim = c(1,8), 
            main = "GWAS in RepeatABEL: body mass",
            pch=20,
            col=c("#af8dc3","#7fbf7b"))
plot2
abline(h= c(-log10(0.05), -log10(0.05/snps)), lty= 2, col='black')  #Bonferroni correction 1.178328e-06
# This is a conservative P-value as it assumes that all markers are independent.

# The summary for the SNPs, which show most significant association can be produced with
results.GWAS_bm<-descriptives.scan(GWAS4, top=10)
results.GWAS_bm

plot3<-plot(GWAS6, ylim = c(1,8), 
            main = "GWAS in RepeatABEL: head-bill length",
            pch=20,
            col=c("#af8dc3","#7fbf7b"))
plot3
abline(h= c(-log10(0.05), -log10(0.05/snps)), lty= 2, col='black')  #Bonferroni correction 1.178328e-06

# The summary for the SNPs, which show most significant association can be produced with
results.GWAS_hb<-descriptives.scan(GWAS6, top=10)
results.GWAS_hb

##

#############################################################################
##

#############################################################################
##

#############################################################################

# again with less SNPS!

# hihi_confirmed_12022020_noZ
## turned into
# phen_gwaa_523_20Feb       ### previously the sex corrected

library(RepeatABEL)
# A model including polygenic effects,
# permanent environmental effects, and e.g. Year effect as random
Mod_1 <- preFitModel(FldgTarsus ~ sex + ClutchSize, random=~1|id + 1|Year + 1|Dam.ID + 1|Social.Sire.ID + 1|MonthYear,
                    genabel.data = phen_gwaa_523_20Feb, phenotype.data = PD,
                    corStruc=list( id=list("GRM") , Year=list("Ind") , 
                                   Dam.ID=list("Ind") , Social.Sire.ID=list("Ind") , MonthYear=list("Ind")) )
GWAS_1 <- rGLS(FldgTarsus ~ sex + ClutchSize, genabel.data = phen_gwaa_523_20Feb,
              phenotype.data = PD, V = Mod_1$V)   #If V is not specified then a model including random polygenic effects 
# and permanent environmental effects is fitted (using the hglm package) to compute V

# GWAS1@call$hglm = Mod2$fitted.hglm
est.hglm <- Mod_1$fitted.hglm
cat("Genotypic and permanent env. variance components:","\n",
    est.hglm$varRanef,", resp.","\n",
    "The residual variance is", est.hglm$varFix,".","\n")
##  computations for the genotypic variance
logVCE <- est.hglm$SummVC2[[1]]
cat("Estimate and SE for the genotypic variance on a natural log scale:",
    "\n", logVCE[1],"and",logVCE[2],", resp.", "\n \n",
    "Confidence interval: [", exp( logVCE[1] - 1.96*logVCE[2] ) , "," ,
    exp( logVCE[1] + 1.96*logVCE[2] ) , "]" , "\n")

## heritability:
cat("Heritability:",
    est.hglm$varRanef[1]/(est.hglm$varFix + sum(est.hglm$varRanef)),
    "\n")
## 0.081375    #without correcting for sire 0.1037163   ## just correcting for ID 0.3634722


lambda(GWAS_1) #1.036452                    #values up to 1.10 are generally considered acceptable for GWAS
estlambda(GWAS_1[, "P1df"], plot=TRUE)




############## FldgMass

Mod_2 <- preFitModel(FldgMass ~ sex + ClutchSize, random=~1|id + 1|Year + 1|Dam.ID + 1|Social.Sire.ID + 1|MonthYear,
                     genabel.data = phen_gwaa_523_20Feb, phenotype.data = PD,
                     corStruc=list( id=list("GRM") , Year=list("Ind") , 
                                    Dam.ID=list("Ind") , Social.Sire.ID=list("Ind") , MonthYear=list("Ind")) )
GWAS_2 <- rGLS(FldgMass ~ sex + ClutchSize, genabel.data = phen_gwaa_523_20Feb,
               phenotype.data = PD, V = Mod_2$V)   #If V is not specified then a model including random polygenic effects 
# and permanent environmental effects is fitted (using the hglm package) to compute V

# GWAS1@call$hglm = Mod2$fitted.hglm
est.hglm <- Mod_2$fitted.hglm
cat("Genotypic and permanent env. variance components:","\n",
    est.hglm$varRanef,", resp.","\n",
    "The residual variance is", est.hglm$varFix,".","\n")
##  computations for the genotypic variance
logVCE <- est.hglm$SummVC2[[1]]
cat("Estimate and SE for the genotypic variance on a natural log scale:",
    "\n", logVCE[1],"and",logVCE[2],", resp.", "\n \n",
    "Confidence interval: [", exp( logVCE[1] - 1.96*logVCE[2] ) , "," ,
    exp( logVCE[1] + 1.96*logVCE[2] ) , "]" , "\n")

## heritability:
cat("Heritability:",
    est.hglm$varRanef[1]/(est.hglm$varFix + sum(est.hglm$varRanef)),
    "\n")
## 0.06487538  #without correcting for sire: 0.09613173  ## just correcting for ID: 0.2306249


lambda(GWAS_2) #0.9827692
estlambda(GWAS_2[, "P1df"], plot=TRUE)

### FldgHead.bill

Mod_3 <- preFitModel(FldgHead.bill ~ sex + ClutchSize, random=~1|id + 1|Year + 1|Dam.ID + 1|Social.Sire.ID + 1|MonthYear,
                    genabel.data = phen_gwaa_523_20Feb, phenotype.data = PD,
                    corStruc=list( id=list("GRM"), Year=list("Ind") , 
                    Dam.ID=list("Ind") , Social.Sire.ID=list("Ind") , MonthYear=list("Ind")))
GWAS_3 <- rGLS(FldgHead.bill ~ sex + ClutchSize, genabel.data = phen_gwaa_523_20Feb,
              phenotype.data = PD, V = Mod_3$V)

# GWAS1@call$hglm = Mod2$fitted.hglm
est.hglm <- Mod_3$fitted.hglm
cat("Genotypic and permanent env. variance components:","\n",
    est.hglm$varRanef,", resp.","\n",
    "The residual variance is", est.hglm$varFix,".","\n")
##  computations for the genotypic variance
logVCE <- est.hglm$SummVC2[[1]]
cat("Estimate and SE for the genotypic variance on a natural log scale:",
    "\n", logVCE[1],"and",logVCE[2],", resp.", "\n \n",
    "Confidence interval: [", exp( logVCE[1] - 1.96*logVCE[2] ) , "," ,
    exp( logVCE[1] + 1.96*logVCE[2] ) , "]" , "\n")

## heritability:
cat("Heritability:",
    est.hglm$varRanef[1]/(est.hglm$varFix + sum(est.hglm$varRanef)),
    "\n")
## 0.07010987    #without correcting for sire 0.08852029      ## just correcting for ID  0.2435047

lambda(GWAS_3) #0.9572672
estlambda(GWAS_3[, "P1df"], plot=TRUE)


############### plot with a more relaxed association threshold
#tarsus
# library(ggplot2)

plot_1<-plot(GWAS_1, ylim = c(1,8), 
            main = "Tarsus length at 21 days",
            pch=20,
            col=c("blue4", "orange3"),
            space=0.2)
plot_1
abline(h= c(-log10(0.05/snps), -log10(1/snps), -log10(0.05)), lty= 2, col=c('black'))  #Bonferroni correction  1.259478e-06

results.GWAS_tl<-descriptives.scan(GWAS_1, top=10)
results.GWAS_tl
## as output table


estlambda(GWAS_1[, "P1df"], plot=TRUE, main = "Q-Q plot of GWAS p-values: Tarsus length")


##mass
plot_2<-plot(GWAS_2, ylim = c(1,8), 
            main = "Body mass at 21 days",
            pch=20,
            col=c("blue4", "orange3"))
plot_2
abline(h= c(-log10(0.05/snps), -log10(1/snps), -log10(0.05)), lty= 2, col=c('black'))  #Bonferroni correction  1.259478e-06

results.GWAS_bm<-descriptives.scan(GWAS_2, top=10)
results.GWAS_bm

estlambda(GWAS_2[, "P1df"], plot=TRUE, main = "Q-Q plot of GWAS p-values: Body mass")


##head.bill
plot_3<-plot(GWAS_3, ylim = c(1,8),
            main = "Head-bill length at 21 days",
            pch=20,
            col=c("blue4", "orange3"))
plot_3
abline(h= c(-log10(0.05/snps), -log10(1/snps), -log10(0.05)), lty= 2, col=c('black'))  #Bonferroni correction  1.259478e-06

results.GWAS_hb<-descriptives.scan(GWAS_3, top=10)
results.GWAS_hb

lambda(GWAS_3)
estlambda(GWAS_3[, "P1df"], plot=TRUE, main = "Q-Q plot of GWAS p-values: Head-bill length",
          proportion = 0.90) ### this was not done :/

## all in one plot
par(mfrow=c(3,1))

plot(GWAS_1, ylim = c(1,8),main="",
             pch=20,
             col=c("blue4", "orange3"))
title("(a) Tarsus length", adj = 0.05, line = -0.1)
abline(h= c(-log10(0.05/snps), -log10(1/snps), -log10(0.05)), lty= 2, col=c('black', 'darkred', 'black'))  #Bonferroni correction  1.259478e-06

plot(GWAS_2, ylim = c(1,8),main="",
     pch=20,
     col=c("blue4", "orange3"))
title("(b) Body mass", adj = 0.05, line = -0.1)
abline(h= c(-log10(0.05/snps), -log10(1/snps), -log10(0.05)), lty= 2, col=c('black', 'darkred', 'black'))  #Bonferroni correction  1.259478e-06

plot(GWAS_3, ylim = c(1,8),main="",
             pch=20,
             col=c("blue4", "orange3"))
title("(c) Head-bill length", adj = 0.05, line = -0.1)
abline(h= c(-log10(0.05/snps), -log10(1/snps), -log10(0.05)), lty= 2, col=c('black', 'darkred', 'black'))  #Bonferroni correction  1.259478e-06


## plot lmbda next to ech other
par(mfrow=c(1,3))
estlambda(GWAS_1[, "P1df"], plot=TRUE, main = "Q-Q plot of GWAS p-values: Tarsus length")
estlambda(GWAS_2[, "P1df"], plot=TRUE, main = "Q-Q plot of GWAS p-values: Body mass")
estlambda(GWAS_3[, "P1df"], plot=TRUE, main = "Q-Q plot of GWAS p-values: Head-bill length")

#### see how correlated the traits are
cor(c(PD$FldgTarsus), c(PD$FldgMass)) #0.6341046
cor(c(PD$FldgTarsus), c(PD$FldgHead.bill)) #0.6212057
cor(c(PD$FldgHead.bill), c(PD$FldgMass)) #0.7255911

#### see how correlated the significant SNPs are
cor(c(GWAS_1@results$P1df), c(GWAS_2@results$P1df)) #0.1349373 
cor(c(GWAS_1@results$P1df), c(GWAS_3@results$P1df)) #0.1627919
cor(c(GWAS_3@results$P1df), c(GWAS_2@results$P1df)) #0.3371452

################ in proper table
#https://www.r-bloggers.com/create-stylish-tables-in-r-using-formattable/

Top_10_SNPs_Tarsus <- results.GWAS_tl[,c(1:7,9)]
Top_10_SNPs_Mass <- results.GWAS_bm[,c(1:7,9)]
Top_10_SNPs_Head <- results.GWAS_hb[,c(1:7,9)]

library(data.table)
library(dplyr)
library(formattable)
library(tidyr)

#Set a few color variables to make our table more visually appealing
customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customRed = "#ff7f7f"

plo1 <- formattable(Top_10_SNPs_Tarsus)
plo2 <- formattable(Top_10_SNPs_Mass)
plo3 <- formattable(Top_10_SNPs_Head)

par(mfrow=c(3,1))
plo1
plo2
plo3

### add MAF
vec <- (summary(phen_gwaa@gtdata[,]))[,8:11]
head(vec)
## make rowname first row
library(data.table)
setDT(vec, keep.rownames = TRUE)[]
head(vec)

# reduce to top 10
top_ten_tarsus <- merge(vec, Top_10_SNPs_Tarsus, by = "rn", all.y = TRUE)

# reorder
setorder(top_ten_tarsus, P1df, na.last=TRUE)
top_ten_tarsus

# rename and reorder columns
names(top_ten_tarsus)[names(top_ten_tarsus)== "rn"] <- "SNP"
names(top_ten_tarsus)[names(top_ten_tarsus)== "Q.2"] <- "MAF"
names(top_ten_tarsus)[names(top_ten_tarsus)== "P.11"] <- "A1A1"
names(top_ten_tarsus)[names(top_ten_tarsus)== "P.12"] <- "A1A2"
names(top_ten_tarsus)[names(top_ten_tarsus)== "P.22"] <- "A2A2"
top_ten_tarsus

ttt <- top_ten_tarsus[,c(1,6,7,9,10,2,3,4,5,11,12,13)]
ttt <- as.data.frame(ttt)

## make table
library(gridExtra)
library(grid)
grid.table(ttt)

#Load the libraries

library(data.table)
library(dplyr)
library(formattable)
library(tidyr)
formattable(ttt, align = c("l", rep("r", NCOL(ttt) - 1)))


### body mass
# reduce to top 10
top_ten_mass <- merge(vec, Top_10_SNPs_Mass, by = "rn", all.y = TRUE)

# reorder
setorder(top_ten_mass, P1df, na.last=TRUE)
top_ten_mass

# rename and reorder columns
names(top_ten_mass)[names(top_ten_mass)== "rn"] <- "SNP"
names(top_ten_mass)[names(top_ten_mass)== "Q.2"] <- "MAF"
names(top_ten_mass)[names(top_ten_mass)== "P.11"] <- "A1A1"
names(top_ten_mass)[names(top_ten_mass)== "P.12"] <- "A1A2"
names(top_ten_mass)[names(top_ten_mass)== "P.22"] <- "A2A2"
top_ten_mass

ttm <- top_ten_mass[,c(1,6,7,9,10,2,3,4,5,11,12,13)]
ttm <- as.data.frame(ttm)

grid.table(ttm)

### head-bill length
# reduce to top 10
top_ten_head <- merge(vec, Top_10_SNPs_Head, by = "rn", all.y = TRUE)

# reorder
setorder(top_ten_head, P1df, na.last=TRUE)
top_ten_head

# rename and reorder columns
names(top_ten_head)[names(top_ten_head)== "rn"] <- "SNP"
names(top_ten_head)[names(top_ten_head)== "Q.2"] <- "MAF"
names(top_ten_head)[names(top_ten_head)== "P.11"] <- "A1A1"
names(top_ten_head)[names(top_ten_head)== "P.12"] <- "A1A2"
names(top_ten_head)[names(top_ten_head)== "P.22"] <- "A2A2"
top_ten_head

tth <- top_ten_head[,c(1,6,7,9,10,2,3,4,5,11,12,13)]
tth <- as.data.frame(tth)

grid.table(tth)

## change chr names
ttt$Chromosome <- as.character(ttt$Chromosome)
ttt$Chromosome[ttt$Chromosome=="31"] <- "1A"

tth$Chromosome <- as.character(tth$Chromosome)
tth$Chromosome[tth$Chromosome=="31"] <- "1A"

## formattable
write.table(ttt, file = "ttt.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(ttm, file = "ttm.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(tth, file = "tth.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)


### venn diagram

library(gplots)
venn(list(Top_10_SNPs_Tarsus[,2], Top_10_SNPs_Mass[,2], Top_10_SNPs_Head[,2]))

# Load library
library(VennDiagram)

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart simple
venn.diagram(
  x = list(Top_20_SNPs_Tarsus[,2], Top_20_SNPs_Mass[,2], Top_20_SNPs_Head[,2]),
  category.names = c("Tarsus length" , "Body mass" , "Head-bill length"),
  filename = '#March7_venn_diagramm.png',
  output=TRUE
)

# https://www.r-graph-gallery.com/14-venn-diagramm.html
library(tidyverse)
library(hrbrthemes)
library(tm)
library(proustr)

# Chart
venn.diagram(
  x = list(Top_20_SNPs_Tarsus[,2], Top_20_SNPs_Mass[,2], Top_20_SNPs_Head[,2]),
  category.names = c("Tarsus length" , "Body mass" , "Head-bill length"),
  filename = '#March7_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  rotation = 1
)

