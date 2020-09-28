[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3517785.svg)](https://doi.org/10.5281/zenodo.3517785)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/markkukuismin/BONE)
![GitHub](https://img.shields.io/github/license/markkukuismin/BONE)

# BONE: Baseline Oriented Network Estimation

Genetic assignment of individuals to known source populations using network estimation tools.

# Install package

The package can be loaded from GitHub. If you are not using RStudio start now!

I used RStudio (version 1.1.453) and Microsoft R Open (version 3.5.1.) to create this package.

Unzip the "BONE.zip" file into a working directory and run the following lines:

```r
library(devtools)
library(roxygen2)

install("BONE")
```
Now you should be able to use all functions and read their descriptions

```r
library(BONE)
library(glmnet)
```

Alternatively you can use functions of BONE using the source files found in the folder "BONE\R" after you have extracted the "BONE.zip" file:

```r
source("impute.R")
source("LASSOSolPath.R")
source("SolPathInference.R")
source("WTAInference.R")
```

BONE depends on the R package "glmnet".

# Example data

In the Data.zip file there are two data tables: a mixture table and a reference table.

The data is a collection of island and mainland populations of house sparrows along the coast of Norway collected in year 2012. Mixture individuals were sampled from the original reference data using the "mixture_draw" function of rubias package (see Moran and Anderson, "Bayesian inference from the conditional genetic stock identification model", Canadian Journal of Fisheries and Aquatic Sciences, 2019, 76:551-560, https://doi.org/10.1139/cjfas-2018-0016).

Unzip files into working directory.

```r
Classes = c(rep("character",4),rep("integer",2000))

MixtureData = read.table("MixtureData.txt",colClasses = Classes,header=T)

ReferenceData = read.table("ReferenceData.txt",colClasses = Classes,header=T)
```

This package works initially with RUBIAS (dplyr data frame) file format. Here is an example of four mixture individuals and one locus:

```
 sample_type repunit collection    indiv V315 V316
     mixture      20       20.1 a353e4fb    4    4
     mixture      20       20.1 5cad80d8    4    2
     mixture      22       22.1 c0fe50a9    4    2
     mixture      22       22.1 27467dae    4    2
```

# Some illustrative experiments

Compute probability of the origin and mixture proportions using the solution path approach:

```r
> str(MixtureData[, 1:6])

'data.frame':	99 obs. of  6 variables:
 $ sample_type: chr  "mixture" "mixture" "mixture" "mixture" ...
 $ repunit    : chr  "20" "20" "22" "22" ...
 $ collection : chr  "20.1" "20.1" "22.1" "22.1" ...
 $ indiv      : chr  "a353e4fb" "5cad80d8" "c0fe50a9" "27467dae" ...
 $ V315       : int  4 4 4 4 2 4 2 4 4 4 ...
 $ V316       : int  4 2 2 2 2 2 2 2 2 2 ...
 
 > str(ReferenceData[, 1:6])
 
'data.frame':	469 obs. of  6 variables:
 $ sample_type: chr  "reference" "reference" "reference" "reference" ...
 $ repunit    : chr  "38" "38" "38" "38" ...
 $ collection : chr  "38.1" "38.1" "38.1" "38.1" ...
 $ indiv      : chr  "ae7297ac" "a1b0350f" "300c3c26" "3bc4f361" ...
 $ V315       : int  4 4 4 4 4 4 4 4 4 4 ...
 $ V316       : int  4 2 4 2 4 4 4 2 2 2 ...
```

Data premodification (compute genotype matrix, impute missing values) are done using the ```r impute ``` function. Run it although you would not need to impute your data. User can also define the genotyping method. If the genotyping method is set as "Fluidigm", the original allele coding is preserved and four "genotypes" are returned (alleles are divided into independent observations at each locus). If missing genotypes are not imputed, missing values are denoted with "18" (Axiom) or "9" (Fluidigm).

For example:

```r
MixtureY = impute(MixtureData, genotyping = "Fluidigm")
MixtureY = MixtureY$Y
```
```
> MixtureY[1:5,1:4]
      a353e4fb 5cad80d8 c0fe50a9 27467dae
V315         4        4        4        4
V316         4        2        2        2
V681         1        1        3        1
V682         3        3        3        1
V1047        4        4        4        2
```
The default genotyping method is Axiom.

```r
MixtureY = impute(MixtureData) # Simple marker mode imputation
ReferenceY = impute(ReferenceData)

MixtureY = MixtureY$Y
ReferenceY = ReferenceY$Y

# Choose only same set of markers:
setdiff(rownames(MixtureY),rownames(ReferenceY))
[1] "V259281"

Markers = intersect(rownames(MixtureY),rownames(ReferenceY))

MixtureY = MixtureY[Markers,]
ReferenceY = ReferenceY[Markers,]

SampleSizeBaselinePop = ncol(ReferenceY)
Y = cbind(ReferenceY,MixtureY)

lambda = seq(0.4,0.02,length.out=40)

MBapprox = LASSOSolPath(Y,lambda,intercept=T,SampleSizeBaselinePop,Baseline=T)

SampleNames=colnames(Y)

NetworkResultsSolpath = SolPathInference(MBapprox,ReferenceData,SampleNames=SampleNames,
                                         SampleSizeBaselinePop,alpha=0.05)
```
or using the "Winner Takes it All" approach:

```r
NetworkResultsWTA = WTAInference(MBapprox,ReferenceData,SampleNames=SampleNames,SampleSizeBaselinePop)
```
```
> NetworkResultsWTA$ProbofOrigin[1:5,]
         20 22 23 24 26 27 28 38 77
a353e4fb  0  0  0  0  1  0  0  0  0
5cad80d8  1  0  0  0  0  0  0  0  0
c0fe50a9  0  1  0  0  0  0  0  0  0
27467dae  0  1  0  0  0  0  0  0  0
3449e4f9  0  0  0  1  0  0  0  0  0
```
```
> round(NetworkResultsSolpath$ProbofOrigin[1:5,],4)
             20     22     23     24     26     27     28     38     77 DifferenceProportion
a353e4fb 0.0511 0.0125 0.0404 0.0351 0.5490 0.2739 0.0036 0.0006 0.0339               0.0125
5cad80d8 0.4720 0.0931 0.0240 0.0038 0.1705 0.2072 0.0004 0.0038 0.0252               0.0267
c0fe50a9 0.0047 0.2933 0.1052 0.1296 0.1225 0.1590 0.0561 0.0014 0.1282               0.0103
27467dae 0.0165 0.3520 0.0787 0.1882 0.2126 0.0819 0.0528 0.0000 0.0173               0.0156
3449e4f9 0.0303 0.3088 0.1769 0.3672 0.0269 0.0252 0.0315 0.0000 0.0332               0.0277
```
```r
> NetworkResultsWTA$MixtureProp
   Pop\tEst.Pi
20 0.03030303
22 0.03030303
23 0.07070707
24 0.11111111
26 0.29292929
27 0.23737374
28 0.10101010
38 0.03030303
77 0.09595960

> NetworkResultsSolpath$MixtureProp
   Pop\tEst.Pi
20 0.02353646
22 0.02261854
23 0.08223077
24 0.11292987
26 0.29707736
27 0.22732527
28 0.11484434
38 0.03478313
77 0.08465427

```
"DifferenceProportion" is the MSE between probability of the origin estimated with the solution path method and expected probability of the origin.

Probability of the origin and mixture proportion estimates are also saved into an external file "NetworkResults" in the working directory.

# Reference

The BONE method is described in:

Kuismin et al. (in press) "Genetic assignment of individuals to known source populations using network estimation tools". *Methods in Ecology and Evolution*.

File "CodeCollection.zip" is a collection of scripts and data sets used to prepare the material in this paper.
