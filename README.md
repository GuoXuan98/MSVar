---
Title: MSVar 1.0.0
Author: Xuan Guo
Date: "2025-12-2"
Contact: guoxuan2020@sinh.ac.cn
---


# Introduction
MSVar is designed for isobaric labeling-based mass spectrometry (ILMS) 
quantitative proteomic profiling to identify hypervariable proteins for a single 
group of samples and differentially variable proteins between different groups
via similar procedures. The latest version of MSVar is always available in the 
[CRAN](https://cran.r-project.org) repository and can thus be easily installed
by typing `install.packages("MSVar")` in an R session.

For other versions of MSVar, select a package from under the `dist` folder,
download it, and type `install.packages("/path/to/the/package", repos = NULL)`.
In this way, you may need to pre-install some dependencies of MSVar. The current
dependencies of the latest MSVar version include locfit (>= 1.5.9), 
scales (>= 0.3.0), and statmod (>= 1.4.34). All these packages are available in
the [CRAN](https://cran.r-project.org) repository. For dependencies of other
MSVar versions, refer to the `Imports` field in the corresponding `DESCRIPTION`
file.

# Format of Input Data

For employing the machinery implemented in MSVar, you need to prepare several
tables that profiles the protein intensity. Large-scale proteomic data sets are
commonly quantified using tandem mass tagging (TMT) method. For TMT-based 
quantitative proteomic data, the following paired two tables provide a example:

Table 1 (raw.intensity):
| protein | Set1_MIX |   T112 |   T113 | Set2_MIX |    T131 |    T135 |
| ------: | -------: | -----: | -----: | -------: | ------: | ------: |
|    A1BG |   197220 | 413850 | 243030 |   227520 |  222190 |  235810 |
|    A1CF |   432210 | 443780 | 421840 |   542740 |  892760 |  590000 |
|     A2M |   939990 | 961500 | 848020 |   737220 | 1029600 | 1432900 |
|    AAAS |   115830 | 110230 |  93278 |   119680 |  107080 |   88204 |

Table 2 (batch.info):
| reference | sample1 | sample2 |
| --------: | ------: | ------: |
|  Set1_MIX |    T112 |    T113 |
|  Set2_MIX |    T131 |    T135 |

Table 1 records raw protein abundance, each row represents a protein and each 
column corresponds to a reference or tissue sample. Table 2 documents the 
correspondence between reference samples and tissue samples.

After preparing these two tables, you can use `proObjFromTMT` function to create
a "proObj" object.

```r
library(MSVar)
proObj <- proObjFromTMT(raw.intensity, batch.info)
```

If the raw abundance data for reference and tissue samples are stored in two
separate tables, you can directly create a "proObj" object using the `proObj`
function. The key step is to specify the correspondence between the reference 
and tissue samples.

```r
library(MSVar)
?proObjFromTMT
?proObj
```

# Technical Variance Estimation

MSVar decomposes the signal variability across a group of ILMS samples into 
technical and biological parts (i.e. $\tau_{ij}^2$ and $\sigma_{ij}^2$), 
corresponding respectively to the technical noise of ILMS experiments and the 
variability of true protein expression in the population represented by the 
sample group.

Technically, it derive estimates of all $\tau_{ij}^2$ by applying the 
sliding-window procedure originally designed by
[MAP](https://doi.org/10.1038/s41421-019-0107-9) and 
[zMAP](https://doi.org/10.1186/s13059-024-03382-9), which exploits the 
mean-variance dependence of ILMS data.

The estimation result can be easily obtained by calling `estTechVar` function:

```r
library(MSVar)
ProObj <- estTechVar(ProObj)
```

This procedure is for estimating the technical noise from MS quantification for 
each protein in each sample. For more detailed parameter settings regarding the 
estimation method implemented in MSVar, type the following code in an R session 
after installing it:

```r
library(MSVar)
?estTechVar
```

# Biological Variance Estimation

With estimated technical variance $\tau_{ij}^2$, MSVar next estimate expression 
mean $m_i$ and technical variance $\sigma_i^2$ separately for each protein by 
applying maximum likelihood estimation to all the observed M-values of it. Then,
the uncertainty of these estimates, denoted by $\hat{m}_i$ and 
$\hat{\sigma}_i^2$, are approximated by using Fisher information matrix.

These estimation results can be gained by calling `estBioVar` function:

```r
library(MSVar)
ProObj <- estBioVar(ProObj)
```

This procedure provides an estimation of $m_i$ and $\sigma_i^2$ for each 
protein. More detailed parameter settings are accessible by typing the following
code in an R session after installing MSVar:

```r
library(MSVar)
?estTechVar
```

# Posterior M-value Derivation

The Bayesian nature of MSVar allows the deduction of posterior estimation of 
true protein expression. A posterior estimate of each $M^\prime_{ij}$ can be
derived by finding the maximum a *posteriori* (MAP) of $M^\prime_{ij}|M_{ij}$.

The posterior estimate of each $M^\prime_{ij}$ can be acquired by calling `PostM` 
function:

```r
library(MSVar)
ProObj <- PostM(ProObj)
```

# Hypervariable Analyses

MSVar identifies hypervariable proteins (HVPs) for a single group of samples by comparing the biological variance estimates $\sigma_i^2$ of all proteins to a pre-specified constant, and the above applications of MSVar were all based on 
the default setting of this constant, which amounted to an average fluctuation 
of 20% around mean expression. 

The following code conducts hypothesis testing procedures for identifying HVPs
under default parameter settings:

```r
library(MSVar)
HVP <- VarTest(ProObj)
```

For more information regarding the testing procedure implemented in MSVar, type
the following code in an R session after installing it:

```r
library(MSVar)
?VarTest
```

# Differentially Variable Analyses

The MSVar tests for identifying differentially variable proteins (DVPs) between
different sample groups it is designed to compare the corresponding biological
variance estimates at log scale.

For identifying DVPs between two sample groups (say group 1 and 2), MSVar first 
fits the single-group model separately for each group and create 2 "proObj" 
objects, namely ProObj1 and ProObj2. An additional index is added to the 
associated parameters to indicate their groups of origin. For example, 
$v_{i,1}$ represents log-biological variance of protein $i$ in group 1.

The following code conducts hypothesis testing procedures for identifying DVPs 
with higher variability in group 1 under default parameter settings:

```r
library(MSVar)
DVP <- DiffVar(ProObj1, ProObj2)
```

For more information regarding the testing procedure implemented in MSVar, type
the following code in an R session after installing it:

```r
library(MSVar)
?DiffVar
```
