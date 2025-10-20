# varSelRFMis
<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.1016/j.apm.2023.06.025.svg)](https://doi.org/10.1016/j.apm.2023.06.025)
<!-- badges: end -->


This is an R package for genes selection using RF with missing values under non-ignorable missingness mechanism.
Website available:<https://SiruRooney.github.io/varSelRFMis/>

For more details, please see:
Wang, S., Qian, G., & Hopper, J. (2023). Integrated logistic ridge regression and random forest for phenotype-genotype association analysis in categorical genomic data containing non-ignorable missing values. Applied Mathematical Modelling, 123, 1-22. [PDF](https://www.sciencedirect.com/science/article/pii/S0307904X23002809)


```{r}
library("devtools");

install_github("SiruRooney/varSelRFMis");

library("varSelRFMis");

SNPsim=varSelRFMis::gene_phe$SNPsim;

phemis=varSelRFMis::gene_phe$phemis;

dectna=detctnafun(SNPsim,phemis);

library(doParallel);

result=varselrfmisfun(SNPsim,phemis,dectna$misind,dectna$detna,iter.varsel=10,initmethod="init",optmethod="cv",dummytran,vsrfmisfun,jimpmisfun,mulimpmisfun,       imputerfmisfun,varSelRF_mis,gensetf,bicglmnet,bdiag_m);
```

