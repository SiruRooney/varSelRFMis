# varSelRFMis
Genes selection using RF with missing values

library("devtools")
install_github("SiruRooney/varSelRFMis")
library("varSelRFMis")

SNPsim=varSelRFMis::gene_phe$SNPsim
phemis=varSelRFMis::gene_phe$phemis

dectna=detctnafun(SNPsim,phemis)
library(doParallel)
result=varselrfmisfun(SNPsim,phemis,dectna$misind,dectna$detna,iter.varsel=10,initmethod="init",optmethod="cv",dummytran,vsrfmisfun,jimpmisfun,mulimpmisfun,       imputerfmisfun,varSelRF_mis,gensetf,bicglmnet,bdiag_m)


