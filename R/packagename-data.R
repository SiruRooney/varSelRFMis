#'
#' A list contains the simulated genotype data with missing values, phenotype data with missing values
#'
#'@docType data
#'@name gene_phe
#'@format a list
#'
#'
NULL
load("E:/varselRFmisnonign210818/simulationgeneration/Simulation_genotypeSNP_analyplusmissing3.RData")
gene_phe=list(SNPsim=SNPsim,phe=phe,phemis=phemis)
rm(SNPsim)
rm(phemis)
rm(phe)

