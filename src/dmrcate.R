EPICdmrcate <- function(combat_beta,
                        design,
                        annotation = c(array = "IlluminaHumanMethylationEPIC",
                                       annotation ="ilm10b4.hg19"),
                        analysis.type="differential",
                        coef="variable",
                        fdr=0.05,
                        lambda=500,
                        C=2,
                        pcutoff=fdr,
                        min.cpgs=2,...){
	require("DMRcate")
	dm <- cpg.annotate(datatype="array",
                     object=combat_beta,
                     what="Beta",
                     annotation=annotation,
                     design=design,
                     arraytype="EPIC",
                     analysis.type=analysis.type,
                     fdr=fdr,
                     coef=coef,
                     method="robust",...)
	DMRs <- dmrcate(dm,
                  lambda=lambda,
                  C=C,
                  pcutoff=pcutoff,
                  min.cpgs=min.cpgs)
	return(list(DMRs=DMRs, cpg_object=dm))
}

# # Generate your model matrix
# ## Fill in your variable of interest (var), extra info and covariates
# var <- "jsex"
# covariates <- c(#"jsex",
#                 "position_on_array",
#                 "CD4T",
#                 "CD8T",
#                 "Bcell",
#                 "Gran",
#                 "Mono",
#                 "NK",
#                 "nRBC")
#
# phenotypic_data <- p_data
# p <- phenotypic_data %>%
#      drop_na(any_of(append(covariates, var)))
#
# combat_beta_sub <- combat_beta[,colnames(combat_beta) %in% p[,1]] # Make sure the p$ is the column name for your samples = here it is swsid
#
# mod <- model.matrix(as.formula(paste("~", var, "+", paste(covariates, collapse="+"), sep="")), data=p)
#
# #filtered_beta_10 <- combat_beta_sub
#
# # Run the analysis:
# dmroutput <- EPICdmrcate(combat_beta, mod, coef=2, fdr= 0.05, pcutoff= 0.05)
# results.ranges <- extractRanges(dmroutput, genome="hg19")
## mod is the model matrix that you use for the analysis
## coef is the column index of the variable of interest in your model matrix
