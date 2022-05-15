# missMethyl
mm_it <- function(table, missmethyl_path, var, var_extra){
if("Bacon_FDR" %in% colnames(table) & any(table$Bacon_FDR <0.25)){
message("Conducting missmethyl on baconed results")
sig_cpgs <- table %>%
            filter(Bacon_FDR <0.25) %>%
            pull(row.names) %>%
            as.character()
# all the CpGs you want to analyse (i.e. top 100 / FDR <0.25) as a character vector - remember to arrange by FDR / bacon_FDR if you want to slice off the top 100!
all_cpgs <- table %>%
            pull(row.names) %>%
            as.character() # every CpG in the analysis as a character vector

g <- missMethyl::gometh(sig.cpg=sig_cpgs, all.cpg=all_cpgs, collection=c("GO","KEGG"), array.type="EPIC") # choose GO or KEGG
go_results <- topGSA(g, n=Inf)
saveRDS(go_results, paste0(missmethyl_path, var, "_", var_extra, ".rds"))

} else if(!"Bacon_FDR" %in% colnames(table) & any(table$adj.P.Val <0.25)) {
  message("Conducting missmethyl on results")
  sig_cpgs <- table %>%
              filter(adj.P.Val <0.25) %>%
              pull(row.names) %>%
              as.character()
  # all the CpGs you want to analyse (i.e. top 100 / FDR <0.25) as a character vector - remember to arrange by FDR / bacon_FDR if you want to slice off the top 100!
  all_cpgs <- table %>%
              pull(row.names) %>%
              as.character() # every CpG in the analysis as a character vector

  g <- missMethyl::gometh(sig.cpg=sig_cpgs, all.cpg=all_cpgs, collection=c("GO","KEGG"), array.type="EPIC") # choose GO or KEGG
  go_results <- topGSA(g, n=Inf)
  saveRDS(go_results, paste0(missmethyl_path, var, "_", var_extra, ".rds"))
}
}
