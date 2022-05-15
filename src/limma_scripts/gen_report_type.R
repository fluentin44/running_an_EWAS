choose_report <- function(table, lambda){
report_template <- vector("list")
#if("P.Value" %in% colnames(table) == TRUE & length(which(str_detect(rownames(table),glob2rx("cg*")))) > 1L){
#if(lambda > 1.2){
if("Bacon_FDR" %in% colnames(table) == TRUE){
message("Results required BACON adjustment")
#report_template$table <- rownames_to_column(table, var="row.names")
report_template$table <- table
names(report_template$table)<-  c("row.names",
                                 "logFC",
                                 "AveExpr",
                                 "t",
                                 "P.Value",
                                 "adj.P.Val",
                                 "B",
                                 "Bacon_P",
                                 "Bacon_FDR",
                                 "ucsc_refgene_name",
                                 "ucsc_refgene_group",
                                 "relation_to_island",
                                 "chr",
                                 "pos",
                                 "strand",
                                 "name",
                                 "regulatory_feature_group",
                                 "ucsc_refgene_accession"
                                 )
report_template$report_type <- "regression_report_limma_bacon"
} else{
message("Results did not require BACON adjustment")
report_template$table <- table
names(report_template$table)<-  c("row.names",
                                  "logFC",
                                  "AveExpr",
                                  "t",
                                  "P.Value",
                                  "adj.P.Val",
                                  "B",
                                  "ucsc_refgene_name",
                                  "ucsc_refgene_group",
                                  "relation_to_island",
                                  "chr",
                                  "pos",
                                  "strand",
                                  "name",
                                  "regulatory_feature_group",
                                  "ucsc_refgene_accession"
                                  )
report_template$report_type <- "regression_report_limma"
}
return(report_template)
}
