regression_summary <- function(var, var_extra, covriates, p, lambda, table, lambdas_path, report_type){
if(report_type == "regression_report_limma"){
  result_list <- tibble(
  "variable" = paste0(var),
  "augment" = paste0(var_extra),
  "covs" = paste0(covariates, collapse = ","),
  "samples" = nrow(p),
  "Lambda" = paste0(signif(lambda,3)),
  "FDR<0.25" = sum(table$adj.P.Val < 0.25),
  "FDR<0.20" = sum(table$adj.P.Val < 0.20),
  "FDR<0.15" = sum(table$adj.P.Val < 0.15),
  "FDR<0.10" = sum(table$adj.P.Val < 0.10),
  "FDR<0.05" = sum(table$adj.P.Val < 0.05)
  ) %>%
  print(width = Inf) %>%
  readr::write_tsv(lambdas_path,
  na = "",
  col_names = TRUE,
  append = TRUE)
} else{
  result_list <- tibble(
  "variable" = paste0(var),
  "augment" = paste0(var_extra),
  "covs" = paste0(covariates, collapse = ","),
  "bacon" = c("Pre-bacon", "Post-bacon"),
  "samples" = nrow(p),
  "Lambda" = c(paste0(lambda$pre_Bacon),paste0(signif(lambda$post_Bacon, 3))),
  "FDR<0.25" = c(sum(table$adj.P.Val < 0.25),sum(table$Bacon_FDR < 0.25)),
  "FDR<0.20" = c(sum(table$adj.P.Val < 0.20),sum(table$Bacon_FDR < 0.20)),
  "FDR<0.15" = c(sum(table$adj.P.Val < 0.15),sum(table$Bacon_FDR < 0.15)),
  "FDR<0.10" = c(sum(table$adj.P.Val < 0.10),sum(table$Bacon_FDR < 0.10)),
  "FDR<0.05" = c(sum(table$adj.P.Val < 0.05),sum(table$Bacon_FDR < 0.05)),
  ) %>%
  print(width = Inf) %>%
  readr::write_tsv(lambdas_path,
  na = "",
  col_names = TRUE,
  append = TRUE
  )
}
}
