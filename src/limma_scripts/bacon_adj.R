bacon_adj <- function(table, adj_method){
lambda <- median(qchisq(table$P.Value, df=1, lower.tail=F))/qchisq(0.5, df=1)
#lambda <- 1.5
if (lambda>1.2){
  bc <- bacon(teststatistics=table$t)
  ps <- as.data.frame(pval(bc))
  table$Bacon_P <- ps$V1
  table$Bacon_FDR <- p.adjust(ps$V1, method=adj_method)
  table <- merge(table, annot, by = "row.names")
  table <- table %>% tibble::as_tibble() %>% dplyr::arrange(Bacon_FDR)
  lambda <- list(pre_Bacon=lambda, post_Bacon=median(qchisq(table$Bacon_P, df=1, lower.tail=F))/qchisq(0.5, df=1))
}else{
  table <- merge(table, annot, by = "row.names")
  table <- table %>% tibble::as_tibble() %>% dplyr::arrange(adj.P.Val)
}
return(list(results_table = table, lambdas = lambda))
}