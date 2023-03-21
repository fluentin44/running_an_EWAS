run.ewas <- function(var,
                     var_extra,
                     combat_data,
                     p,
                     results_path=limma_save_path,
                     covariates=NULL,
                     test=F,
                     sva=F,
                     first_run=T,
                     no_svs=NULL
                     ){

if(!sva){

message(paste0("starting analysis of ", var, ", ", var_extra))
message(paste0("Non-SVA first run"))

p <-
  p %>%
  drop_na(any_of(append(covariates, var)))

if(!grep("id", colnames(p), ignore.case=T) == 1){
  stop("First column is not ID or column does not exist")
}

combat_beta_sub <- combat_data[,colnames(combat_data) %in% p[,1]] # Make sure the p$ is the column name for your samples = here it is swsid

if(nrow(p)>ncol(combat_beta_sub)){
  p <- p[p[,1] %in% colnames(combat_beta_sub), ]
}

#ensures that the order of samples (SWS IDs) in the combat beta match the order in the phenotypic file
message(summary(colnames(combat_beta_sub)==p[,1])) #check order and return True/False
combat_beta_sub <- combat_beta_sub[,match(p[,1], colnames(combat_beta_sub))] #re-orders samples
message(summary(colnames(combat_beta_sub)==p[,1])) #check order and return True/False, should now be True

if(test) {
  combat <-
    combat_beta_sub[1:10000, ]
  message("Test run")
} else{
  combat <-
    combat_beta_sub
  message("Full run")
}

message(paste0("Dimension of combat beta sub = ", nrow(combat)))
message(paste0("number of non-na samples for ", var," ", var_extra, " = ", nrow(p)))

if(!is.null(covariates)){
  message("covariates exist")
  mod <- model.matrix(as.formula(paste("~", var, "+", paste(covariates, collapse="+"), sep="")), data=p)
  message(paste("Model = ~", var, "+", paste(covariates, collapse="+"), sep=""))
} else{
  message("covariates dont exist")
  mod <- model.matrix(as.formula(paste0("~", var)), data=p)
  message(paste0("Model = ~", var))
}

print(head(mod))

set.seed(12345)
fit <- lmFit(combat, mod, method = "robust") # Use this line for full run, hash out for testing
fit_2 <- eBayes(fit)
table <- 
limma::topTable(fit_2, 
coef = 2, 
number = Inf,
confint=TRUE)

analysis <- bacon_adj(table, adj_method = "BH")
table <- analysis$results_table
lambda <- analysis$lambdas
covariates <- paste0(covariates, collapse = ",")

table %>%
  print() %>%
  readr::write_tsv(
    results_path,
    na = "",
    col_names = TRUE
  )

unaltered_table <- table
report <- choose_report(table, lambda)
table <- report$table
report_type <- report$report_type

# Summary table
result_list <- regression_summary(var,
                                  var_extra,
                                  covriates,
                                  p,
                                  lambda,
                                  table,
                                  lambdas_path,
                                  report_type)

# Render report - parameters are classified by whats defined in function
gen_reportee(
  voi                        = var,
  raw_results                = table,
  report_template_filename   = report_type,
  output_file_name           = output_file_name,
  covariate_names            = covariates,
  results_summary            = result_list,
  date_started               = date_run_started,
  p_data                     = p
)

# vivs<-list(table = table, report_type = report_type)
print(mod)

} else if(sva & first_run) {

  message(paste0("Getting SVA values for the analysis of ", var, ", ", var_extra))
  message(paste0("SVA first run"))

  p <-
    p %>%
    drop_na(any_of(append(covariates, var)))

  combat_beta_sub <- combat_data[,colnames(combat_data) %in% p[,1]] # Make sure the p$ is the column name for your samples = here it is swsid

  if(nrow(p)>ncol(combat_beta_sub)){
    p <- p[p[,1] %in% colnames(combat_beta_sub), ]
  }

  #ensures that the order of samples (SWS IDs) in the combat beta match the order in the phenotypic file
  print(summary(colnames(combat_beta_sub)==p[,1])) #check order and return True/False
  combat_beta_sub <- combat_beta_sub[,match(p[,1], colnames(combat_beta_sub))] #re-orders samples
  print(summary(colnames(combat_beta_sub)==p[,1])) #check order and return True/False, should now be True

  if(test) {
    combat <-
      combat_beta_sub[1:10000, ]
    message("Test run")
  } else{
    combat <-
      combat_beta_sub
    message("Full run")
  }

  message(paste0("Dimension of combat beta sub = ", nrow(combat)))
  message(paste0("number of non-na samples for ", var," ", var_extra, " = ", nrow(p)))


  mod <- model.matrix(as.formula(paste("~", var, "+", paste(covariates, collapse="+"), sep="")), data=p)
  mod0 <- model.matrix(as.formula(paste("~", paste(covariates, collapse="+"), sep="")), data=p)

  set.seed(12345)

  sv_n= sva(combat, # Use this block for full_run, hash out for testing
            mod=mod,
            mod0=mod0,
            controls = NULL,
            method = "irw",
            vfilter = NULL,
            numSVmethod = "be")

  print(paste0("number of surrogate variables = ", sv_n$n.sv))

} else if(sva & !first_run) {

  message(paste0("starting the SVA analysis of ", var, ", ", var_extra))
  message(paste0("SVA non-first run"))

  p <-
    p %>%
    drop_na(any_of(append(covariates, var)))

  combat_beta_sub <- combat_data[,colnames(combat_data) %in% p[,1]] # Make sure the p$ is the column name for your samples = here it is swsid

  if(nrow(p)>ncol(combat_beta_sub)){
    p <- p[p[,1] %in% colnames(combat_beta_sub), ]
  }

  #ensures that the order of samples (SWS IDs) in the combat beta match the order in the phenotypic file
  print(summary(colnames(combat_beta_sub)==p[,1])) #check order and return True/False
  combat_beta_sub <- combat_beta_sub[,match(p[,1], colnames(combat_beta_sub))] #re-orders samples
  print(summary(colnames(combat_beta_sub)==p[,1])) #check order and return True/False, should now be True

  if(test) {
    combat <-
      combat_beta_sub[1:10000, ]
    message("Test run")
  } else{
    combat <-
      combat_beta_sub
    message("Full run")
  }

  message(paste0("Dimension of combat beta sub = ", nrow(combat)))
  message(paste0("number of non-na samples for ", var," ", var_extra, " = ", nrow(p)))

  mod <- model.matrix(as.formula(paste("~", var, "+", paste(covariates, collapse="+"), sep="")), data=p)
  mod0 <- model.matrix(as.formula(paste("~", paste(covariates, collapse="+"), sep="")), data=p)

  set.seed(12345)

  sv_n= sva(combat, # Use this block for full_run, hash out for testing
            mod=mod,
            mod0=mod0,
            controls = NULL,
            method = "irw",
            vfilter = NULL,
            numSVmethod = "be")

  print(paste0("Total number of calculated surrogate variables = ", sv_n$n.sv))

  if(is.null(no_svs) & sv_n$n.sv > 0){
    message(paste0("Integrating all ", ncol(sv_n$sv), " surrogate variables"))
    mod1 <- cbind(mod, sv_n$sv[,1:ncol(sv_n$sv)])

    og <- colnames(mod1)[!colnames(mod1) == ""]
    bab <- append(paste0(og), paste0("SV", 1:ncol(sv_n$sv)))
    colnames(mod1) <- bab

    covariates <- append(paste0(covariates, collapse = ","), paste0("SV", 1:ncol(sv_n$sv)))
} else if(no_svs > 0 & sv_n$n.sv > 0){
    message(paste0("Integrating ", no_svs, " surrogate variables"))
    mod1 <- cbind(mod, sv_n$sv[,1:no_svs])

    og <- colnames(mod1)[!colnames(mod1) == ""]
    bab <- append(paste0(og), paste0("SV", 1:no_svs))
    colnames(mod1) <- bab

    covariates <- append(paste0(covariates, collapse = ","), paste0("SV", 1:no_svs))
  }


#  og <- colnames(mod1)[!colnames(mod1) == ""]
#  bab <- append(paste0(og), paste0("SV", 1:no_svs))
#  colnames(mod1) <- bab

  print(colnames(mod1))

  #print(paste0("number of surrogate variables in model = ", sum(colnames(mod1) == "", na.rm=TRUE)))

  set.seed(12345)
  fit <- lmFit(combat, mod1, method = "robust") # Use this line for full run, hash out for testing
  fit_2 <- eBayes(fit)
  table <- 
  limma::topTable(fit_2, 
  coef = 2, 
  number = Inf, 
  confint=TRUE)

  analysis <- bacon_adj(table, adj_method = "BH")
  table <- analysis$results_table
  lambda <- analysis$lambdas

  table %>%
    print() %>%
    readr::write_tsv(
      results_path,
      na = "",
      col_names = TRUE
    )
  unaltered_table <- table 
  report <- choose_report(table, lambda)
  table <- report$table
  report_type <- report$report_type

  # Summary table
  result_list <- regression_summary(var,
                                    var_extra,
                                    covriates,
                                    p,
                                    lambda,
                                    table,
                                    lambdas_path,
                                    report_type)

  # Render report - parameters are classified by whats defined in function
  gen_reportee(
    voi                        = var,
    raw_results                = table,
    report_template_filename   = report_type,
    output_file_name           = output_file_name,
    covariate_names            = covariates,
    results_summary            = result_list,
    date_started               = date_run_started,
    p_data                     = p
  )
  # vivs<-list(table = table, report_type = report_type)
  mod <- mod1
  #print(mod)
}
#return(mod)
return(list(mod=mod, 
            p=p, 
            combat=combat, 
            table=table,
            unaltered_table=unaltered_table, 
            lambda=lambda, 
            fit=fit, 
            fit_2=fit_2,
            var=var,
            table=table,
            report_type=report_type,
            output_file_name=output_file_name,
            covariates=covariates,
            result_list=result_list,
            date_run_started=date_run_started))
}
