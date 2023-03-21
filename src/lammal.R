# pheno_data: Data frame with all phenotype data and covariates used in the analysis, with the first column as sample IDs
# meth_matrix: Data frame of beta values, with rows as probes and columns as samples. Column names need to be sample IDs
# adj_method: Method for adjustment for multiple testing (e.g. BH, bonferroni, FDR), passed onto p.adjust()
# outcome and predictor: One of these should be 'methylation', the other should be the variable of interest, depending on the analysis required
# plots: logical for if you want regression diagnostics plots saved for top 5 dmCpGs (will be saved to working directory)
# seed: Number to begin random number generation for replicability, passed onto set.seed()
# coef: If your predictor is a factor and does not equal "methylation", coef needs to specifiy the level of the factor that you want results for (cannot=1 (reference) level)

suppressMessages(library(robust))
suppressMessages(library(MASS))
suppressMessages(library(lmtest))
suppressMessages(library(sandwich))
suppressMessages(library(ggfortify))
suppressMessages(library(progress))
suppressMessages(library(bacon))
suppressMessages(library(cowplot))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(png))
suppressMessages(library(parallel))
suppressMessages(library(tidyr))
suppressMessages(library(reshape2))

# meth_matrix <- beta
# pheno_data <- p_data
# outcome <- "stotfatwh"
# predictor <- "methylation"
# covariates <- c("jsex", "sdxaage", "smokpreg", "sdxaht", "Bcell", "CD4T", "CD8T", "Mono", "NK", "nRBC")
# robust <- TRUE
# adj_method <- "BH"
# plots <- TRUE
# seed <- 48208742
# coef <- NULL
# plot_file_name <- "EWAS_lmtest_plots"
# nCores <- 40
# x <- c("cg00433159", "cg02712922", "cg18566189", "cg25000382", "cg27539927", "cg01702246", "cg24307114", "cg09164239", "cg07644413", "cg01702246")
# meth_matrix <- rbind(meth_matrix, beta[rownames(beta) %in% x,])

# source("../EWAS_lmtest_par.R")
# system.time(results <- EWAS_lmtest_par(beta, p_data, "stotfatwh", predictor="methylation", covariates, robust=TRUE, adj_method="BH", plots=TRUE, seed=1000, coef=NULL, plot_file_name="EWAS_lmtest_plots", nCores=64))

EWAS_lmtest <- function(meth_matrix, pheno_data, outcome, predictor="methylation", covariates, robust=TRUE, adj_method="BH", plots=TRUE, seed=1000, coef=NULL, plot_file_name="EWAS_lmtest_plots", nCores=1){
	
	cores_avail <- detectCores()
	if (nCores>cores_avail){
		message(paste(nCores, " cores requested but only ", cores_avail, " are available, using ", cores_avail, " cores instead for processing", sep=""))
		cores <- cores_avail
	} else if (nCores<=cores_avail){
		message(paste("Using ", nCores, " available cores for processing", sep=""))
		cores <- nCores
	}

	out <- ifelse(outcome=="methylation", predictor, outcome)

	pheno_data <- pheno_data[!is.na(pheno_data[,out]),]
	meth_matrix <- meth_matrix[,colnames(meth_matrix) %in% pheno_data[,1]]
	if (sum(pheno_data[,1] == colnames(meth_matrix), na.rm = T)== dim(meth_matrix)[2]) {
		message("Matched phenotype data")
	}else{
		message("Phenotype data and methylation data are not matched, re-ordering phenotype data to match methylation data")
		pheno_data <- pheno_data[match(colnames(meth_matrix), pheno_data[,1]),]
		if (sum(pheno_data[,1] == colnames(meth_matrix), na.rm = T)== dim(meth_matrix)[2]) {
			message("Re-ordered phenotype data is now matched")
		}
	}
	tmp_res <- list()
	plots_list <- list()

	if (outcome!="methylation"){
		if (robust & (is.numeric(pheno_data[,outcome][!is.na(pheno_data[,outcome])]))){
			type <- 1
		} else if (robust & is.factor(pheno_data[,outcome][!is.na(pheno_data[,outcome])])){
			type <- 2
		} else if (!robust & is.numeric(pheno_data[,outcome][!is.na(pheno_data[,outcome])])){
			type <- 3
		} else if (!robust & is.factor(pheno_data[,outcome][!is.na(pheno_data[,outcome])])){
			type <- 4
		}
	} else if (outcome=="methylation"){
		type <- ifelse(robust, 1, 3)
	}

	if (type==1){
		message("Outcome is continuous and robust=TRUE, carrying out robust linear regression")
		tmp_res <- mclapply(rownames(meth_matrix), function(x) {
			if (predictor=="methylation"){
				eq <- as.formula(paste(outcome, "~", x, "+", paste(covariates, collapse="+"), sep=""))
			}else if (outcome=="methylation"){
				eq <- as.formula(paste(x, "~", predictor, "+", paste(covariates, collapse="+"), sep=""))
			} else if (predictor != "methylation" | outcome != "methylation"){
				stop("One of 'predictor' or 'outcome' needs to be 'methylation', please change")
			}
			tmp <- as.data.frame(cbind(pheno_data, as.numeric(meth_matrix[x,])))
			colnames(tmp)[dim(tmp)[2]] <- x
			tmp <- tmp[order(tmp[,x], decreasing=TRUE),]
			tmp$diff <- NA
			for (j in 1:(nrow(tmp)-1)){
				tmp$diff[j] <- abs(tmp[j,x] - tmp[j+1,x])
			}
			outlier_present <- ifelse(sum(tmp$diff>0.1, na.rm=TRUE)>=1, TRUE, FALSE)
			if (outlier_present){
				if (length(which(tmp$diff>0.1))>1){
					tmp_2 <- tmp
				} else {
					if (abs(which(tmp$diff>0.1)-nrow(tmp))<abs(which(tmp$diff>0.1)-1)){
						if(which(tmp$diff>0.1)-(nrow(tmp)-1)>0){
							tmp_2 <- tmp[-seq(which(tmp$diff>0.1)+1, nrow(tmp)),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						} else {
							tmp_2 <- tmp[-(which(tmp$diff>0.1)+1),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						}
					} else if (abs(which(tmp$diff>0.1)-nrow(tmp))>abs(which(tmp$diff>0.1)-1)){
						if((which(tmp$diff>0.1)-1)>0){
							tmp_2 <- tmp[-seq(1, which(tmp$diff>0.1)),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						} else {
							tmp_2 <- tmp[-(which(tmp$diff>0.1)),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						}
					}
				}
			}
			set.seed(seed)
			mod <- rlm(eq, data=tmp)
			if (predictor != "methylation"){
				if (is.factor(pheno_data[,predictor])){
					cf <- coeftest(mod, vcov=vcovHC(mod, type="HC0"))
					cf <- as.data.frame(t(cf[coef,]))
					cis <- coefci(mod, vcov=vcovHC(mod, type="HC0"))
					cf$`Lower_95%CI` <- cis[coef,1]
					cf$`Upper_95%CI` <- cis[coef,2]
				} else {
					cf <- coeftest(mod, vcov=vcovHC(mod, type="HC0"))
					cf <- as.data.frame(t(cf[2,]))
					cis <- coefci(mod, vcov=vcovHC(mod, type="HC0"))
					cf$`Lower_95%CI` <- cis[2,1]
					cf$`Upper_95%CI` <- cis[2,2]
				}
			} else {
				cf <- coeftest(mod, vcov=vcovHC(mod, type="HC0"))
				cf <- as.data.frame(t(cf[2,]))
				cis <- coefci(mod, vcov=vcovHC(mod, type="HC0"))
				cf$`Lower_95%CI` <- cis[2,1]
				cf$`Upper_95%CI` <- cis[2,2]
			}
			cf$FLAG <- ifelse((outlier_present & length(which(tmp$diff>0.1))>1), "FLAG-check", ifelse((outlier_present & length(which(tmp$diff>0.1))==1), "FLAG", "Not flagged"))
			if (outlier_present & length(which(tmp$diff>0.1)==1)){
				set.seed(seed)
				mod <- rlm(eq, data=tmp_2)
				if (predictor != "methylation"){
					if (is.factor(pheno_data[,predictor])){
						cf2 <- coeftest(mod, vcov=vcovHC(mod, type="HC0"))
						cf2 <- as.data.frame(t(cf2[coef,]))
						cis <- coefci(mod, vcov=vcovHC(mod, type="HC0"))
						cf2$`Lower_95%CI` <- cis[coef,1]
						cf2$`Upper_95%CI` <- cis[coef,2]
					} else {
						cf2 <- coeftest(mod, vcov=vcovHC(mod, type="HC0"))
						cf2 <- as.data.frame(t(cf2[2,]))
						cis <- coefci(mod, vcov=vcovHC(mod, type="HC0"))
						cf2$`Lower_95%CI` <- cis[2,1]
						cf2$`Upper_95%CI` <- cis[2,2]
					}
				} else {
					cf2 <- coeftest(mod, vcov=vcovHC(mod, type="HC0"))
					cf2 <- as.data.frame(t(cf2[2,]))
					cis <- coefci(mod, vcov=vcovHC(mod, type="HC0"))
					cf2$`Lower_95%CI` <- cis[2,1]
					cf2$`Upper_95%CI` <- cis[2,2]
				}
				cf$new <- tidyr::unite(cf2, new, sep=",")
			} else {
				cf$new <- NA
			}
			cf
		}, mc.cores=cores)
	} else if (type==2){
		message("Outcome is a factor and robust=TRUE, carrying out robust binary logistic regression")
		tmp_res <- mclapply(rownames(meth_matrix), function(x) {
			if (progress_bar){pb$tick(1)}
			if (predictor=="methylation"){
				eq <- as.formula(paste(outcome, "~", x, "+", paste(covariates, collapse="+"), sep=""))
			}else if (outcome=="methylation"){
				stop("Methylation cannot be a factor, please change to a numeric")
			} else if (predictor != "methylation" | outcome != "methylation"){
				stop("One of 'predictor' or 'outcome' needs to be 'methylation', please change")
			}
			tmp <- as.data.frame(cbind(pheno_data, as.numeric(meth_matrix[x,])))
			colnames(tmp)[dim(tmp)[2]] <- x
			tmp <- tmp[order(tmp[,x], decreasing=TRUE),]
			tmp$diff <- NA
			for (j in 1:(nrow(tmp)-1)){
				tmp$diff[j] <- abs(tmp[j,x] - tmp[j+1,x])
			}
			outlier_present <- ifelse(sum(tmp$diff>0.1, na.rm=TRUE)>=1, TRUE, FALSE)
			if (outlier_present){
				if (length(which(tmp$diff>0.1))>1){
					tmp_2 <- tmp
				} else {
					if (abs(which(tmp$diff>0.1)-nrow(tmp))<abs(which(tmp$diff>0.1)-1)){
						if(which(tmp$diff>0.1)-(nrow(tmp)-1)>0){
							tmp_2 <- tmp[-seq(which(tmp$diff>0.1)+1, nrow(tmp)),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						} else {
							tmp_2 <- tmp[-(which(tmp$diff>0.1)+1),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						}
					} else if (abs(which(tmp$diff>0.1)-nrow(tmp))>abs(which(tmp$diff>0.1)-1)){
						if((which(tmp$diff>0.1)-1)>0){
							tmp_2 <- tmp[-seq(1, which(tmp$diff>0.1)),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						} else {
							tmp_2 <- tmp[-(which(tmp$diff>0.1)),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						}
					}
				}
			}
			set.seed(seed)
			mod <- robustbase::glmrob(eq, family="binomial", data=tmp)
			cf <- as.data.frame(t(coef(summary(mod))[2,]))
			cis <- confint.default(mod)[2,]
			cf$`Lower_95%CI` <- cis[1]
			cf$`Upper_95%CI` <- cis[2]
			cf$FLAG <- ifelse((outlier_present & length(which(tmp$diff>0.1))>1), "FLAG-check", ifelse((outlier_present & length(which(tmp$diff>0.1))==1), "FLAG", "Not flagged"))
			if (outlier_present & length(which(tmp$diff>0.1)==1)){
				set.seed(seed)
				mod <- robustbase::glmrob(eq, family="binomial", data=tmp_2)
				cf2 <- as.data.frame(t(coef(summary(mod))[2,]))
				cis <- confint.default(mod)[2,]
				cf2$`Lower_95%CI` <- cis[1]
				cf2$`Upper_95%CI` <- cis[2]
				cf$new <- tidyr::unite(cf2, new, sep=",")
			} else {
				cf$new <- NA
			}
			cf
		}, mc.cores=cores )
	} else if (type==3){
		message("Outcome is continuous and robust=FALSE, carrying out ordinary least squares regression")
		tmp_res <- lapply(rownames(meth_matrix), function(x) {
			if (progress_bar){pb$tick(1)}
			if (predictor=="methylation"){
				eq <- as.formula(paste(outcome, "~", x, "+", paste(covariates, collapse="+"), sep=""))
			}else if (outcome=="methylation"){
				eq <- as.formula(paste(x, "~", predictor, "+", paste(covariates, collapse="+"), sep=""))
			} else if (predictor != "methylation" | outcome != "methylation"){
				stop("One of 'predictor' or 'outcome' needs to be 'methylation', please change")
			}
			tmp <- as.data.frame(cbind(pheno_data, as.numeric(meth_matrix[x,])))
			colnames(tmp)[dim(tmp)[2]] <- x
			tmp <- tmp[order(tmp[,x], decreasing=TRUE),]
			tmp$diff <- NA
			for (j in 1:(nrow(tmp)-1)){
				tmp$diff[j] <- abs(tmp[j,x] - tmp[j+1,x])
			}
			outlier_present <- ifelse(sum(tmp$diff>0.1, na.rm=TRUE)>=1, TRUE, FALSE)
			if (outlier_present){
				if (length(which(tmp$diff>0.1))>1){
					tmp_2 <- tmp
				} else {
					if (abs(which(tmp$diff>0.1)-nrow(tmp))<abs(which(tmp$diff>0.1)-1)){
						if(which(tmp$diff>0.1)-(nrow(tmp)-1)>0){
							tmp_2 <- tmp[-seq(which(tmp$diff>0.1)+1, nrow(tmp)),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						} else {
							tmp_2 <- tmp[-(which(tmp$diff>0.1)+1),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						}
					} else if (abs(which(tmp$diff>0.1)-nrow(tmp))>abs(which(tmp$diff>0.1)-1)){
						if((which(tmp$diff>0.1)-1)>0){
							tmp_2 <- tmp[-seq(1, which(tmp$diff>0.1)),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						} else {
							tmp_2 <- tmp[-(which(tmp$diff>0.1)),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						}
					}
				}
			}
			set.seed(seed)
			mod <- lm(eq, data=tmp)
			if (predictor != "methylation"){
				if (is.factor(pheno_data[,predictor])){
					cf <- as.data.frame(t(coef(summary(mod))[coef,]))
					cis <- confint(mod)[coef,]
					cf$`Lower_95%CI` <- cis[1]
					cf$`Upper_95%CI` <- cis[2]
				} else {
					cf <- as.data.frame(t(coef(summary(mod))[2,]))
					cis <- confint(mod)[2,]
					cf$`Lower_95%CI` <- cis[1]
					cf$`Upper_95%CI` <- cis[2]
				}
			} else {
				cf <- as.data.frame(t(coef(summary(mod))[2,]))
				cis <- confint(mod)[2,]
				cf$`Lower_95%CI` <- cis[1]
				cf$`Upper_95%CI` <- cis[2]
			}
			cf$FLAG <- ifelse((outlier_present & length(which(tmp$diff>0.1))>1), "FLAG-check", ifelse((outlier_present & length(which(tmp$diff>0.1))==1), "FLAG", "Not flagged"))
			if (outlier_present & length(which(tmp$diff>0.1)==1)){
				mod <- lm(eq, data=tmp)
				if (predictor != "methylation"){
					if (is.factor(pheno_data[,predictor])){
						cf2 <- as.data.frame(t(coef(summary(mod))[coef,]))
						cis <- confint(mod)[coef,]
						cf2$`Lower_95%CI` <- cis[1]
						cf2$`Upper_95%CI` <- cis[2]
					} else {
						cf <- as.data.frame(t(coef(summary(mod))[2,]))
						cis <- confint(mod)[2,]
						cf2$`Lower_95%CI` <- cis[1]
						cf2$`Upper_95%CI` <- cis[2]
					}
				} else {
					cf2 <- as.data.frame(t(coef(summary(mod))[2,]))
					cis <- confint(mod)[2,]
					cf2$`Lower_95%CI` <- cis[1]
					cf2$`Upper_95%CI` <- cis[2]
				}
				cf$new <- tidyr::unite(cf2, new, sep=",")
			} else {
				cf$new <- NA
			}
			cf
		} )
	} else if (type==4){
		message("Outcome is a factor and robust=FALSE, carrying out binary logistic regression")
		tmp_res <- lapply(rownames(meth_matrix), function(x) {
			if (progress_bar){pb$tick(1)}
			if (predictor=="methylation"){
				eq <- as.formula(paste(outcome, "~", x, "+", paste(covariates, collapse="+"), sep=""))
			}else if (outcome=="methylation"){
				stop("Methylation cannot be a factor, please change to a numeric")
			} else if (predictor != "methylation" | outcome != "methylation"){
				stop("One of 'predictor' or 'outcome' needs to be 'methylation', please change")
			}
			tmp <- as.data.frame(cbind(pheno_data, as.numeric(meth_matrix[x,])))
			colnames(tmp)[dim(tmp)[2]] <- x
			tmp <- tmp[order(tmp[,x], decreasing=TRUE),]
			tmp$diff <- NA
			for (j in 1:(nrow(tmp)-1)){
				tmp$diff[j] <- abs(tmp[j,x] - tmp[j+1,x])
			}
			outlier_present <- ifelse(sum(tmp$diff>0.1, na.rm=TRUE)>=1, TRUE, FALSE)
			if (outlier_present){
				if (length(which(tmp$diff>0.1))>1){
					tmp_2 <- tmp
				} else {
					if (abs(which(tmp$diff>0.1)-nrow(tmp))<abs(which(tmp$diff>0.1)-1)){
						if(which(tmp$diff>0.1)-(nrow(tmp)-1)>0){
							tmp_2 <- tmp[-seq(which(tmp$diff>0.1)+1, nrow(tmp)),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						} else {
							tmp_2 <- tmp[-(which(tmp$diff>0.1)+1),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						}
					} else if (abs(which(tmp$diff>0.1)-nrow(tmp))>abs(which(tmp$diff>0.1)-1)){
						if((which(tmp$diff>0.1)-1)>0){
							tmp_2 <- tmp[-seq(1, which(tmp$diff>0.1)),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						} else {
							tmp_2 <- tmp[-(which(tmp$diff>0.1)),]
							tmp_2 <- tmp_2[,-ncol(tmp_2)]
						}
					}
				}
			}
			set.seed(seed)
			mod <- glm(eq, family="binomial", data=tmp)
			cf <- as.data.frame(t(coef(summary(mod))[2,]))
			cis <- suppressMessages(confint(mod)[2,])
			cf$`Lower_95%CI` <- cis[1]
			cf$`Upper_95%CI` <- cis[2]
			cf$FLAG <- ifelse((outlier_present & length(which(tmp$diff>0.1))>1), "FLAG-check", ifelse((outlier_present & length(which(tmp$diff>0.1))==1), "FLAG", "Not flagged"))
			if (outlier_present & length(which(tmp$diff>0.1)==1)){
				set.seed(seed)
				mod <- robustbase::glmrob(eq, family="binomial", data=tmp_2)
				cf2 <- as.data.frame(t(coef(summary(mod))[2,]))
				cis <- confint.default(mod)[2,]
				cf2$`Lower_95%CI` <- cis[1]
				cf2$`Upper_95%CI` <- cis[2]
				cf$new <- tidyr::unite(cf2, new, sep=",")
			} else {
				cf$new <- NA
			}
			cf
		} )
	}
	names(tmp_res) <- rownames(meth_matrix)
	probenames <- names(tmp_res)
	col_names <- c("Estimate", "Std. Error", "statistic", "P.Value", "Lower_95%CI", "Upper_95%CI", "FLAG", "new")
	all_res <- data.frame(matrix(unlist(tmp_res), nrow=length(tmp_res), byrow=T))
	rownames(all_res) <- probenames
	colnames(all_res) <- col_names
	all_res$adj.P.Val <- p.adjust(all_res$P.Value, method=adj_method)
	all_res$FLAG <- factor(all_res$FLAG, levels=c("Not flagged", "FLAG", "FLAG-check"))
	all_res <- all_res[order(all_res$FLAG, all_res$adj.P.Val, all_res$P.Value, decreasing=FALSE),]
	all_res$new <- ifelse(!is.na(all_res$new), paste(rownames(all_res), all_res$new, sep=","), NA)
	flagged_tmp <- all_res[all_res$FLAG %in% c("FLAG", "FLAG-check"),c("new", "FLAG")]
	flag <- flagged_tmp$FLAG
	flagged_results <- as.character(flagged_tmp$new)
	flagged_results <- reshape2::colsplit(flagged_results, ",", c("Probe", col_names[1:6]))
	flagged_results$FLAG <- flag
	all_res <- all_res[,-8]
	flagged_results <- flagged_results[,-9]
	all_res[,1:6] <- lapply(all_res[,1:6], function(x) as.numeric(as.character(x)))

	if (plots){
		plots_list <- list()
		num <- 0
		cgs <- rownames(all_res)[1:5]
		if ((type==1&predictor=="methylation")|(type==3&predictor=="methylation")){
			for (id in cgs){
				num <- num+1
				eq <- as.formula(paste(outcome, "~", id, "+", paste(covariates, collapse="+"), sep=""))
				tmp <- as.data.frame(cbind(pheno_data, as.numeric(meth_matrix[id,])))
				colnames(tmp)[dim(tmp)[2]] <- id
				tmp <- tmp[!is.na(tmp[,outcome]),]
				if (type==1){
					set.seed(seed)
					mod <- try(rlm(eq, data=tmp))
					plots_list[[num]] <- ggplot2::autoplot(mod, which=1:3, ncol=2)
					tmp_p <- ggplot(tmp, aes_string(x=id, y=outcome)) + geom_point() + geom_smooth(method="rlm") + ggtitle("Scatter plot")
					plots_list[[num]] <- cowplot::plot_grid(plots_list[[num]]@plots[[1]], plots_list[[num]]@plots[[2]], plots_list[[num]]@plots[[3]], tmp_p, ncol=2)
				} else if (type==3){
					set.seed(seed)
					mod <- try(lm(eq, data=tmp))
					plots_list[[num]] <- ggplot2::autoplot(mod, which=1:6, ncol=2)
					tmp_p <- ggplot(tmp, aes_string(x=id, y=outcome)) + geom_point() + geom_smooth(method="lm") + ggtitle("Scatter plot")
					plots_list[[num]] <- cowplot::plot_grid(plots_list[[num]]@plots[[1]], plots_list[[num]]@plots[[2]], plots_list[[num]]@plots[[3]], plots_list[[num]]@plots[[4]], plots_list[[num]]@plots[[5]], plots_list[[num]]@plots[[6]],tmp_p, ncol=2)
				}
			}
		} else if ((type==1&is.numeric(pheno_data[,predictor]))|(type==3&is.numeric(pheno_data[,predictor]))){
 			for (id in cgs){
				num <- num+1
				eq <- as.formula(paste(id, "~", predictor, "+", paste(covariates, collapse="+"), sep=""))
				tmp <- as.data.frame(cbind(pheno_data, as.numeric(meth_matrix[id,])))
				colnames(tmp)[dim(tmp)[2]] <- id
				tmp <- tmp[!is.na(tmp[,predictor]),]
				if (type==1){
					set.seed(seed)
					mod <- try(rlm(eq, data=tmp))
					plots_list[[num]] <- ggplot2::autoplot(mod, which=1:3, ncol=2)
					tmp_p <- ggplot(tmp, aes_string(x=predictor, y=id)) + geom_point() + geom_smooth(method="rlm") + ggtitle("Scatter plot")
					plots_list[[num]] <- cowplot::plot_grid(plots_list[[num]]@plots[[1]], plots_list[[num]]@plots[[2]], plots_list[[num]]@plots[[3]], tmp_p, ncol=2)
				} else if (type==3){
					set.seed(seed)
					mod <- try(lm(eq, data=tmp))
					plots_list[[num]] <- ggplot2::autoplot(mod, which=1:6, ncol=2)
					tmp_p <- ggplot(tmp, aes_string(x=predictor, y=id)) + geom_point() + geom_smooth(method="lm") + ggtitle("Scatter plot")
					plots_list[[num]] <- cowplot::plot_grid(plots_list[[num]]@plots[[1]], plots_list[[num]]@plots[[2]], plots_list[[num]]@plots[[3]], plots_list[[num]]@plots[[4]], plots_list[[num]]@plots[[5]], plots_list[[num]]@plots[[6]],tmp_p, ncol=2)
				}
			}
		} else if (type==2|type==4){
			message("Only generating plots for top 5 dmCpGs, binary logistic regression assumptions need to be checked manually")
			for (id in cgs){
				num <- num+1
				tmp <- as.data.frame(cbind(pheno_data, as.numeric(meth_matrix[id,])))
				colnames(tmp)[dim(tmp)[2]] <- id
				tmp <- tmp[!is.na(tmp[,outcome]),]
				tmp_p <- ggplot(tmp, aes_string(x=outcome, y=id)) + geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.1)
				plots_list[[num]] <- tmp_p
			}
		} else if ((type==1&is.factor(pheno_data[,predictor]))|(type==3&is.factor(pheno_data[,predictor]))){
			for (id in cgs){
				num <- num+1
				eq <- as.formula(paste(id, "~", predictor, "+", paste(covariates, collapse="+"), sep=""))
				tmp <- as.data.frame(cbind(pheno_data, as.numeric(meth_matrix[id,])))
				colnames(tmp)[dim(tmp)[2]] <- id
				tmp <- tmp[!is.na(tmp[,predictor]),]
				tmp_p <- ggplot(tmp, aes_string(x=predictor, y=id)) + geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.1)
				set.seed(seed)
				mod <- try(rlm(eq, data=tmp))
				if (type==1){
					plots_list[[num]] <- ggplot2::autoplot(mod, which=1:3, ncol=2)
					plots_list[[num]] <- cowplot::plot_grid(plots_list[[num]]@plots[[1]], plots_list[[num]]@plots[[2]], plots_list[[num]]@plots[[3]], tmp_p, ncol=2)
				} else if (type==3){
					plots_list[[num]] <- ggplot2::autoplot(mod, which=1:6, ncol=2)
					plots_list[[num]] <- cowplot::plot_grid(plots_list[[num]]@plots[[1]], plots_list[[num]]@plots[[2]], plots_list[[num]]@plots[[3]], plots_list[[num]]@plots[[4]], plots_list[[num]]@plots[[5]], plots_list[[num]]@plots[[6]],tmp_p, ncol=2)
				}
			}
		}
		phist <- ggplot(all_res, aes(x=P.Value)) + geom_histogram(bins=100, boundary=0)
		qqp <- ggplot(all_res) +
				geom_point(aes(y=-log10(sort(P.Value)), x=-log10(ppoints(length(all_res$P.Value)))), shape = 1, size = 3) +
				geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
				geom_line(aes(x=-log10(ppoints(length(all_res$P.Value))), y=-log10(qbeta(p=(1-0.95)/2, shape1=1:length(all_res$P.Value), shape2=length(all_res$P.Value):1))), linetype = 2) +
				geom_line(aes(x=-log10(ppoints(length(all_res$P.Value))), y=-log10(qbeta(p=(1+0.95)/2, shape1=1:length(all_res$P.Value), shape2=length(all_res$P.Value):1))), linetype = 2) +
				xlab(expression(paste("Expected -log"[10], plain(P)))) + ylab(expression(paste("Observed -log"[10], plain(P))))

		options(bitmapType="cairo")
		png("QQ_plot.png", res=1200, width=5.25, height=5.25, units="in", pointsize=4)
		par(mar= c(5, 5, 2, 2), xaxs="i", yaxs="i", cex.axis=2, cex.lab=2)
		print(qqp)
		dev.off()
		qqp2 <- grid::rasterGrob(png::readPNG("QQ_plot.png", native = FALSE), interpolate = FALSE)

		height <- ifelse(type==3, 11, 7)
		pdf(paste(plot_file_name, ".pdf", sep=""), height=height)
		gridExtra::grid.arrange(qqp2)
		print(phist)
		for (k in 1:5) {
			print(plots_list[[k]])
		}
		dev.off()
		file.remove("QQ_plot.png")
	}

	lambda <- median(qchisq(all_res$P.Value, df=1, lower.tail=F))/qchisq(0.5, df=1)
	if (lambda>1.2){
		bc <- bacon(teststatistics=all_res$statistic)
		ps <- as.data.frame(pval(bc))
		all_res$Bacon_P <- ps$V1
		all_res$Bacon_FDR <- p.adjust(ps$V1, method=adj_method)
		lambda <- list(pre_Bacon=lambda, post_Bacon=median(qchisq(all_res$Bacon_P, df=1, lower.tail=F))/qchisq(0.5, df=1))
	}
	if (plots){
		return(list(results=all_res, flagged_results=flagged_results, lambda=lambda, plots=list(qqp, phist, plots_list)))
	} else {
		return(list(results=all_res, flagged_results=flagged_results, lambda=lambda))
	}
}
