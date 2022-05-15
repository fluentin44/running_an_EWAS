# Notes:
## FIRST, before you get into r - Navigate to the project folder - cd /scratch/<YOUR USERNAME>/starting_an_ewas/
## Make sure your combat data and phenotypic data are both in the "data" folder
## Anything between arrows needs to be changed, i.e. <YOUR USERNAME> and remember to move the arrows too!

# Setting working directory *EDIT THIS SECTION* -------------------------------------------------------------------------------

setwd("/mainfs/scratch/<YOUR USERNAME>/running_an_ewas/") # need to change this to your wd

# Setting up some options ----------------------------------------------------------------------------------------

capabilities()
options(tibble.width = Inf)
options(bitmapType="cairo")

# Load packages --------------------------------------------------------------------------------------------------

suppressPackageStartupMessages(library(ChAMP))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(bacon))
suppressPackageStartupMessages(library(DMRcate))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICmanifest))

# Source scripts --------------------------------------------------------------------------------------------------

source("src/limma_scripts/gen_report_type.R")
source("src/limma_scripts/gen_lambdas.R")
source("src/limma_scripts/gen_report.R")
source("src/limma_scripts/bacon_adj.R")
source("src/limma_scripts/run_ewas.R")
source("src/dmrcate.R")
source("src/miss_methyl.R")

# Source scripts --------------------------------------------------------------------------------------------------

## Import samplesheets, annot, combatted meth values and nullify pool_ID
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annot <- annot[,c(22,24,19,1:4,31,23)]
annot[1:5,1:5]

## Set up vars
date_run_started <- Sys.Date()

dir.create("./results/lambdas_and_hits/", recursive = T)
dir.create(paste0("./results/", date_run_started, "/"), recursive = T)
dir.create(paste0("./results/", date_run_started, "/", "dmrcate", "/"))
dir.create(paste0("./results/", date_run_started, "/", "missmethyl", "/"))
model_save_path <- paste0("./results/", date_run_started, "/")
dmrcate_path <- paste0("./results/", date_run_started, "/", "dmrcate", "/")
missmethyl_path <- paste0("./results/", date_run_started, "/", "missmethyl", "/")
lambdas_path <- paste0("./results/lambdas_and_hits/", date_run_started, "_", "lambdas_and_hits.txt")

# Import betas and pheotypic data *EDIT THIS SECTION* ---------------------------------------------------------------------------------

# Three options depending what filetype your combat beta is:
combat_beta <- read.table("data/<PATH TO YOUR COMBAT DATA FILE WITH FILE EXTENSION>", sep=",", header=T)         # a csv combat beta
combat_beta <- read.csv("data/<PATH TO YOUR COMBAT DATA FILE WITH FILE EXTENSION>", sep=",", header=T)           # a csv combat beta
load("<PATH TO COMBAT BETA>")                                                           # An Rdata combat beta - this will just load an object - wll need to find out what its called with ls()

class(combat_beta)                                                                      # Needs to say "matrix" if not, contact me

p_data <- read.table("data/<PATH TO YOUR PHENOTYPIC DATA FILE WITH FILE EXTENSION>", sep=",", header=T)          # What is your phenotypic data called

str(p_data)                                              # Can let you look at your p_data, what format columns are in
glimpse(p_data)                                          # Another option to do the same thing

p_data <-
  p_data %>%
  select(<A>,
         <B>,                                            # X, Y and Z are any variables of interest or covariates
         <C>                                             # Sometimes is useful to make a smaller table if you have a lot of pheotypic data
         ) %>%
  mutate(across(c(jsex, position_on_array), as_factor),  # Your phneotypic data columns NEEDS to be in the right format (numeric/integer, character, factor) two ways to do that are 72 & 73
         across(c(1), as.character)                      # FOR EXAMPLE: Line 72 will convert jsex and position_on_array to factors, and line 73 will convert the first column to character
         )




# Differential methylation analysis *EDIT THIS SECTION* -------------------------------------------------------------------------------

var <- c(                                # Need to enter your variable of interest, i.e. the subject of the regression - can put in multiple or only one.
 "<A>",
 "<B>",
 "<C>"
)
var_extra <- "base_model"                # PLEASE CHANGE - This can be a little bit of info to distinguish one report from another, i.e. 'no RBCs'
covariates <- c("jsex",                  # PLEASE CHANGE - This is all the covariates in your regression (left common examples in but replace as needed!!)
                "position_on_array")

# Regression script -------------------------------------------------------------------------------------------------------

for(var in var){                         # This will successuvely work through every variable in var and run cde between the curly brckets { }

limma_save_path <- paste0(model_save_path, date_run_started, "_", var, "_", var_extra, "_results.txt")
output_file_name <- paste0(date_run_started, "-", var, "_", var_extra)

ewas_objects <-
run.ewas(var,
         var_extra,
         combat_beta,
         p=p_data,
         covariates=covariates,
         results_path=limma_save_path,
         test=T
         #sva=F,
         #first_run=F,
         #no_svs=NULL
         )

# DMRcate - This code will look for DMRs and save the output as an object for an FDR <0.05 (which can be changed)
dmroutput <- EPICdmrcate(combat_beta=ewas_objects$combat, design=ewas_objects$mod, coef=2, fdr=0.05)
#results.ranges <- extractRanges(dmroutput, genome="hg19")
saveRDS(dmroutput, paste0(dmrcate_path, var, "_", var_extra, ".rds"))

#missMethyl - This will run missMethyl if you have hits FDR<0.25
mm_it(table=ewas_objects$table,
      missmethyl_path=missmethyl_path,
      var=var,
      var_extra=var_extra)
}
