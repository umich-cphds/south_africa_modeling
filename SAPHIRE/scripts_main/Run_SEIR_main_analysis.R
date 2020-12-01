## IMPORTANT: Please set code_root variable properly. 
## code_root should be set to the directory
code_root=""
setwd('')
library(BayesianTools)
library(vioplot)
library(corrplot)
library(readr)
library(cairoDevice)
library(tidyverse)
##
source(paste0(code_root, "R/fun_SEIRpred.R"))
source(paste0(code_root, "R/fun_SEIRsimu.R"))
source(paste0(code_root, "R/fun_SEIRfitting.R"))
source(paste0(code_root, "R/init_cond.R"))
source(paste0(code_root, "R/fun_R0estimate.R"))

init_sets_list=get_init_sets_list(r0 = 0.10)###initial ascertained rate
set.seed('20200601')
SEIRfitting(init_sets_list, randomize_startValue = T,
            run_id = "analysis_0.10", output_ret = T, skip_MCMC=F)

