rm(list=ls())
# stress and fetal growth analysis 

load("/Users/danagoin/Documents/Research projects/CiOB-ECHO/Projects/Stress and birth weight/data/stress_fetal_growth_imp")

# g-computation for individual exposures and sum 
# first check positivity 
library(mice)
library(dplyr)
library(corrplot)
library(RColorBrewer)
library(ggplot2)
library(msm)


stressor_measures <- exposures <- c("econ_cat", "food_cat", "job_highstrain_missing", "poor_perceived_nhbd", 
                     "low_perceived_socialstat", 
                     "caregiving_missing", "high_stressful_events", "unplanned_pregnancy")


# --------------------------------------------------------------------------------------------------------------------------
# estimate individual exposure propensity scores and record which people have propensity scores within the area of support
# --------------------------------------------------------------------------------------------------------------------------

pscore <- function(i, j) {
df_t <- mice::complete(imp, i)
  x <- predict(with(df_t, glm(get(stressor_measures[j]) ~  mat_age + mat_race_eth + mat_edu + parity + marital, family="binomial")), type="response")
 return(x)
}

include <- list()

for (j in 1:length(stressor_measures)) {
  
# calculate propensity scores for each imputed data set 
pscores_list <- lapply(1:30, function(x) pscore(x, j))
pscores <- do.call(cbind, pscores_list) 

# make matrix of exposure levels for each imputed data set 
observed_ls <- lapply(1:30, function(x) mice::complete(imp, x) %>% select(stressor_measures[j]))
observed_e <- do.call(cbind, observed_ls)

# find the minimum of propensity scores for people who were exposed
# and the maximum of propensity scores for people who were unexposed 
# to ensure there is adequate support 

pscores_support_ls <- lapply(1:30, function(x) cbind(min(pscores[,x][observed_e[,x]==1]), max(pscores[,x][observed_e[,x]==0])))
pscores_support <- do.call(rbind, pscores_support_ls)

# identify whether propensity scores are within the area of support for each imputed data set 
include_ls <- lapply(1:30, function(x) as.numeric(pscores[,x] >= pscores_support[x,1] & pscores[,x]<= pscores_support[x,2]))

include[[j]] <- do.call(cbind, include_ls)
}

names(include) <- stressor_measures


# --------------------------------------------------------------------------------------------------------------------------
# estimate association term birth weight for individual exposures  
# --------------------------------------------------------------------------------------------------------------------------

gcomp_mi <- function(m, exposure) {

df_p <- mice::complete(imp, m)

# include those who fell within the area of support for each imputed data set 
df_p <- df_p %>% filter(include[[exposure]][,m]==1)

df_p$exposure <- df_p[[exposure]]

fit <- glm(baby_wt_grams ~ exposure + mat_age + mat_edu + mat_race_eth + parity + marital + pp_bmi, data=df_p)
nobs <- nobs(fit)

# calculate the ATE via just the beta coefficient -- this works because it's a linear model 
ate <- coef(fit)["exposure"] 

cov <- vcov(fit)
cov <- cov["exposure","exposure"]

results <- data.frame(cbind(ate = ate, var = cov, N=nobs)) 
return(results)
}

# --------------------------------------------------------------------------------------------------------------------------
# combine results across imputations using Rubin's combining rules and write results to file    
# --------------------------------------------------------------------------------------------------------------------------

all_results <- data.frame(matrix(ncol=5, nrow=length(exposures)))

for (i in 1:length(exposures)) {

# variance calculations from https://amstat.tandfonline.com/doi/pdf/10.1080/01621459.1996.10476908?needAccess=true
# var = average of complete data variances + variance of complete data means 
ate_M_ls <- lapply(1:30, function(x) gcomp_mi(x, exposures[i]))
ate_M_df <- data.frame(do.call(rbind, ate_M_ls))
names(ate_M_df) <- c("ate", "var", "N")

results <- data.frame(cbind(rd = mean(ate_M_df$ate), var = mean(ate_M_df$var) + ((30+1)/30)*var(ate_M_df$ate)))
results$lb <- results$rd - 1.96*sqrt(results$var)
results$ub <- results$rd + 1.96*sqrt(results$var)
results$N <- mean(ate_M_df$N)

all_results[i,] <- results 
}

rownames(all_results) <- exposures
colnames(all_results) <-  c("rd","var","lb","ub","N")

write.csv(all_results, file="/Users/danagoin/Documents/Research projects/CiOB-ECHO/Projects/Stress and birth weight/results/R1/stress_bw_indiv_R1.csv")


# --------------------------------------------------------------------------------------------------------------------------
# estimate mixture propensity scores -- propensity scores for each exposure while controlling for the others 
# --------------------------------------------------------------------------------------------------------------------------

# currently exposures are all stressors. decide if you want to focus on specific ones 

pscore_mix <- function(i, j) {
  df_t <- mice::complete(imp, i)
  exposure <- exposures[j]
  
  coexposures <- exposures[-(match(exposure, exposures))]
  
  # regress exposure of interest on other exposures plus covariates 
  formula <- paste(exposure, "~",   paste(coexposures, collapse="+ "), "+ mat_age + mat_edu + mat_race_eth + parity + marital")
  
  
  x <- predict(glm(formula=formula, data=df_t, family="binomial"), type="response")
  return(x)
}


include_mix <- list()

for (j in 1:length(exposures)) {
  
  # calculate propensity scores for each imputed data set 
  pscores_list <- lapply(1:30, function(x) pscore_mix(x, j))
  pscores_mix <- do.call(cbind, pscores_list) 
  
  # make matrix of exposure levels for each imputed data set 
  observed_ls <- lapply(1:30, function(x) mice::complete(imp, x) %>% select(stressor_measures[j]))
  observed_e <- do.call(cbind, observed_ls)
  
  # find the minimum of propensity scores for people who were exposed
  # and the maximum of propensity scores for people who were unexposed 
  # to ensure there is adequate support 
  
  pscores_support_ls <- lapply(1:30, function(x) cbind(min(pscores_mix[,x][observed_e[,x]==1]), max(pscores_mix[,x][observed_e[,x]==0])))
  pscores_support_mix <- do.call(rbind, pscores_support_ls)
  
  # identify whether propensity scores are within the area of support for each imputed data set 
  include_ls <- lapply(1:30, function(x) as.numeric(pscores_mix[,x] >= pscores_support_mix[x,1] & pscores_mix[,x]<= pscores_support_mix[x,2]))
  
  include_mix[[j]] <- do.call(cbind, include_ls)
}

names(include_mix) <- stressor_measures

# --------------------------------------------------------------------------------------------------------------------------
# estimate G-computation for birth weight for individual exposures while controlling for the others 
# --------------------------------------------------------------------------------------------------------------------------


gcomp_mi_mix <- function(m, exposure) {
  
  df_p <- mice::complete(imp, m)
  
  df_p <- df_p %>% filter(include_mix[[exposure]][,m]==1)
  
  df_p$exposure <- df_p[[exposure]]
  
  
  coexposures <- exposures[-(match(exposure, exposures))]
  
  
  formula <- paste("baby_wt_grams ~ exposure +",   paste(coexposures, collapse="+ "), "+ mat_age + mat_edu + mat_race_eth + parity + marital + pp_bmi")
   
  fit <- glm(formula, data=df_p)
  nobs <- nobs(fit)
  
  # calculate the ATE via the beta coefficient -- this works because it's a linear model 
  ate <- coef(fit)["exposure"] 
  
  cov <- vcov(fit)
  cov <- cov["exposure","exposure"]
  
  results <- data.frame(cbind(ate = ate, var = cov, N=nobs)) 
  return(results)
}


# --------------------------------------------------------------------------------------------------------------------------
# combine results across imputations using Rubin's combining rules and write results to file  
# --------------------------------------------------------------------------------------------------------------------------

all_results_mix <- data.frame(matrix(ncol=5, nrow=length(exposures)))

for (i in 1:length(exposures)) {
# variance calculations from https://amstat.tandfonline.com/doi/pdf/10.1080/01621459.1996.10476908?needAccess=true
# var = average of complete data variances + variance of complete data means 
ate_M_ls <- lapply(1:30, function(x) gcomp_mi_mix(x, exposures[i]))
ate_M_df <- data.frame(do.call(rbind, ate_M_ls))
names(ate_M_df) <- c("ate", "var", "N")

results <- data.frame(cbind(rd = mean(ate_M_df$ate), var = mean(ate_M_df$var) + ((30+1)/30)*var(ate_M_df$ate)))
results$lb <- results$rd - 1.96*sqrt(results$var)
results$ub <- results$rd + 1.96*sqrt(results$var)
results$N <- mean(ate_M_df$N)

all_results_mix[i, ] <- results
}
rownames(all_results_mix) <- exposures 
colnames(all_results_mix) <- c("rd","var","lb","ub","N")


write.csv(all_results_mix, file="/Users/danagoin/Documents/Research projects/CiOB-ECHO/Projects/Stress and birth weight/results/R1/stress_bw_mix_R1.csv")


# --------------------------------------------------------------------------------------------------------------------------
# estimate G-computation for birth weight for joint effect, all exposures at once
# --------------------------------------------------------------------------------------------------------------------------

# don't have anyone who experience all or none of the stressors who have proensity scores in the area of common support 
# are there specific combinations we would be interested in? 
# maybe econ_cat with unplanned_pregnancy, food_cat?
# caregiving with jobstrain? 
# individual associations while controlling for the other exposures  

# consider looking at pairs of exposures -- include interactions? 
# use delta method rather than bootstrap for birth weight analyses 
# -- but also look into why you are getting 0 variance when you bootstrap. issue with sample size?

gcomp_mi_joint <- function(m, ej) {
  
  df_p <- mice::complete(imp, m)
  
  df_p <- df_p %>% filter(include_mix[[exposures_joint[[ej]][1]]][,m]==1 & include_mix[[exposures_joint[[ej]][2]]][,m]==1)
  
  coexposures <- exposures[-(match(exposures_joint[[ej]], exposures))]
  
  
  formula <- paste("baby_wt_grams ~ ", paste0(exposures_joint[[ej]], collapse="*"), " +", paste(coexposures, collapse=" + "), "+ mat_age + mat_edu + mat_race_eth + parity + marital + pp_bmi")
  
  fit <- glm(formula, data=df_p)
  nobs <- nobs(fit)
  
  # calculate the ATE via sum of beta coefficients -- this works because it's a linear model 
  ate <- sum(coef(fit)[c(exposures_joint[[ej]], paste0(exposures_joint[[ej]], collapse=":"))])

  # calculate the SE via the delta method 
    
  cov <- vcov(fit)
  cov <- cov[c(exposures_joint[[ej]], paste0(exposures_joint[[ej]], collapse=":")),c(exposures_joint[[ej]], paste0(exposures_joint[[ej]], collapse=":"))]
  
  dm_var <- deltamethod(~ x1 + x2 + x3, 
              mean = coef(fit)[c(exposures_joint[[ej]], paste0(exposures_joint[[ej]], collapse=":"))], cov = cov, ses=T)^2
  
  
  results <- data.frame(cbind(ate = ate, var = dm_var, N=nobs)) 
  return(results)
}

# --------------------------------------------------------------------------------------------------------------------------
# combine results across imputations using Rubin's combining rules and write results to file  
# --------------------------------------------------------------------------------------------------------------------------

exposures_joint <- list(c("food_cat","unplanned_pregnancy"), 
                        c("high_stressful_events","unplanned_pregnancy"), 
                        c("food_cat", "high_stressful_events"), 
                        c("job_highstrain_missing", "caregiving_missing"), 
                        c("job_highstrain_missing", "high_stressful_events"), 
                        c("high_stressful_events", "caregiving_missing"), 
                        c("food_cat", "caregiving_missing"), 
                        c("job_highstrain_missing", "food_cat"), 
                        c("job_highstrain_missing", "unplanned_pregnancy"), 
                        c("unplanned_pregnancy", "caregiving_missing"), 
                        c("econ_cat","food_cat"),
                        c("econ_cat","unplanned_pregnancy"),
                        c("econ_cat","high_stressful_events"),
                        c("econ_cat","job_highstrain_missing"),
                        c("econ_cat", "caregiving_missing"), 
                        c("poor_perceived_nhbd","econ_cat"), 
                        c("poor_perceived_nhbd", "food_cat"), 
                        c("poor_perceived_nhbd", "job_highstrain_missing"), 
                        c("poor_perceived_nhbd", "low_perceived_socialstat"), 
                        c("poor_perceived_nhbd", "caregiving_missing"), 
                        c("poor_perceived_nhbd", "high_stressful_events"), 
                        c("poor_perceived_nhbd", "unplanned_pregnancy"), 
                        c("low_perceived_socialstat", "econ_cat"), 
                        c("low_perceived_socialstat", "food_cat"), 
                        c("low_perceived_socialstat", "job_highstrain_missing"), 
                        c("low_perceived_socialstat", "caregiving_missing"), 
                        c("low_perceived_socialstat", "high_stressful_events"), 
                        c("low_perceived_socialstat", "unplanned_pregnancy"))



  # variance calculations from https://amstat.tandfonline.com/doi/pdf/10.1080/01621459.1996.10476908?needAccess=true
  # var = average of complete data variances + variance of complete data means 

all_results_joint_ls <- list()

for (ej in 1:length(exposures_joint)) {
  ate_M_ls <- lapply(1:30, function(x) gcomp_mi_joint(x, ej))
  ate_M_df <- data.frame(do.call(rbind, ate_M_ls))
  names(ate_M_df) <- c("ate", "var", "N")
  
  results <- data.frame(cbind(rd = mean(ate_M_df$ate), var = mean(ate_M_df$var) + ((30+1)/30)*var(ate_M_df$ate)))
  results$lb <- results$rd - 1.96*sqrt(results$var)
  results$ub <- results$rd + 1.96*sqrt(results$var)
  
  results$N <- mean(ate_M_df$N)

  all_results_joint_ls[[ej]] <- results 
}

all_results_joint <- data.frame(do.call(rbind, all_results_joint_ls)) 
rownames(all_results_joint) <- unlist(lapply(1:length(exposures_joint), function(x) paste0(exposures_joint[[x]], collapse=":")))



write.csv(all_results_joint, file="/Users/danagoin/Documents/Research projects/CiOB-ECHO/Projects/Stress and birth weight/results/R1/stress_bw_joint_R1.csv")

  # delta method SE by hand 
  # by hand 
  # want the sum of both exposures plus interaction 
  # so partial derivative of each x with respect to the sum is just 1 for each 
  #grad <- c(1,1,1)
  
  # covariance matrix 
  
  #cov <- vcov(fit)
  #cov <- cov[c(exposures_joint, paste0(exposures_joint, collapse=":")),c(exposures_joint, paste0(exposures_joint, collapse=":"))]
  
  #vG <- t(grad) %*% cov %*% grad
  
  #sqrt(vG)
  