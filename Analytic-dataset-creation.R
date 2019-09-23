
rm(list=ls())

library(dplyr)
library(readxl)
library(ggplot2)
library(stringr)
library(tidyverse)
library(haven)

# prepare data for multiple dimensions of stress and fetal growth analysis 

# read in CiOB survey data

df <- read.csv("/Users/danagoin/Documents/Research projects/CiOB-ECHO/Data/CiOB2 data/questionnaire.csv")
# this will be read in later with stress measures coded by Stephanie 
df <- df  %>% select(-cmnty_ladder)

# read in biospecimen collection data 

df_b <- read.csv("/Users/danagoin/Documents/Research projects/CiOB-ECHO/Data/CiOB2 data/biospecimen.csv")
df_b <- df_b %>% select(ppt_id, secondtri_date)


# read in medical record abstraction data 

df_mr <- read.csv("/Users/danagoin/Documents/Research projects/CiOB-ECHO/Data/CiOB2 data/medicalrecordabstraction.csv")

# fetal growth created variables 
# this data set came from Stephanie, to be consistent with how fetal growth measures have previously been defined 
df_f <- read.csv("/Users/danagoin/Documents/Research projects/CiOB-ECHO/Data/CiOB2 data/fetal growth.csv")

df_f <- df_f %>% select(ID, ga_weeks_comb, SGA_10, baby_wt_grams, lbw_cat, ga_cat.4)
df_f$ppt_id <- df_f$ID

df_f$ptb <- ifelse(df_f$ga_cat.4=="Preterm",1,ifelse(df_f$ga_cat.4=="Full Term",0,NA))
df_f$sga <- ifelse(df_f$SGA_10=="1: Small for GA",1, ifelse(df_f$SGA_10=="0: Normal for GA",0,NA))

df <- left_join(df, df_f)
df <- left_join(df, df_b)
df <- left_join(df, df_mr)

# create gestational age at time of second trimester visit 

df$gestwks_secondtri <- round(difftime(strptime(df$secondtri_date, format="%Y-%m-%d"), strptime(df$lmp_dt_mr, format="%Y-%m-%d"), units = "weeks"),0)
# ask if it's possible that people were seen before 12 or after 28 weeks for their second trimester visit 
df$gestwks_secondtri <- ifelse(df$gestwks_secondtri<12 | df$gestwks_secondtri>28, NA, df$gestwks_secondtri)
df$gestwks_secondtri <- ifelse(is.na(df$gestwks_secondtri), round(df$ga_weeks_comb - floor(difftime(strptime(df$dob_c_mr, format="%Y-%m-%d"), strptime(df$secondtri_date, format="%Y-%m-%d"), units = "weeks")), 0), df$gestwks_secondtri)

summary(df$gestwks_secondtri)
#View(df %>% select(ppt_id, lmp_dt_mr, secondtri_date, dob_c_mr, gestwks_secondtri))

# recode key variables  

df$mat_edu <- ifelse(df$edu_m>=97, NA, df$edu_m)                   
df$mat_edu <- factor(df$mat_edu, levels=c(0,1,2,3,4,5), 
                     labels=c("Less than high school","High school grad","Some college","College grad","Master's degree","Doctoral degree"))


df$mat_race <- ifelse(df$race_m>=97,NA, df$race_m)

df$mat_eth <- ifelse(df$latina_m>=97,NA, df$latina_m)

df$mat_race_eth <- ifelse(df$mat_eth==1, 1, 
                          ifelse(df$mat_eth==0 & (df$mat_race==1 | df$mat_race==2), 2, 
                                 ifelse(df$mat_eth==0 & df$mat_race==3, 3, 
                                        ifelse(df$mat_eth==0 & df$mat_race==4, 4, 
                                               ifelse(df$mat_eth==0 & (df$mat_race==5 | df$mat_race==6), 5, NA)))))

df$mat_race_eth <- factor(df$mat_race_eth, levels=c(1,2,3,4,5), 
                          labels=c("Latina","Asian/PI", "Black","White","Other or multiple"))


df$marital <- ifelse(df$marital_stat>=97,NA, df$marital_stat)

df$marital <- ifelse(df$marital==4, 3, df$marital)
df$marital <- factor(df$marital, levels=c(1,3,5), labels=c("Married","Widowed, separated, or divorced","Never married"))

df$medi_cal <- ifelse(df$medi_cal_m>=97,NA, df$medi_cal_m)

df$parity <- ifelse(df$parity_mr==9, NA, ifelse(df$parity_mr>3, 3, df$parity_mr))
df$parity <- factor(df$parity, levels=c(0,1,2,3), labels=c("0","1","2","3+"))

df$num_child_at_home <- ifelse(df$child_n==0 & !is.na(df$child_n), 0, 
                               ifelse(df$child_live_n>3 & !is.na(df$child_live_n), 3, df$child_live_n))

# if missing child_live_n but not missing child_n, replace with value from child_n 
df$num_child_at_home <- ifelse(is.na(df$num_child_at_home) & !is.na(df$child_n) & df$child_n>3 & df$child_n<95, 3, 
                               ifelse(is.na(df$num_child_at_home) & df$child_n==99, NA, 
                                      ifelse(is.na(df$num_child_at_home) & !is.na(df$child_n), df$child_n, df$num_child_at_home)))

df$num_child_at_home <- factor(df$num_child_at_home, levels=c(0,1,2,3), labels=c("0","1","2","3+"))

# infant sex: 0 = female, 1 = male 
df$infant_sex <- ifelse(df$sex_c_mr==1, 1, 
                        ifelse(df$sex_c_mr==0, 0, NA))
df$infant_sex <- factor(df$infant_sex, levels=c(0,1), labels = c("Female","Male"))

# recode variables you need 
# age at delivery from medical record
df$mat_age <- ifelse(df$age_dlvry_mr==9999, NA, df$age_dlvry_mr)
# if missing, use mother's date of birth and date of child's birth from the medical record and round to nearest year 
df$mat_age <- ifelse(is.na(df$mat_age), 
                     round(difftime(strptime(df$dob_c_mr, format="%Y-%m-%d"), strptime(df$dob_m_mr, format="%Y-%m-%d"), units = "days")/365.25,0), df$mat_age)
# if still missing, use mother's date of birth and date of second trimester visit and round to nearest year 
df$mat_age <- ifelse(is.na(df$mat_age), 
                     round(difftime(strptime(df$secondtri_date, format="%Y-%m-%d"), strptime(df$dob_m_mr, format="%Y-%m-%d"), units = "days")/365.25,0), df$mat_age)
# should have no missing after this 
sum(is.na(df$mat_age))

# create income categories 

df$income_hh <- ifelse(df$income_hh %in% c(97,98,99), NA, df$income_hh)
df$income_20k <- ifelse(df$income_20k %in% c(97,98,99), NA, df$income_20k)
df$income_40k <- ifelse(df$income_40k %in% c(97,98,99), NA, df$income_40k)

df$hh_income_cat1 <- ifelse(df$income_hh==1,1, 
                            ifelse(df$income_hh==2,1,
                                   ifelse(df$income_hh==3, 1, 
                                          ifelse(df$income_hh==4,1,
                                                 ifelse(df$income_hh==5,1, 
                                                        ifelse(df$income_hh==6,1, 
                                                               ifelse(df$income_hh==7,1,
                                                                      ifelse(df$income_hh==8,1, 0))))))))

df$hh_income_cat1 <-  ifelse(((df$income_20k==0 & !is.na(df$income_20k)) | (df$income_40k==0 & !is.na(df$income_40k))), 1, df$hh_income_cat1)



df$hh_income_cat2 <- ifelse(df$income_hh==9,2, 
                            ifelse(df$income_hh==10,2, 
                                   ifelse(df$income_hh==11,2, 
                                          ifelse(df$income_hh==12,2, 
                                                 ifelse(df$income_hh==13,2,
                                                        ifelse(df$income_hh==14,2,
                                                               ifelse(df$income_hh==15,2, 0)))))))

df$hh_income_cat2 <- ifelse((df$income_80k==0 & !is.na(df$income_80k)), 2, df$hh_income_cat2)


df$hh_income_cat3 <-  ifelse(df$income_hh==16,3, 
                             ifelse(df$income_hh==17,3, 
                                    ifelse(df$income_hh==18,3,
                                           ifelse(df$income_hh==19,3,
                                                  ifelse(df$income_hh==20,3,0)))))

df$hh_income_cat <- ifelse(df$hh_income_cat1==1,1, 
                           ifelse(df$hh_income_cat2==2, 2, 
                                  ifelse(df$hh_income_cat3==3, 3, NA)))



df$hh_income_cat <- factor(df$hh_income_cat, levels=c(1,2,3), labels=c("<40,000","$40,000-$79,999","$80,000+"))

# create bmi measures
# weight
df$pp_weight_lbs <- ifelse(df$wt_preprg_mr == 9999, NA, 
                             ifelse(df$wt_preprg_mr<80, NA, df$wt_preprg_mr))
# what are weights that are infeasible? less than 80 lbs? 


# height 
# recode 9999 to 0 for feet, because this corresponds usually to height being recorded in inches only
df$height_ft_mr <- ifelse(df$height_ft_mr==9999, 0, df$height_ft_mr)

# recode 9999 to NA for inches, because this corresponds to missing both feet and inches
df$height_in_mr <- ifelse(df$height_in_mr==9999, NA, df$height_in_mr)

# change height_mr to character to enable regular expressions / string operations 
df$height_mr <- as.character(df$height_mr)

# if the values for height in inches and feet are missing, replace with values from height_mr 
df$height_ft_mr <- ifelse(is.na(df$height_ft_mr) & str_detect(df$height_mr,"([0-9]')"), str_extract(str_extract(df$height_mr, "([0-9]')"), "([0-9])"), 
                            ifelse(is.na(df$height_ft_mr) & str_detect(df$height_mr, "([0-9] ft)"), str_extract(str_extract(df$height_mr,"([0-9] ft.)"), "([0-9])"), df$height_ft_mr))

df$height_in_mr <- ifelse(is.na(df$height_in_mr) & str_detect(df$height_mr,"([0-9]\")"), str_extract(str_extract(df$height_mr, "([0-9]\")"), "([0-9])"), 
                            ifelse(is.na(df$height_in_mr) & str_detect(df$height_mr,"([0-9][0-9]\")"), str_extract(str_extract(df$height_mr, "([0-9][0-9]\")"), "([0-9][0-9])"), 
                                   ifelse(is.na(df$height_in_mr) & str_detect(df$height_mr,"([0-9]'')"), str_extract(str_extract(df$height_mr, "([0-9]'')"), "([0-9])"), 
                                          ifelse(is.na(df$height_in_mr) & str_detect(df$height_mr,"([0-9][0-9]'')"), str_extract(str_extract(df$height_mr, "([0-9][0-9]'')"), "([0-9][0-9])"), 
                                                 ifelse(is.na(df$height_in_mr) & str_detect(df$height_mr, "([0-9] in)"), str_extract(str_extract(df$height_mr,"([0-9] in)"), "([0-9])"), 
                                                        ifelse(is.na(df$height_in_mr) & str_detect(df$height_mr, "([0-9][0-9] in)"), str_extract(str_extract(df$height_mr,"([0-9][0-9] in)"), "([0-9][0-9])"), df$height_in_mr))))))


# calculate height in inches 
df$height_ft_mr <- as.numeric(df$height_ft_mr)
df$height_in_mr <- as.numeric(df$height_in_mr)

df <- df %>% mutate(pp_height_in = height_ft_mr*12 + height_in_mr)

df$height_num <- ifelse(df$height_mr=="9999",NA,as.numeric(df$height_mr))

# for those recorded height in cm, capture and convert to inches
df$height_mr_cm <- ifelse(grepl("cm",df$height_mr, ignore.case = T), str_replace(df$height_mr,"cm",""), NA)
# anything greater than 120 centimeters, which is just less than 4 ft, is likely recorded in cm 
df$height_mr_cm <- ifelse(is.na(df$height_mr_cm) & df$height_num>120, df$height_num, df$height_mr_cm)
df$height_mr_cm <- as.numeric(df$height_mr_cm)
# height less than 120 cm is implausible, recode to missing
df$height_mr_cm <- ifelse(df$height_mr_cm<120, NA, df$height_mr_cm)  


df$pp_height_in <- ifelse(is.na(df$pp_height_in), df$height_mr_cm/2.54, df$pp_height_in)

# recode prepregnancy BMI if missing from height and weight calculated variables 
df$pp_bmi <- ifelse(df$bmi_preprg_mr==9999,NA,df$bmi_preprg_mr)
df$pp_bmi <- ifelse(is.na(df$pp_bmi), 703*df$pp_weight_lbs/(df$pp_height_in^2), df$pp_bmi)



# merge on measures of stress 

df_s <- read_sas("/Users/danagoin/Documents/Research projects/CiOB-ECHO/Projects/Stress and birth weight/data/stress_clean.sas7bdat")

df_s <- df_s %>% select(ppt_id, cmnty_ladder, perceived_stress, stressful_life_events_missing, cesd_cat, econ_cat, 
                        job_highstrain_missing, discrimination, unplanned_pregnancy, caregiving_missing, food_cat, 
                        collective_efficacy_cat, neighborhood_safety_cat, neighborhood_dissatisfaction_cat, neighborhood_disorder_cat)

# missing food_cat_labs --- just use food_cat instead
# join with rest of data 
df_s$ppt_id <- as.integer(df_s$ppt_id)

df <- left_join(df, df_s) 

# financial strain -- econ_cat 
# food insecurity -- food_cat
# high strain job -- job_highstrain_missing
# poor perceived neighborhood quality -- collective_efficacy_cat ==1 or neighborhood_safety_cat==1 | neighborhood_dissatisfaction_cat==1 | neighborhood_disorder_cat==1 
# unplanned pregnancy -- unplanned_pregnancy 
# caregiving -- caregiving_missing 
# low perceived social status -- cmty_ladder
# high perceived stress -- perceived_stress 
# experieincing discrimination -- discrimination 
# stressful life events -- stressful_life_events_missing 
# depressive symptoms -- cesd_cat -- 2 or 3 points on 6 or more questions  

# will need to categorize these
# then in imputation, will need to make sure that the cutoff is coded correctly 

# create poor perceived neighborhood quality indicator 
# if any of the indicators are 1, code as 1 
# if all are 0, code as 0 
# if any are missing but some are 1, code as 1 
# if any are missing but rest are 0, code as NA 
df$poor_perceived_nhbd <- ifelse((is.na(df$collective_efficacy_cat) | is.na(df$neighborhood_safety_cat) | is.na(df$neighborhood_dissatisfaction_cat) | is.na(df$neighborhood_disorder_cat))
                                 & rowSums(df %>% select(collective_efficacy_cat, neighborhood_safety_cat, neighborhood_dissatisfaction_cat, neighborhood_disorder_cat), na.rm=T)==0, NA, 
                              ifelse(df$collective_efficacy_cat==1  | df$neighborhood_safety_cat==1  | df$neighborhood_dissatisfaction_cat==1 | df$neighborhood_disorder_cat==1, 1, 0))


# low perceived social status 
# if they rated themself below a 5, coded as 1 
# if they rated themself above or equal to 5, code as 0 
df$low_perceived_socialstat <- ifelse(df$cmnty_ladder<5, 1, ifelse(df$cmnty_ladder>=5, 0, NA))

# high perceived stress 
# cutoff of >=9 to denote a high level of perceived stress 
df$high_perceived_stress <- ifelse(df$perceived_stress>=9,1, ifelse(df$perceived_stress<9, 0, NA))

# stressful life events 
# 2 or more stressful life events were considered a high level of stress
# may want to reconsider the cutoff to 3 -- many people have experienced 2 
df$high_stressful_events <- ifelse(df$stressful_life_events_missing>=2, 1, ifelse(df$stressful_life_events_missing<2,0, NA))



df_m <- df %>% select(ppt_id, ptb, sga, ga_weeks_comb, baby_wt_grams, econ_cat, food_cat, job_highstrain_missing, poor_perceived_nhbd, 
                      unplanned_pregnancy, caregiving_missing, low_perceived_socialstat, high_perceived_stress, discrimination, high_stressful_events, cesd_cat, 
                      mat_age, mat_edu, mat_race_eth, hh_income_cat, medi_cal, marital, parity, pp_bmi, infant_sex)




save(df_m, file="/Users/danagoin/Documents/Research projects/CiOB-ECHO/Projects/Stress and birth weight/data/stress_fetal_growth_complete_case")

library(mice)

x <- md.pattern(df_m)
tail(x)

imp0 <- mice(df_m,  seed=978674653, max=0)
meth <- imp0$method
meth["ptb"] <- "~I(ga_weeks_comb<37)"

imp <- mice(df_m, seed=978674653, maxit=20, m=30, meth=meth)

# check that ga_weeks_comb and ptb have the correct relationship 
mice::complete(imp,3)[is.na(df_m$ptb), ]
#imp$imp$ptb


plot(imp, c("ptb","sga", "ga_weeks_comb"))
plot(imp, c("baby_wt_grams","mat_race_eth"))
plot(imp, c("mat_edu", "hh_income_cat","medi_cal"))
plot(imp, c("marital","parity","pp_bmi"))
plot(imp, c("infant_sex", "econ_cat", "food_cat"))
plot(imp, c("job_highstrain_missing", "poor_perceived_nhbd", "unplanned_pregnancy"))
plot(imp, c("caregiving_missing", "low_perceived_socialstat", "high_perceived_stress"))
plot(imp, c("discrimination", "high_stressful_events", "cesd_cat"))

#stripplot(imp, pch=20, cex=1.2)
#xyplot(imp, mat_race_eth ~ mat_edu, subset =.imp==1, cex=1.2, pch=20, jitter.x=T)

densityplot(imp, scales=list(x=list(relation="free")))
densityplot(imp, ~ptb + sga + pp_bmi)
densityplot(imp, ~pp_bmi)

# complete cases data set 
dfcc <- df_m[complete.cases(df_m),] 

dfcc  <- cbind(dfcc, model.matrix(~mat_race_eth + mat_edu + hh_income_cat + marital + parity + infant_sex - 1, data=dfcc))

sumstats <- function(var) {
  if(class(with(dfcc, get(var)))=="numeric") {
    x <- dfcc %>% summarise(N = sum(!is.na(get(var))), min = min(get(var)), max=max(get(var)), median=median(get(var)), mean=mean(get(var)), sd = sqrt(var(get(var))))
    rownames(x) <- var
  }
  else {
    x <- NULL
  }
  return(x)
}
ls <- lapply(1:length(colnames(dfcc)), function(x) sumstats(colnames(dfcc)[x]))
table <- do.call(rbind, ls)





# for imputed data 
imp.data <- mice::complete(imp, "long")
imp.data  <- cbind(imp.data, model.matrix(~mat_race_eth + mat_edu + hh_income_cat + marital + parity + infant_sex - 1, data=imp.data))

isumstats <- function(var) {
  if(class(with(imp.data, get(var)))=="numeric") {
    x <- data.frame(t(colMeans(imp.data %>% group_by(.imp) %>% summarise(N = sum(!is.na(get(var))), min = min(get(var)), max=max(get(var)), median=median(get(var)), mean=mean(get(var)), sd = sqrt(var(get(var)))))))
    rownames(x) <- var
  }
  else {
    x <- NULL
  }
  return(x)
}

ils <- lapply(1:length(colnames(imp.data)), function(x) isumstats(colnames(imp.data)[x]))
itable <- do.call(rbind, ils)

# by imputed versus non-imputed values 

#MM <- data.frame(model.matrix(~mat_race_eth + mat_edu + hh_income_cat + marital + parity + infant_sex - 1, data=df_m))
#MM <- MM[match(rownames(df_m),rownames(MM)),]
# find which variables are missing, replace others with their observed values 
#MM[!is.na(df_m$infant_sex),] <- model.matrix(~mat_race_eth + mat_edu + hh_income_cat + marital + parity ~ , data=df_m)

#MM[,"b"] <- dat$b
#rownames(MM) <- rownames(dat)

#df_t  <- cbind(df_m, model.matrix(~mat_race_eth + mat_edu + hh_income_cat + marital + parity + infant_sex - 1, data=df_m))

# calculate descriptive statistics separately for each variable 
msumstats <- function(var) {
  if(class(with(df_m, get(var)))=="numeric") {
    x <- data.frame(t(colMeans(df_m %>% summarise(N = sum(!is.na(get(var))), min = min(get(var), na.rm=T), max=max(get(var), na.rm=T), median=median(get(var), na.rm=T), mean=mean(get(var), na.rm=T), sd = sqrt(var(get(var), na.rm=T))))))
    rownames(x) <- var
  }
  else {
    x <- NULL
  }
  return(x)
}

mls <- lapply(1:length(colnames(df_m)), function(x) msumstats(colnames(df_m)[x]))
mtable <- do.call(rbind, mls)


# save imputed data set 
save(imp, file="/Users/danagoin/Documents/Research projects/CiOB-ECHO/Projects/Stress and birth weight/data/stress_fetal_growth_imp")

