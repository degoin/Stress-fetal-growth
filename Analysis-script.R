
# stress and fetal growth analysis 

load("/Users/danagoin/Documents/Research projects/CiOB-ECHO/Projects/Stress and birth weight/data/stress_fetal_growth_imp")

stress_measures <- c("econ_cat", "food_cat", "job_highstrain_missing", "poor_perceived_nhbd", 
                     "low_perceived_socialstat", "discrimination",  "high_perceived_stress", 
                     "caregiving_missing", "high_stressful_events", "cesd_cat", "unplanned_pregnancy")


outcome_measures <- c("sga, baby_wt_grams")



bivar <- function(outcome, exposure) {
x <- summary(pool(with(imp, glm(get(outcome) ~ get(exposure)))))["get(exposure)",c("estimate","std.error","statistic","p.value")]
return(x)
}

ls <- lapply(1:length(stress_measures), function(x) bivar("sga",stress_measures[x]))

results <- data.frame(do.call(rbind, ls))
rownames(results) <- stress_measures 
