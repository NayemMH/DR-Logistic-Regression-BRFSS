library(haven)
library(tidyverse)
library(tictoc)
library(ResourceSelection)  
library(pROC)               
library(speedglm)          
library(drglm)             
library(car)                
library(pryr)               
library(mice)               
library(gt)                 



brfss_cleaned <- brfss |>
  filter(across(
    c("DIABETE4","EXERANY2","_MICHD","ASTHMA3",
      "ADDEPEV3","CHCKDNY2","MARITAL","_EDUCAG",
      "DRNKANY6","_AGEG5YR","_BMI5CAT","_INCOMG1",
      "_RFSMOK3","GENHLTH","HAVARTH4","SEXVAR"),
    ~ !(.x %in% c(7, 9)))) |>
  filter(`_AGEG5YR` >= 5) |>
  mutate(across(
    c("DIABETE4","EXERANY2","_MICHD","ASTHMA3",
      "ADDEPEV3","CHCKDNY2","MARITAL","_EDUCAG",
      "DRNKANY6","_AGEG5YR","_BMI5CAT","_INCOMG1",
      "_RFSMOK3","GENHLTH","HAVARTH4","SEXVAR"),
    as.factor )) |>
  mutate(DIABETE4 = as.factor(case_when(
    DIABETE4 == "1" ~ "Yes",
    DIABETE4 == "3" ~ "No",
    TRUE ~ NA_character_ ))) |>
  drop_na()


set.seed(42)
brfss_cleaned <- brfss_cleaned[sample(nrow(brfss_cleaned)), ]

desc_overall <- brfss_cleaned |>
  summarise(
    n           = n(),
    pct_female  = mean(SEXVAR == "2") * 100,
    pct_bmi_ob  = mean(`_BMI5CAT` %in% c("3","4")) * 100,
    pct_smoker  = mean(`_RFSMOK3` == "1") * 100,
    pct_noex    = mean(EXERANY2 == "2") * 100,
    pct_cvd     = mean(`_MICHD` == "1") * 100,
    pct_fairpoor= mean(GENHLTH %in% c("4","5")) * 100,
    pct_college = mean(`_EDUCAG` == "4") * 100,
    pct_highinc = mean(`_INCOMG1` == "6") * 100
  )

desc_by_diab <- brfss_cleaned |>
  group_by(DIABETE4) |>
  summarise(
    n           = n(),
    pct_female  = mean(SEXVAR == "2") * 100,
    pct_bmi_ob  = mean(`_BMI5CAT` %in% c("3","4")) * 100,
    pct_smoker  = mean(`_RFSMOK3` == "1") * 100,
    pct_noex    = mean(EXERANY2 == "2") * 100,
    pct_cvd     = mean(`_MICHD` == "1") * 100,
    pct_fairpoor= mean(GENHLTH %in% c("4","5")) * 100,
    pct_college = mean(`_EDUCAG` == "4") * 100,
    pct_highinc = mean(`_INCOMG1` == "6") * 100
  )

print(desc_overall)
print(desc_by_diab)

# Chi-square / t-test p-values
pvals <- sapply(
  c("SEXVAR","_BMI5CAT","_RFSMOK3","EXERANY2",
    "_MICHD","GENHLTH","_EDUCAG","_INCOMG1"),
  function(v) chisq.test(table(brfss_cleaned[[v]],
                               brfss_cleaned$DIABETE4))$p.value
)
print(pvals)   


# glm()

mem_before_glm <- pryr::mem_used()
tic("glm()")
fit_glm <- glm(DIABETE4 ~ ., family = binomial(), data = brfss_cleaned)
time_glm <- toc(log = TRUE, quiet = TRUE)
mem_glm  <- pryr::mem_used() - mem_before_glm


# speedglm()

mem_before_sg <- pryr::mem_used()
tic("speedglm()")
fit_sg <- speedglm(DIABETE4 ~ ., family = binomial(), data = brfss_cleaned)
time_sg <- toc(log = TRUE, quiet = TRUE)
mem_sg  <- pryr::mem_used() - mem_before_sg


k_values <- c(25, 50, 100)
dnr_results <- list()

for (k in k_values) {
  mem_before_k <- pryr::mem_used()
  tic(paste0("D&R K=", k))
  fit_k <- drglm(
    DIABETE4 ~ .,
    family      = "binomial",
    k           = k,
    fitfunction = "glm",
    data        = brfss_cleaned
  )
  time_k <- toc(log = TRUE, quiet = TRUE)
  mem_k  <- pryr::mem_used() - mem_before_k
  
  dnr_results[[paste0("K", k)]] <- list(
    fit  = fit_k,
    time = time_k$toc - time_k$tic,
    mem  = mem_k
  )
  cat(sprintf("D&R K=%d | time: %.1f s | mem: %.2f GB\n",
              k, time_k$toc - time_k$tic, mem_k / 1e9))
}

coef_glm  <- summary(fit_glm)$coefficients
coef_sg   <- summary(fit_sg)$coefficients
coef_dnr50 <- dnr_results$K50$fit



# Hosmer-Lemeshow goodness of fit
hl_test <- hoslem.test(fit_glm$y, fitted(fit_glm), g = 10)
print(hl_test)   

# AUC / ROC 
phat     <- fitted(fit_glm)
roc_obj  <- roc(brfss_cleaned$DIABETE4, phat)
print(ci.auc(roc_obj))

# VIF
vif_vals <- car::vif(fit_glm)
print(vif_vals)


# Table tab:computation_brfss)

time_summary <- data.frame(
  Method          = c("glm()", "speedglm()",
                      "D&R (K=25)", "D&R (K=50)", "D&R (K=100)"),
  Time_sec        = c(time_glm$toc - time_glm$tic,
                      time_sg$toc  - time_sg$tic,
                      dnr_results$K25$time,
                      dnr_results$K50$time,
                      dnr_results$K100$time),
  Mem_GB          = c(mem_glm / 1e9,
                      mem_sg  / 1e9,
                      dnr_results$K25$mem / 1e9,
                      dnr_results$K50$mem / 1e9,
                      dnr_results$K100$mem / 1e9)
)
time_summary$Speedup <- round(time_summary$Time_sec[1] /
                                time_summary$Time_sec, 2)
time_summary$Mem_reduction <- round(
  (1 - time_summary$Mem_GB / time_summary$Mem_GB[1]) * 100, 0)
print(time_summary)


# SENSITIVITY ANALYSIS: multiple K values
k_sens <- c(10, 25, 50, 75, 100, 150)
sens_coefs <- lapply(k_sens, function(k) {
  fit_k <- drglm(
    DIABETE4 ~ ., family = "binomial", k = k,
    fitfunction = "glm", data = brfss_cleaned
  )
  coef(fit_k)
})
names(sens_coefs) <- paste0("K=", k_sens)

# Maximum deviation from K=50 baseline across all parameters
baseline <- sens_coefs[["K=50"]]
max_devs <- sapply(sens_coefs, function(c_k)
  max(abs(c_k - baseline))
)
cat("Max deviation from K=50 baseline:\n")
print(round(max_devs, 4))   


# MULTIPLE IMPUTATION SENSITIVITY ANALYSIS (mice on a random 500,000-obs subsample, 20 imputations)

set.seed(123)
idx_sub <- sample(nrow(brfss_cleaned), 500000)

# Raw (pre-factored) numeric subset for mice
brfss_raw_sub <- brfss |>
  filter(across(
    c("DIABETE4","EXERANY2","_MICHD","ASTHMA3",
      "ADDEPEV3","CHCKDNY2","MARITAL","_EDUCAG",
      "DRNKANY6","_AGEG5YR","_BMI5CAT","_INCOMG1",
      "_RFSMOK3","GENHLTH","HAVARTH4","SEXVAR"),
    ~ !(.x %in% c(7, 9))
  )) |>
  filter(`_AGEG5YR` >= 5) |>
  mutate(DIABETE4 = case_when(
    DIABETE4 == 1 ~ 1,
    DIABETE4 == 3 ~ 0,
    TRUE ~ NA_real_
  )) |>
  slice(idx_sub)

tic("Multiple imputation (mice)")
imp <- mice::mice(brfss_raw_sub, m = 20, maxit = 5,
                  method = "pmm", seed = 123, printFlag = FALSE)
toc()

fits_imp <- lapply(1:imp$m, function(i) {
  dat_i <- mice::complete(imp, i)
  glm(mi_formula, family = binomial(), data = dat_i)
})

pool_imp    <- mice::pool(fits_imp)
summary_imp <- summary(pool_imp)

# Compare to complete-case estimates on same subset
fit_cc_sub <- glm(DIABETE4 ~ ., family = binomial(),
                  data = brfss_cleaned[idx_sub, ])

mean_abs_diff <- mean(abs(
  coef(fit_cc_sub) - summary_imp$estimate
), na.rm = TRUE)
cat("Mean absolute difference (CC vs MI):", round(mean_abs_diff, 4),
    "log-odds units\n")  


# HOMOGENEITY CHECKS ACROSS D&R SUBSETS

n_total <- nrow(brfss_cleaned)
K_check <- 50
subset_ids <- ceiling(seq_len(n_total) / (n_total / K_check))

if ("IYEAR" %in% colnames(brfss_cleaned)) {
  year_props_overall <- prop.table(table(brfss_cleaned$IYEAR))
  year_props_by_sub  <- tapply(
    brfss_cleaned$IYEAR, subset_ids,
    function(y) prop.table(table(y))
  )
  max_temporal_dev <- max(sapply(year_props_by_sub, function(p) {
    common <- intersect(names(p), names(year_props_overall))
    max(abs(p[common] - year_props_overall[common]))
  }))
  cat("Max temporal deviation across subsets:", max_temporal_dev, "\n")
  cat("Within 2%:", max_temporal_dev < 0.02, "\n")
}

# Demographic homogeneity: sex distribution per subset
sex_overall <- prop.table(table(brfss_cleaned$SEXVAR))
sex_by_sub  <- tapply(brfss_cleaned$SEXVAR, subset_ids,
                      function(x) prop.table(table(x)))
max_demo_dev <- max(sapply(sex_by_sub, function(p) {
  common <- intersect(names(p), names(sex_overall))
  max(abs(p[common] - sex_overall[common]))
}))
cat("Max demographic (sex) deviation across subsets:", max_demo_dev, "\n")
cat("Within 2%:", max_demo_dev < 0.02, "\n")

