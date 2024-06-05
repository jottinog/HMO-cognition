# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Jonatan Ottino González for GitHub (May 29th 2024)                                  #
# Manuscript Number: AJCN-D-24-00564                                                  #
# Title: Consumption of Different Combinations of Human Milk Oligosaccharides         #
# in the First 6 Months of Infancy is Positively Associated with Early Cognition      #
# at 2 Years of Age in Latino Children                                                #  
# Journal: The American Journal of Clinical Nutrition                                 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#0# Load packages
pacman::p_load(tidyverse, performance, data.table, plsRglm, summarytools, vtable, effectsize)

#1# Set paths, covariates (vi) and load data
datadir <- 'Users/path/to/data'

vi <- c('mom_age', 'ses', 'prepreg_bmi_kgm2', 'mode_of_delivery', 'bf_frequency_1month', 'child_sex')     

df <- fread(paste0(datadir, '/projects/nutrition/df_analysis.csv')) %>% 
  select(., study_id, matches(vi), everything()) %>% 
  as.data.frame() 

st(df, vars = c(vi, 'child_age', 'secretor_status', 'bsid_cog_ss'))     # table of summary statistics

#2# Prepare X and Y data and perform PLS
response_variable <- df$bsid_cog_ss    # define response variable vector (Y)

predictor_matrix <- df %>% select(., ends_with('nmol_ml_1month')) %>% mutate_all(., ~as.numeric(.))    # save predictors

model <- plsR(response_variable, predictor_matrix, scaleX = T, scaleY = T, nt = dim(predictor_matrix)[2], typeVC = 'none')

#3# Test intercept-models against HMO components 
for (i in 1:dim(model$tt)[2]) {
  df[paste0('comp', i)] <- model$tt[,i]
}

# Function to validate HMO components
get_model_metrics <- function(component, data) {
  # Create the full model with the specified component
  full_formula <- as.formula(paste("response_variable ~", component))
  full_model <- lm(full_formula, data = data)
  # Create the intercept-only model
  intercept_formula <- as.formula("response_variable ~ 1")
  intercept_model <- lm(intercept_formula, data = data)
  # Compare the full model to the intercept-only model
  anova_result <- anova(intercept_model, full_model)
  # Extract F-statistic, p-value, and Residual Standard Error
  f_statistic <- round(anova_result$F[2], 3)
  dof <- anova_result[2,1] 
  p_value <- round(anova_result$Pr[2], 4)
  # Extract the Residual Standard Error from the full model
  rse <- summary(full_model)$sigma
  # Convert RSE to RMSE
  n <- length(full_model$residuals)
  rmse <- round(sqrt(mean(full_model$residuals^2)), 3)
  # Calculate R-squared
  r_squared <- round(summary(full_model)$r.squared, 3)
  # Return a data frame with the results
  return(data.frame(
    comp.name = component,
    F_statistic = f_statistic,
    dof = dof,
    rmse = rmse,
    r_squared = r_squared,
    p_value = p_value
  ))
}

predictor_variables <- paste0("comp", 1:19)     # generate character vector w/comp names

result_df <- data.frame()

for (predictor in predictor_variables) {
  # Call the function to validate HMO components and bind the result to the results data frame
  result <- get_model_metrics(predictor, df)
  result_df <- rbind(result_df, result)
}

print(head(result_df))     # check comparison of intercept-only models vs HMO components

#4# append to the nuisance covariate list (vi) significanat HMO components
(comp.vi <- c(vi, result_df$comp.name[result_df$p_value < 0.05]))  

#5# Variable selection w/bootstrap (optimization of HMO component)
set.seed(4994)

boot.model.1m <- bootpls(model, typeboot = 'plsmodel', R = 250)       #  bootstrapped mean weight of each feature across all components

plots.confints.bootpls(confints.bootpls(boot.model.1m, indice = 1:dim(predictor_matrix)[2] + 1), type = 'BCa', las = 3, mar = c(12,1,1,1)) # pick features that whose interval does not include/pass through 0

data <- as.data.frame(confints.bootpls(boot.model.1m, indice = 2:20)[,7:8]) %>% 
  rownames_to_column(., var = 'Variable') %>% 
  rename(Lower_Limit = V1, Upper_Limit = V2)

check_confidence_intervals <- function(data) {
  variables_without_zero <- character()
  for (i in 1:nrow(data)) {
    variable <- data$Variable[i]
    lower_limit <- data$Lower_Limit[i]
    upper_limit <- data$Upper_Limit[i]
    if (lower_limit > 0 || upper_limit < 0) {
      variables_without_zero <- c(variables_without_zero, variable)
    }
  }
  return(variables_without_zero)
}

(result <- check_confidence_intervals(data))

head(boot.predictor.matrix <- predictor_matrix %>% select(., all_of(result)))

(boot.model.1m <- plsR(response_variable, boot.predictor.matrix, scaleX = T, scaleY = T, nt = dim(predictor_matrix)[2], typeVC = 'none'))

(result <- names(boot.model.1m$ww[,1])[abs(boot.model.1m$ww[,1]) > 0.10])      # keep only > .10 loadings 

head(boot.predictor.matrix <- predictor_matrix %>% select(., all_of(result)))

(boot.model.1m <- plsR(response_variable, boot.predictor.matrix, scaleX = T, scaleY = T, nt = dim(boot.predictor.matrix)[2], typeVC = 'none'))

for (i in 1:dim(boot.model.1m$tt)[2]) {
  df[paste0('boot.comp', i)] <- boot.model.1m$tt[,i]
}

(predictor_variables <- names(df %>% select(starts_with('boot.comp'))))

result_df <- data.frame()

for (predictor in predictor_variables) {
  # Call the function and bind the result to the results data frame
  result <- get_model_metrics(predictor, df)
  result_df <- rbind(result_df, result)
}

print(result_df)

(boot.comp.vi <- c(vi, result_df$comp.name[result_df$p_value < 0.05]))   # specific list of covariates (vi) + boot.comp for later supplementary analysis

#~~~~~~~~~~~~~~#
#~~NULL MODEL~~#     #7# HMO components already generated. Now test first null model performance.
#~~~~~~~~~~~~~~#

summary(null_model <- lm(as.formula(paste('bsid_cog_ss', paste(vi, collapse = '+'), sep = '~')), df))  
olsrr::ols_test_normality(null_model)   # check normality of residuals
olsrr::ols_test_breusch_pagan(null_model)   # check heteroskedasticity
any(olsrr::ols_coll_diag(null_model)$vif_t$VIF > 2)   # no VIF > 2 (multicollinearity)

#6# function to look, flag and remove for outliers if normality or heteroskedasticity assumptions are violated
if(olsrr::ols_test_outlier(null_model)$bonferroni_p_val < 0.05 &
   olsrr::ols_test_normality(null_model)[[1]]$p.value < 0.05 | 
   olsrr::ols_test_breusch_pagan(null_model)[[4]] < 0.05) {
  outliers <- as.numeric(rownames(olsrr::ols_test_outlier(null_model)))  
  df.no.outlier <- df[-outliers, ]         # if outlier present, df.no.outlier is a version of df without outlier
  paste0('outlier present: ', outliers)
} else {
  df.no.outlier <- df                      # if no outlier present, df.no.outlier == df
  paste0('no outliers')
}

#7# re-run null model in new df (with no outlier) and check assumptions again
summary(null_model.new <- lm(as.formula(paste('bsid_cog_ss', paste(vi, collapse = '+'), sep = '~')), df.no.outlier))
olsrr::ols_test_normality(null_model.new)
olsrr::ols_test_breusch_pagan(null_model.new)
sqrt(mean(null_model.new$residuals^2))              # RMSE
confint(null_model.new)                             # compute CI 95% for all covariates in null model for reporting

# # # # # #
# MODEL 1 #     #11# fit/compare final model (only difference with null model is the addition of HMO comps)
# # # # # #

anova(null_model.new, model1 <- lm(as.formula(paste('bsid_cog_ss', paste(setdiff(comp.vi, 'comp2'), collapse = '+'), sep = '~')), df.no.outlier))  # test comp1
cohens_f(null_model.new, model2 = model1, alternative = 'two.sided')   # compute effect size
(p.val1 <- as.numeric(anova(null_model.new, model1)[[6]][2]))   # store p-value for later adjustement

olsrr::ols_test_normality(model1)            # residuals are normal
olsrr::ols_test_breusch_pagan(model1)        # homoskedasticity is normal
any(olsrr::ols_coll_diag(model1)$vif_t$VIF > 2)   # no VIF > 2

summary(model1)                         # assess the behavior of each covariate in the model
sqrt(mean(model1$residuals^2))          # RMSE
summary(model1)[[9]]                    # get the adj-R2
confint(model1, 'comp1')                # compute 95% CI of comp1
car::avPlot(model1, 'comp1')

# # # # # #
# MODEL 2 #
# # # # # #

anova(null_model.new, model2 <- lm(as.formula(paste('bsid_cog_ss', paste(setdiff(comp.vi, 'comp1'), collapse = '+'), sep = '~')), df.no.outlier))  # test comp2
cohens_f(null_model.new, model2 = model2, alternative = 'two.sided')
(p.val2 <- as.numeric(anova(null_model.new, model2)[[6]][2]))

olsrr::ols_test_normality(model2)            
olsrr::ols_test_breusch_pagan(model2)        
any(olsrr::ols_coll_diag(model2)$vif_t$VIF > 2)   

summary(model2)
sqrt(mean(model2$residuals^2))         
summary(model2)[[9]]
confint(model2, 'comp2')
car::avPlot(model2, 'comp2')

# Select the best performing model
models <- list(model1, model2)

results_model <- as.data.frame(matrix(nrow = length(models), ncol = 2))

colnames(results_model) <- c('comp', 'pval')

for (i in seq_along(models)) {
  results_model[i,1] <- setdiff(names(models[[i]]$model), names(null_model$model))
  results_model[i,2] <- anova(null_model.new, models[[i]])$Pr[2]
}

(best_comp <- results_model$comp[which.min(results_model$pval)])

# # # # # # # #
#  BOOT MODEL #
# # # # # # # #

anova(null_model.new, boot.final_model <- lm(as.formula(paste('bsid_cog_ss', paste(boot.comp.vi, collapse = '+'), sep = '~')), df.no.outlier))  
(p.val3 <- as.numeric(anova(null_model.new, boot.final_model)[[6]][2]))
cohens_f(null_model.new, model2 = boot.final_model, alternative = 'two.sided')

olsrr::ols_test_normality(boot.final_model)
olsrr::ols_test_breusch_pagan(boot.final_model)
any(olsrr::ols_coll_diag(boot.final_model)$vif_t$VIF > 2)   

summary(boot.final_model)
sqrt(mean(boot.final_model$residuals^2))          
summary(boot.final_model)[[9]]
confint(boot.final_model, 'boot.comp1')

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~ 6 - M O N T H  A N A L Y S I S #~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

(vi2 <- c(vi, 'bf_frequency_6month', best_comp, 'boot.comp1'))

df2 <- df %>% 
  select(., vi2, child_age, secretor_status, bsid_cog_ss, ends_with('nmol_ml_6month')) %>% 
  drop_na()

st(df2, vars = c(vi2, 'child_age', 'secretor_status', 'bsid_cog_ss'))     

#8# Prepare X and Y data and perform PLS
response_variable <- df2$bsid_cog_ss    

predictor_matrix <- df2 %>% select(., ends_with('nmol_ml_6month')) %>% mutate_all(., ~as.numeric(.))    

model <- plsR(response_variable, predictor_matrix, scaleX = T, scaleY = T, nt = dim(predictor_matrix)[2], typeVC = 'none')

#9# Test intercept-models against HMO components 
for (i in 1:dim(model$tt)[2]) {
  df2[paste0('comp_6m', i)] <- model$tt[,i]
}

predictor_variables <- paste0("comp_6m", 1:19)    

result_df <- data.frame()

for (predictor in predictor_variables) {
  # Call the function to validate HMO components and bind the result to the results data frame
  result <- get_model_metrics(predictor, df2)
  result_df <- rbind(result_df, result)
}

print(head(result_df))     

#10# append to the nuisance covariate list (vi) significanat HMO components
(comp.vi2 <- c(vi2[vi2 != 'boot.comp1'], result_df$comp.name[result_df$p_value < 0.05]))  

#11# Variable selection w/bootstrap (optimization of HMO component)
set.seed(4994)

boot.model.6m <- bootpls(model, typeboot = 'plsmodel', R = 250)       

plots.confints.bootpls(confints.bootpls(boot.model.6m, indice = 1:dim(predictor_matrix)[2] + 1), type = 'BCa', las = 3, mar = c(12,1,1,1)) 

data <- as.data.frame(confints.bootpls(boot.model.6m, indice = 2:20)[,7:8]) %>% 
  rownames_to_column(., var = 'Variable') %>% 
  rename(Lower_Limit = V1, Upper_Limit = V2)

check_confidence_intervals <- function(data) {
  variables_without_zero <- character()
  for (i in 1:nrow(data)) {
    variable <- data$Variable[i]
    lower_limit <- data$Lower_Limit[i]
    upper_limit <- data$Upper_Limit[i]
    if (lower_limit > 0 || upper_limit < 0) {
      variables_without_zero <- c(variables_without_zero, variable)
    }
  }
  return(variables_without_zero)
}

(result <- check_confidence_intervals(data))

head(boot.predictor.matrix <- predictor_matrix %>% select(., all_of(result)))

(boot.model.6m <- plsR(response_variable, boot.predictor.matrix, scaleX = T, scaleY = T, nt = dim(predictor_matrix)[2], typeVC = 'none'))

(result <- names(boot.model.6m$ww[,1])[abs(boot.model.6m$ww[,1]) > 0.10])      

head(boot.predictor.matrix <- predictor_matrix %>% select(., all_of(result)))

(boot.model.6m <- plsR(response_variable, boot.predictor.matrix, scaleX = T, scaleY = T, nt = dim(boot.predictor.matrix)[2], typeVC = 'none'))

for (i in 1:dim(boot.model.6m$tt)[2]) {
  df2[paste0('boot.comp_6m', i)] <- boot.model.6m$tt[,i]
}

(predictor_variables <- names(df2 %>% select(starts_with('boot.comp_6m'))))

result_df <- data.frame()

for (predictor in predictor_variables) {
  # Call the function and bind the result to the results data frame
  result <- get_model_metrics(predictor, df2)
  result_df <- rbind(result_df, result)
}

print(result_df)

(boot.comp.vi2 <- c(vi2[vi2 != 'comp1'], result_df$comp.name[result_df$p_value < 0.05]))   # specific list of covariates (vi) + boot.comp for later supplementary analysis

#~~~~~~~~~~~~~~#
#~~NULL MODEL~~#     #7# HMO components already generated. Now test first null model performance.
#~~~~~~~~~~~~~~#

summary(null_model <- lm(as.formula(paste('bsid_cog_ss', paste(setdiff(comp.vi2, c('comp_6m1', 'comp_6m2')), collapse = '+'), sep = '~')), df2))  
olsrr::ols_test_normality(null_model)   
olsrr::ols_test_breusch_pagan(null_model)   
any(olsrr::ols_coll_diag(null_model)$vif_t$VIF > 2)  

#12# function to look, flag and remove for outliers if normality or heteroskedasticity assumptions are violated
if(olsrr::ols_test_outlier(null_model)$bonferroni_p_val < 0.05 &
   olsrr::ols_test_normality(null_model)[[1]]$p.value < 0.05 | 
   olsrr::ols_test_breusch_pagan(null_model)[[4]] < 0.05) {
  outliers <- as.numeric(rownames(olsrr::ols_test_outlier(null_model)))  
  df.no.outlier <- df2[-outliers, ]         
  paste0('outlier present: ', outliers)
} else {
  df.no.outlier <- df2                      
  paste0('no outliers')
}

#13# re-run null model in new df (with no outlier) and check assumptions again
summary(null_model.new <- lm(as.formula(paste('bsid_cog_ss', paste(setdiff(comp.vi2, c('comp_6m1', 'comp_6m2')), collapse = '+'), sep = '~')), df.no.outlier))
sqrt(mean(null_model.new$residuals^2))
olsrr::ols_test_normality(null_model.new)
olsrr::ols_test_breusch_pagan(null_model.new)

sqrt(mean(null_model.new$residuals^2))
confint(null_model.new)            

# # # # # #
# MODEL 1 #     #11# fit/compare final model (only difference with null model is the addition of HMO comps)
# # # # # #

anova(null_model.new, model1 <- lm(as.formula(paste('bsid_cog_ss', paste(setdiff(comp.vi2, 'comp_6m2'), collapse = '+'), sep = '~')), df.no.outlier)) 
cohens_f(null_model.new, model2 = model1, alternative = 'two.sided')  
(p.val4 <- as.numeric(anova(null_model.new, model1)[[6]][2]))   

olsrr::ols_test_normality(model1)           
olsrr::ols_test_breusch_pagan(model1)       
any(olsrr::ols_coll_diag(model1)$vif_t$VIF > 2)   

summary(model1)       
sqrt(mean(model1$residuals^2))
summary(model1)[[9]]                   
confint(model1, 'comp_6m1')                

# # # # # #
# MODEL 2 #
# # # # # #

anova(null_model.new, model2 <- lm(as.formula(paste('bsid_cog_ss', paste(setdiff(comp.vi2, 'comp_6m1'), collapse = '+'), sep = '~')), df.no.outlier))  
cohens_f(null_model.new, model2 = model2, alternative = 'two.sided')
(p.val5 <- as.numeric(anova(null_model.new, model2)[[6]][2]))

olsrr::ols_test_normality(model2)            
olsrr::ols_test_breusch_pagan(model2)        
any(olsrr::ols_coll_diag(model2)$vif_t$VIF > 2)   

summary(model2)
sqrt(mean(model2$residuals^2))
summary(model2)[[9]]
confint(model2, 'comp_6m2')

# Select the best performing model
models <- list(model1, model2)

results_model <- as.data.frame(matrix(nrow = length(models), ncol = 2))

colnames(results_model) <- c('comp', 'pval')

for (i in seq_along(models)) {
  results_model[i,1] <- setdiff(names(models[[i]]$model), names(null_model$model))
  results_model[i,2] <- anova(null_model.new, models[[i]])$Pr[2]
}

(best_comp <- results_model$comp[which.min(results_model$pval)])

# # # # # # # #
#  BOOT MODEL #
# # # # # # # #

summary(null_model.new <- lm(as.formula(paste('bsid_cog_ss', paste(setdiff(boot.comp.vi2, 'boot.comp_6m1'), collapse = '+'), sep = '~')), df2))  
anova(null_model.new, boot.final_model <- lm(as.formula(paste('bsid_cog_ss', paste(boot.comp.vi2, collapse = '+'), sep = '~')), df.no.outlier))  
(p.val6 <- as.numeric(anova(null_model.new, boot.final_model)[[6]][2]))
cohens_f(null_model.new, model2 = boot.final_model, alternative = 'two.sided')

olsrr::ols_test_normality(boot.final_model)
olsrr::ols_test_breusch_pagan(boot.final_model)
any(olsrr::ols_coll_diag(boot.final_model)$vif_t$VIF > 2)   

summary(boot.final_model)
sqrt(mean(boot.final_model$residuals^2))
summary(boot.final_model)[[9]]
confint(boot.final_model, 'boot.comp1')

#14 and FINAL# Adjust all p-values 
(raw.pvals <- c(p.val1, p.val2, p.val3, p.val4, p.val5, p.val6))
(p.adjust(raw.pvals, method = 'holm'))


