#########################################
# by model regressions
#########################################

# setting path
path <- "/Users/ledolino/Documents/LSA/mirt version/"

# loading data and model results
load(paste0(path, "data/scoreddata.RData"))


#####
# q7 (binary)

# function to fit models for each imputation model x dimension across all dataframes
fit_models <- function(list_of_pvs, model, dim, predictor_formulas) {
  # store results for each predictor configuration
  results_per_predictor <- list()
  
  # iterate through each predictor formula
  for (predictor_formula in predictor_formulas) {
    models_list <- vector("list", length(list_of_pvs))
    
    # loop over each dataframe in the list (each dataframe is one imputation)
    for (i in seq_along(list_of_pvs)) {
      data <- list_of_pvs[[i]]
      outcome_name <- sprintf("%s.%d.Dim%d", model, i, dim)  # Construct the outcome variable name
      full_formula <- as.formula(paste(outcome_name, "~", predictor_formula))
      
      # fit the model if the outcome variable exists in the dataframe
      if (outcome_name %in% names(data)) {
        models_list[[i]] <- lm(full_formula, data = data)
      } else {
        models_list[[i]] <- NULL
      }
    }
    
    # remove NULL entries from the list
    models_list <- Filter(Negate(is.null), models_list)
    
    # apply Rubin's rules to pool results from models across imputations
    if (length(models_list) > 0) {
      pooled <- pool(models_list)
      results_per_predictor[[predictor_formula]] <- summary(pooled)
    } else {
      results_per_predictor[[predictor_formula]] <- NULL
    }
  }
  
  return(results_per_predictor)
}

# define the imputation models, dimensions, and predictor configurations
imputation_models <- c("EMP", "SCH", "ONEPC", "TWOPC", "THREEPC", "FOURPC", "FIVEPC", "SIXPC", "SEVPC", "EIGHTPC",
                       "NINEPC", "TENPC", "ELEVPC")
dimensions <- 1:3
predictors <- c("q7")

# storage for final results
q7_storage <- list()

# loop through each imputation model and dimension combination
for (model in imputation_models) {
  for (dim in dimensions) {
    key <- paste(model, "Dim", dim, sep = ".")
    # Fit models for the current combination across all dataframes
    q7_storage[[key]] <- fit_models(list_of_pvs, model, dim, predictors)
  }
}

# Summarize into a table the means and p-values of intercept and q7
q7_results <- data.frame(model = character(36), dim = character(36), par = character(36), 
                      est = numeric(36), std.error = numeric(36), p.value = numeric(36))

for(i in 1:39) {
  # Get model and dim names for the current index
  model_name <- sub("\\..*", "", names(q7_storage)[i])
  dim_name <- sub(".*\\.(Dim\\.\\d+)$", "\\1", names(q7_storage)[i])
  
  # Row indices for intercept and q7
  row_intercept <- 2 * i - 1
  row_q7 <- 2 * i
  
  # Fill in the row for "intercept"
  q7_results[row_intercept, "model"] <- model_name
  q7_results[row_intercept, "dim"] <- dim_name
  q7_results[row_intercept, "par"] <- "intercept"
  q7_results[row_intercept, "est"] <- q7_storage[[names(q7_storage)[i]]][["q7"]][["estimate"]][1]
  q7_results[row_intercept, "std.error"] <- q7_storage[[names(q7_storage)[i]]][["q7"]][["std.error"]][1]
  q7_results[row_intercept, "p.value"] <- q7_storage[[names(q7_storage)[i]]][["q7"]][["p.value"]][1]
  
  # Fill in the row for "q7"
  q7_results[row_q7, "model"] <- model_name
  q7_results[row_q7, "dim"] <- dim_name
  q7_results[row_q7, "par"] <- "q7"
  q7_results[row_q7, "est"] <- q7_storage[[names(q7_storage)[i]]][["q7"]][["estimate"]][2]
  q7_results[row_q7, "std.error"] <- q7_storage[[names(q7_storage)[i]]][["q7"]][["std.error"]][2]
  q7_results[row_q7, "p.value"] <- q7_storage[[names(q7_storage)[i]]][["q7"]][["p.value"]][2]
}

# append with regressions on true value
mod_thm <- lm(thm~q7,bqdum)
mod_thr <- lm(thr~q7,bqdum)
mod_ths <- lm(ths~q7,bqdum)

q7_results <- rbind(
  q7_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "intercept",
    est = coef(summary(mod_thm))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thm))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thm))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q7",
    est = coef(summary(mod_thm))["q72", "Estimate"],
    std.error = coef(summary(mod_thm))["q72", "Std. Error"],
    p.value = coef(summary(mod_thm))["q72", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for thr (Dim.2)
q7_results <- rbind(
  q7_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "intercept",
    est = coef(summary(mod_thr))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thr))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thr))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q7",
    est = coef(summary(mod_thr))["q72", "Estimate"],
    std.error = coef(summary(mod_thr))["q72", "Std. Error"],
    p.value = coef(summary(mod_thr))["q72", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for ths (Dim.3)
q7_results <- rbind(
  q7_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "intercept",
    est = coef(summary(mod_ths))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_ths))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_ths))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q7",
    est = coef(summary(mod_ths))["q72", "Estimate"],
    std.error = coef(summary(mod_ths))["q72", "Std. Error"],
    p.value = coef(summary(mod_ths))["q72", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

View(q7_results)
#####

#####
# q4 (continuous)

predictors <- c("q4")
q4_storage <- list()

for (model in imputation_models) {
  for (dim in dimensions) {
    key <- paste(model, "Dim", dim, sep = ".")
    q4_storage[[key]] <- fit_models(list_of_pvs, model, dim, predictors)
  }
}

q4_results <- data.frame(model = NULL, dim = NULL, par = NULL, 
                         est = NULL, std.error = NULL, p.value = NULL)

for(i in 1:39) {
  model_name <- sub("\\..*", "", names(q4_storage)[i])
  dim_name <- sub(".*\\.(Dim\\.\\d+)$", "\\1", names(q4_storage)[i])
  
  row_intercept <- 2 * i - 1
  row_q4 <- 2 * i
  
  q4_results[row_intercept, "model"] <- model_name
  q4_results[row_intercept, "dim"] <- dim_name
  q4_results[row_intercept, "par"] <- "intercept"
  q4_results[row_intercept, "est"] <- q4_storage[[names(q4_storage)[i]]][["q4"]][["estimate"]][1]
  q4_results[row_intercept, "std.error"] <- q4_storage[[names(q4_storage)[i]]][["q4"]][["std.error"]][1]
  q4_results[row_intercept, "p.value"] <- q4_storage[[names(q4_storage)[i]]][["q4"]][["p.value"]][1]
  
  q4_results[row_q4, "model"] <- model_name
  q4_results[row_q4, "dim"] <- dim_name
  q4_results[row_q4, "par"] <- "q4"
  q4_results[row_q4, "est"] <- q4_storage[[names(q4_storage)[i]]][["q4"]][["estimate"]][2]
  q4_results[row_q4, "std.error"] <- q4_storage[[names(q4_storage)[i]]][["q4"]][["std.error"]][2]
  q4_results[row_q4, "p.value"] <- q4_storage[[names(q4_storage)[i]]][["q4"]][["p.value"]][2]
}

# append with regressions on true value
mod_thm <- lm(thm~q4,bqdum)
mod_thr <- lm(thr~q4,bqdum)
mod_ths <- lm(ths~q4,bqdum)

q4_results <- rbind(
  q4_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "intercept",
    est = coef(summary(mod_thm))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thm))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thm))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q4",
    est = coef(summary(mod_thm))["q4", "Estimate"],
    std.error = coef(summary(mod_thm))["q4", "Std. Error"],
    p.value = coef(summary(mod_thm))["q4", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for thr (Dim.2)
q4_results <- rbind(
  q4_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "intercept",
    est = coef(summary(mod_thr))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thr))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thr))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q4",
    est = coef(summary(mod_thr))["q4", "Estimate"],
    std.error = coef(summary(mod_thr))["q4", "Std. Error"],
    p.value = coef(summary(mod_thr))["q4", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for ths (Dim.3)
q4_results <- rbind(
  q4_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "intercept",
    est = coef(summary(mod_ths))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_ths))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_ths))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q4",
    est = coef(summary(mod_ths))["q4", "Estimate"],
    std.error = coef(summary(mod_ths))["q4", "Std. Error"],
    p.value = coef(summary(mod_ths))["q4", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

View(q4_results)
#####

#####
# q5 (continuous)

predictors <- c("q5")
q5_storage <- list()

for (model in imputation_models) {
  for (dim in dimensions) {
    key <- paste(model, "Dim", dim, sep = ".")
    q5_storage[[key]] <- fit_models(list_of_pvs, model, dim, predictors)
  }
}

q5_results <- data.frame(model = NULL, dim = NULL, par = NULL, 
                         est = NULL, std.error = NULL, p.value = NULL)

for(i in 1:39) {
  model_name <- sub("\\..*", "", names(q5_storage)[i])
  dim_name <- sub(".*\\.(Dim\\.\\d+)$", "\\1", names(q5_storage)[i])
  
  row_intercept <- 2 * i - 1
  row_q5 <- 2 * i
  
  q5_results[row_intercept, "model"] <- model_name
  q5_results[row_intercept, "dim"] <- dim_name
  q5_results[row_intercept, "par"] <- "intercept"
  q5_results[row_intercept, "est"] <- q5_storage[[names(q5_storage)[i]]][["q5"]][["estimate"]][1]
  q5_results[row_intercept, "std.error"] <- q5_storage[[names(q5_storage)[i]]][["q5"]][["std.error"]][1]
  q5_results[row_intercept, "p.value"] <- q5_storage[[names(q5_storage)[i]]][["q5"]][["p.value"]][1]
  
  q5_results[row_q5, "model"] <- model_name
  q5_results[row_q5, "dim"] <- dim_name
  q5_results[row_q5, "par"] <- "q5"
  q5_results[row_q5, "est"] <- q5_storage[[names(q5_storage)[i]]][["q5"]][["estimate"]][2]
  q5_results[row_q5, "std.error"] <- q5_storage[[names(q5_storage)[i]]][["q5"]][["std.error"]][2]
  q5_results[row_q5, "p.value"] <- q5_storage[[names(q5_storage)[i]]][["q5"]][["p.value"]][2]
}

# append with regressions on true value
mod_thm <- lm(thm~q5,bqdum)
mod_thr <- lm(thr~q5,bqdum)
mod_ths <- lm(ths~q5,bqdum)

q5_results <- rbind(
  q5_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "intercept",
    est = coef(summary(mod_thm))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thm))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thm))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q5",
    est = coef(summary(mod_thm))["q5", "Estimate"],
    std.error = coef(summary(mod_thm))["q5", "Std. Error"],
    p.value = coef(summary(mod_thm))["q5", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for thr (Dim.2)
q5_results <- rbind(
  q5_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "intercept",
    est = coef(summary(mod_thr))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thr))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thr))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q5",
    est = coef(summary(mod_thr))["q5", "Estimate"],
    std.error = coef(summary(mod_thr))["q5", "Std. Error"],
    p.value = coef(summary(mod_thr))["q5", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for ths (Dim.3)
q5_results <- rbind(
  q5_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "intercept",
    est = coef(summary(mod_ths))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_ths))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_ths))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q5",
    est = coef(summary(mod_ths))["q5", "Estimate"],
    std.error = coef(summary(mod_ths))["q5", "Std. Error"],
    p.value = coef(summary(mod_ths))["q5", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

View(q5_results)
######

######
# q6 (continuous)

predictors <- c("q6")
q6_storage <- list()

for (model in imputation_models) {
  for (dim in dimensions) {
    key <- paste(model, "Dim", dim, sep = ".")
    q6_storage[[key]] <- fit_models(list_of_pvs, model, dim, predictors)
  }
}

q6_results <- data.frame(model = NULL, dim = NULL, par = NULL, 
                         est = NULL, std.error = NULL, p.value = NULL)

for(i in 1:39) {
  model_name <- sub("\\..*", "", names(q6_storage)[i])
  dim_name <- sub(".*\\.(Dim\\.\\d+)$", "\\1", names(q6_storage)[i])
  
  row_intercept <- 2 * i - 1
  row_q6 <- 2 * i
  
  q6_results[row_intercept, "model"] <- model_name
  q6_results[row_intercept, "dim"] <- dim_name
  q6_results[row_intercept, "par"] <- "intercept"
  q6_results[row_intercept, "est"] <- q6_storage[[names(q6_storage)[i]]][["q6"]][["estimate"]][1]
  q6_results[row_intercept, "std.error"] <- q6_storage[[names(q6_storage)[i]]][["q6"]][["std.error"]][1]
  q6_results[row_intercept, "p.value"] <- q6_storage[[names(q6_storage)[i]]][["q6"]][["p.value"]][1]
  
  q6_results[row_q6, "model"] <- model_name
  q6_results[row_q6, "dim"] <- dim_name
  q6_results[row_q6, "par"] <- "q6"
  q6_results[row_q6, "est"] <- q6_storage[[names(q6_storage)[i]]][["q6"]][["estimate"]][2]
  q6_results[row_q6, "std.error"] <- q6_storage[[names(q6_storage)[i]]][["q6"]][["std.error"]][2]
  q6_results[row_q6, "p.value"] <- q6_storage[[names(q6_storage)[i]]][["q6"]][["p.value"]][2]
}

# append with regressions on true value
mod_thm <- lm(thm~q6,bqdum)
mod_thr <- lm(thr~q6,bqdum)
mod_ths <- lm(ths~q6,bqdum)

q6_results <- rbind(
  q6_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "intercept",
    est = coef(summary(mod_thm))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thm))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thm))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q6",
    est = coef(summary(mod_thm))["q6", "Estimate"],
    std.error = coef(summary(mod_thm))["q6", "Std. Error"],
    p.value = coef(summary(mod_thm))["q6", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for thr (Dim.2)
q6_results <- rbind(
  q6_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "intercept",
    est = coef(summary(mod_thr))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thr))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thr))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q6",
    est = coef(summary(mod_thr))["q6", "Estimate"],
    std.error = coef(summary(mod_thr))["q6", "Std. Error"],
    p.value = coef(summary(mod_thr))["q6", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for ths (Dim.3)
q6_results <- rbind(
  q6_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "intercept",
    est = coef(summary(mod_ths))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_ths))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_ths))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q6",
    est = coef(summary(mod_ths))["q6", "Estimate"],
    std.error = coef(summary(mod_ths))["q6", "Std. Error"],
    p.value = coef(summary(mod_ths))["q6", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

View(q6_results)
#####

#####
# q8 (4 categories)

predictors <- c("q8")
q8_storage <- list()

for (model in imputation_models) {
  for (dim in dimensions) {
    key <- paste(model, "Dim", dim, sep = ".")
    q8_storage[[key]] <- fit_models(list_of_pvs, model, dim, predictors)
  }
}

q8_results <- data.frame(model = NULL, dim = NULL, par = NULL, 
                         est = NULL, std.error = NULL, p.value = NULL)

for(i in 1:39) {
  model_name <- sub("\\..*", "", names(q8_storage)[i])
  dim_name <- sub(".*\\.(Dim\\.\\d+)$", "\\1", names(q8_storage)[i])
  
  row_intercept <- 4 * i - 3
  row_q8_2 <- 4 * i - 2
  row_q8_3 <- 4 * i - 1
  row_q8_4 <- 4 * i
  
  q8_results[row_intercept, "model"] <- model_name
  q8_results[row_intercept, "dim"] <- dim_name
  q8_results[row_intercept, "par"] <- "intercept"
  q8_results[row_intercept, "est"] <- q8_storage[[names(q8_storage)[i]]][["q8"]][["estimate"]][1]
  q8_results[row_intercept, "std.error"] <- q8_storage[[names(q8_storage)[i]]][["q8"]][["std.error"]][1]
  q8_results[row_intercept, "p.value"] <- q8_storage[[names(q8_storage)[i]]][["q8"]][["p.value"]][1]
  
  q8_results[row_q8_2, "model"] <- model_name
  q8_results[row_q8_2, "dim"] <- dim_name
  q8_results[row_q8_2, "par"] <- "q8_2"
  q8_results[row_q8_2, "est"] <- q8_storage[[names(q8_storage)[i]]][["q8"]][["estimate"]][2]
  q8_results[row_q8_2, "std.error"] <- q8_storage[[names(q8_storage)[i]]][["q8"]][["std.error"]][2]
  q8_results[row_q8_2, "p.value"] <- q8_storage[[names(q8_storage)[i]]][["q8"]][["p.value"]][2]
  
  q8_results[row_q8_3, "model"] <- model_name
  q8_results[row_q8_3, "dim"] <- dim_name
  q8_results[row_q8_3, "par"] <- "q8_3"
  q8_results[row_q8_3, "est"] <- q8_storage[[names(q8_storage)[i]]][["q8"]][["estimate"]][3]
  q8_results[row_q8_3, "std.error"] <- q8_storage[[names(q8_storage)[i]]][["q8"]][["std.error"]][3]
  q8_results[row_q8_3, "p.value"] <- q8_storage[[names(q8_storage)[i]]][["q8"]][["p.value"]][3]
  
  q8_results[row_q8_4, "model"] <- model_name
  q8_results[row_q8_4, "dim"] <- dim_name
  q8_results[row_q8_4, "par"] <- "q8_4"
  q8_results[row_q8_4, "est"] <- q8_storage[[names(q8_storage)[i]]][["q8"]][["estimate"]][4]
  q8_results[row_q8_4, "std.error"] <- q8_storage[[names(q8_storage)[i]]][["q8"]][["std.error"]][4]
  q8_results[row_q8_4, "p.value"] <- q8_storage[[names(q8_storage)[i]]][["q8"]][["p.value"]][4]
}

# append with regressions on true value
mod_thm <- lm(thm~q8,bqdum)
mod_thr <- lm(thr~q8,bqdum)
mod_ths <- lm(ths~q8,bqdum)

q8_results <- rbind(
  q8_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "intercept",
    est = coef(summary(mod_thm))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thm))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thm))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q8_2",
    est = coef(summary(mod_thm))["q82", "Estimate"],
    std.error = coef(summary(mod_thm))["q82", "Std. Error"],
    p.value = coef(summary(mod_thm))["q82", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q8_3",
    est = coef(summary(mod_thm))["q83", "Estimate"],
    std.error = coef(summary(mod_thm))["q83", "Std. Error"],
    p.value = coef(summary(mod_thm))["q83", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q8_4",
    est = coef(summary(mod_thm))["q84", "Estimate"],
    std.error = coef(summary(mod_thm))["q84", "Std. Error"],
    p.value = coef(summary(mod_thm))["q84", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for thr (Dim.2)
q8_results <- rbind(
  q8_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "intercept",
    est = coef(summary(mod_thr))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thr))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thr))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q8_2",
    est = coef(summary(mod_thr))["q82", "Estimate"],
    std.error = coef(summary(mod_thr))["q82", "Std. Error"],
    p.value = coef(summary(mod_thr))["q82", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q8_3",
    est = coef(summary(mod_thr))["q83", "Estimate"],
    std.error = coef(summary(mod_thr))["q83", "Std. Error"],
    p.value = coef(summary(mod_thr))["q83", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q8_4",
    est = coef(summary(mod_thr))["q84", "Estimate"],
    std.error = coef(summary(mod_thr))["q84", "Std. Error"],
    p.value = coef(summary(mod_thr))["q84", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for ths (Dim.3)
q8_results <- rbind(
  q8_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "intercept",
    est = coef(summary(mod_ths))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_ths))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_ths))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q8_2",
    est = coef(summary(mod_ths))["q82", "Estimate"],
    std.error = coef(summary(mod_ths))["q82", "Std. Error"],
    p.value = coef(summary(mod_ths))["q82", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q8_3",
    est = coef(summary(mod_ths))["q83", "Estimate"],
    std.error = coef(summary(mod_ths))["q83", "Std. Error"],
    p.value = coef(summary(mod_ths))["q83", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q8_4",
    est = coef(summary(mod_ths))["q84", "Estimate"],
    std.error = coef(summary(mod_ths))["q84", "Std. Error"],
    p.value = coef(summary(mod_ths))["q84", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

View(q8_results)
#####

######
# q9 (4 categories)

predictors <- c("q9")
q9_storage <- list()

for (model in imputation_models) {
  for (dim in dimensions) {
    key <- paste(model, "Dim", dim, sep = ".")
    q9_storage[[key]] <- fit_models(list_of_pvs, model, dim, predictors)
  }
}

q9_results <- data.frame(model = NULL, dim = NULL, par = NULL, 
                         est = NULL, std.error = NULL, p.value = NULL)

for(i in 1:39) {
  model_name <- sub("\\..*", "", names(q9_storage)[i])
  dim_name <- sub(".*\\.(Dim\\.\\d+)$", "\\1", names(q9_storage)[i])
  
  row_intercept <- 4 * i - 3
  row_q9_2 <- 4 * i - 2
  row_q9_3 <- 4 * i - 1
  row_q9_4 <- 4 * i
  
  q9_results[row_intercept, "model"] <- model_name
  q9_results[row_intercept, "dim"] <- dim_name
  q9_results[row_intercept, "par"] <- "intercept"
  q9_results[row_intercept, "est"] <- q9_storage[[names(q9_storage)[i]]][["q9"]][["estimate"]][1]
  q9_results[row_intercept, "std.error"] <- q9_storage[[names(q9_storage)[i]]][["q9"]][["std.error"]][1]
  q9_results[row_intercept, "p.value"] <- q9_storage[[names(q9_storage)[i]]][["q9"]][["p.value"]][1]
  
  q9_results[row_q9_2, "model"] <- model_name
  q9_results[row_q9_2, "dim"] <- dim_name
  q9_results[row_q9_2, "par"] <- "q9_2"
  q9_results[row_q9_2, "est"] <- q9_storage[[names(q9_storage)[i]]][["q9"]][["estimate"]][2]
  q9_results[row_q9_2, "std.error"] <- q9_storage[[names(q9_storage)[i]]][["q9"]][["std.error"]][2]
  q9_results[row_q9_2, "p.value"] <- q9_storage[[names(q9_storage)[i]]][["q9"]][["p.value"]][2]
  
  q9_results[row_q9_3, "model"] <- model_name
  q9_results[row_q9_3, "dim"] <- dim_name
  q9_results[row_q9_3, "par"] <- "q9_3"
  q9_results[row_q9_3, "est"] <- q9_storage[[names(q9_storage)[i]]][["q9"]][["estimate"]][3]
  q9_results[row_q9_3, "std.error"] <- q9_storage[[names(q9_storage)[i]]][["q9"]][["std.error"]][3]
  q9_results[row_q9_3, "p.value"] <- q9_storage[[names(q9_storage)[i]]][["q9"]][["p.value"]][3]
  
  q9_results[row_q9_4, "model"] <- model_name
  q9_results[row_q9_4, "dim"] <- dim_name
  q9_results[row_q9_4, "par"] <- "q9_4"
  q9_results[row_q9_4, "est"] <- q9_storage[[names(q9_storage)[i]]][["q9"]][["estimate"]][4]
  q9_results[row_q9_4, "std.error"] <- q9_storage[[names(q9_storage)[i]]][["q9"]][["std.error"]][4]
  q9_results[row_q9_4, "p.value"] <- q9_storage[[names(q9_storage)[i]]][["q9"]][["p.value"]][4]
}

# append with regressions on true value
mod_thm <- lm(thm~q9,bqdum)
mod_thr <- lm(thr~q9,bqdum)
mod_ths <- lm(ths~q9,bqdum)

q9_results <- rbind(
  q9_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "intercept",
    est = coef(summary(mod_thm))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thm))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thm))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q9_2",
    est = coef(summary(mod_thm))["q92", "Estimate"],
    std.error = coef(summary(mod_thm))["q92", "Std. Error"],
    p.value = coef(summary(mod_thm))["q92", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q9_3",
    est = coef(summary(mod_thm))["q93", "Estimate"],
    std.error = coef(summary(mod_thm))["q93", "Std. Error"],
    p.value = coef(summary(mod_thm))["q93", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q9_4",
    est = coef(summary(mod_thm))["q94", "Estimate"],
    std.error = coef(summary(mod_thm))["q94", "Std. Error"],
    p.value = coef(summary(mod_thm))["q94", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for thr (Dim.2)
q9_results <- rbind(
  q9_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "intercept",
    est = coef(summary(mod_thr))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thr))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thr))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q9_2",
    est = coef(summary(mod_thr))["q92", "Estimate"],
    std.error = coef(summary(mod_thr))["q92", "Std. Error"],
    p.value = coef(summary(mod_thr))["q92", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q9_3",
    est = coef(summary(mod_thr))["q93", "Estimate"],
    std.error = coef(summary(mod_thr))["q93", "Std. Error"],
    p.value = coef(summary(mod_thr))["q93", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q9_4",
    est = coef(summary(mod_thr))["q94", "Estimate"],
    std.error = coef(summary(mod_thr))["q94", "Std. Error"],
    p.value = coef(summary(mod_thr))["q94", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for ths (Dim.3)
q9_results <- rbind(
  q9_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "intercept",
    est = coef(summary(mod_ths))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_ths))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_ths))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q9_2",
    est = coef(summary(mod_ths))["q92", "Estimate"],
    std.error = coef(summary(mod_ths))["q92", "Std. Error"],
    p.value = coef(summary(mod_ths))["q92", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q9_3",
    est = coef(summary(mod_ths))["q93", "Estimate"],
    std.error = coef(summary(mod_ths))["q93", "Std. Error"],
    p.value = coef(summary(mod_ths))["q93", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q9_4",
    est = coef(summary(mod_ths))["q94", "Estimate"],
    std.error = coef(summary(mod_ths))["q94", "Std. Error"],
    p.value = coef(summary(mod_ths))["q94", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

View(q9_results)
#####

#####
# q10 (4 categories)

predictors <- c("q10")
q10_storage <- list()

for (model in imputation_models) {
  for (dim in dimensions) {
    key <- paste(model, "Dim", dim, sep = ".")
    q10_storage[[key]] <- fit_models(list_of_pvs, model, dim, predictors)
  }
}

q10_results <- data.frame(model = NULL, dim = NULL, par = NULL, 
                         est = NULL, std.error = NULL, p.value = NULL)

for(i in 1:39) {
  model_name <- sub("\\..*", "", names(q10_storage)[i])
  dim_name <- sub(".*\\.(Dim\\.\\d+)$", "\\1", names(q10_storage)[i])
  
  row_intercept <- 4 * i - 3
  row_q10_2 <- 4 * i - 2
  row_q10_3 <- 4 * i - 1
  row_q10_4 <- 4 * i
  
  q10_results[row_intercept, "model"] <- model_name
  q10_results[row_intercept, "dim"] <- dim_name
  q10_results[row_intercept, "par"] <- "intercept"
  q10_results[row_intercept, "est"] <- q10_storage[[names(q10_storage)[i]]][["q10"]][["estimate"]][1]
  q10_results[row_intercept, "std.error"] <- q10_storage[[names(q10_storage)[i]]][["q10"]][["std.error"]][1]
  q10_results[row_intercept, "p.value"] <- q10_storage[[names(q10_storage)[i]]][["q10"]][["p.value"]][1]
  
  q10_results[row_q10_2, "model"] <- model_name
  q10_results[row_q10_2, "dim"] <- dim_name
  q10_results[row_q10_2, "par"] <- "q10_2"
  q10_results[row_q10_2, "est"] <- q10_storage[[names(q10_storage)[i]]][["q10"]][["estimate"]][2]
  q10_results[row_q10_2, "std.error"] <- q10_storage[[names(q10_storage)[i]]][["q10"]][["std.error"]][2]
  q10_results[row_q10_2, "p.value"] <- q10_storage[[names(q10_storage)[i]]][["q10"]][["p.value"]][2]
  
  q10_results[row_q10_3, "model"] <- model_name
  q10_results[row_q10_3, "dim"] <- dim_name
  q10_results[row_q10_3, "par"] <- "q10_3"
  q10_results[row_q10_3, "est"] <- q10_storage[[names(q10_storage)[i]]][["q10"]][["estimate"]][3]
  q10_results[row_q10_3, "std.error"] <- q10_storage[[names(q10_storage)[i]]][["q10"]][["std.error"]][3]
  q10_results[row_q10_3, "p.value"] <- q10_storage[[names(q10_storage)[i]]][["q10"]][["p.value"]][3]
  
  q10_results[row_q10_4, "model"] <- model_name
  q10_results[row_q10_4, "dim"] <- dim_name
  q10_results[row_q10_4, "par"] <- "q10_4"
  q10_results[row_q10_4, "est"] <- q10_storage[[names(q10_storage)[i]]][["q10"]][["estimate"]][4]
  q10_results[row_q10_4, "std.error"] <- q10_storage[[names(q10_storage)[i]]][["q10"]][["std.error"]][4]
  q10_results[row_q10_4, "p.value"] <- q10_storage[[names(q10_storage)[i]]][["q10"]][["p.value"]][4]
}

# append with regressions on true value
mod_thm <- lm(thm~q10,bqdum)
mod_thr <- lm(thr~q10,bqdum)
mod_ths <- lm(ths~q10,bqdum)

q10_results <- rbind(
  q10_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "intercept",
    est = coef(summary(mod_thm))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thm))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thm))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q10_2",
    est = coef(summary(mod_thm))["q102", "Estimate"],
    std.error = coef(summary(mod_thm))["q102", "Std. Error"],
    p.value = coef(summary(mod_thm))["q102", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q10_3",
    est = coef(summary(mod_thm))["q103", "Estimate"],
    std.error = coef(summary(mod_thm))["q103", "Std. Error"],
    p.value = coef(summary(mod_thm))["q103", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q10_4",
    est = coef(summary(mod_thm))["q104", "Estimate"],
    std.error = coef(summary(mod_thm))["q104", "Std. Error"],
    p.value = coef(summary(mod_thm))["q104", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for thr (Dim.2)
q10_results <- rbind(
  q10_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "intercept",
    est = coef(summary(mod_thr))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thr))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thr))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q10_2",
    est = coef(summary(mod_thr))["q102", "Estimate"],
    std.error = coef(summary(mod_thr))["q102", "Std. Error"],
    p.value = coef(summary(mod_thr))["q102", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q10_3",
    est = coef(summary(mod_thr))["q103", "Estimate"],
    std.error = coef(summary(mod_thr))["q103", "Std. Error"],
    p.value = coef(summary(mod_thr))["q103", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q10_4",
    est = coef(summary(mod_thr))["q104", "Estimate"],
    std.error = coef(summary(mod_thr))["q104", "Std. Error"],
    p.value = coef(summary(mod_thr))["q104", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for ths (Dim.3)
q10_results <- rbind(
  q10_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "intercept",
    est = coef(summary(mod_ths))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_ths))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_ths))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q10_2",
    est = coef(summary(mod_ths))["q102", "Estimate"],
    std.error = coef(summary(mod_ths))["q102", "Std. Error"],
    p.value = coef(summary(mod_ths))["q102", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q10_3",
    est = coef(summary(mod_ths))["q103", "Estimate"],
    std.error = coef(summary(mod_ths))["q103", "Std. Error"],
    p.value = coef(summary(mod_ths))["q103", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q10_4",
    est = coef(summary(mod_ths))["q104", "Estimate"],
    std.error = coef(summary(mod_ths))["q104", "Std. Error"],
    p.value = coef(summary(mod_ths))["q104", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

View(q10_results)
#####

#####
# q11 (4 categories)

predictors <- c("q11")
q11_storage <- list()

for (model in imputation_models) {
  for (dim in dimensions) {
    key <- paste(model, "Dim", dim, sep = ".")
    q11_storage[[key]] <- fit_models(list_of_pvs, model, dim, predictors)
  }
}

q11_results <- data.frame(model = NULL, dim = NULL, par = NULL, 
                         est = NULL, std.error = NULL, p.value = NULL)

for(i in 1:39) {
  model_name <- sub("\\..*", "", names(q11_storage)[i])
  dim_name <- sub(".*\\.(Dim\\.\\d+)$", "\\1", names(q11_storage)[i])
  
  row_intercept <- 4 * i - 3
  row_q11_2 <- 4 * i - 2
  row_q11_3 <- 4 * i - 1
  row_q11_4 <- 4 * i
  
  q11_results[row_intercept, "model"] <- model_name
  q11_results[row_intercept, "dim"] <- dim_name
  q11_results[row_intercept, "par"] <- "intercept"
  q11_results[row_intercept, "est"] <- q11_storage[[names(q11_storage)[i]]][["q11"]][["estimate"]][1]
  q11_results[row_intercept, "std.error"] <- q11_storage[[names(q11_storage)[i]]][["q11"]][["std.error"]][1]
  q11_results[row_intercept, "p.value"] <- q11_storage[[names(q11_storage)[i]]][["q11"]][["p.value"]][1]
  
  q11_results[row_q11_2, "model"] <- model_name
  q11_results[row_q11_2, "dim"] <- dim_name
  q11_results[row_q11_2, "par"] <- "q11_2"
  q11_results[row_q11_2, "est"] <- q11_storage[[names(q11_storage)[i]]][["q11"]][["estimate"]][2]
  q11_results[row_q11_2, "std.error"] <- q11_storage[[names(q11_storage)[i]]][["q11"]][["std.error"]][2]
  q11_results[row_q11_2, "p.value"] <- q11_storage[[names(q11_storage)[i]]][["q11"]][["p.value"]][2]
  
  q11_results[row_q11_3, "model"] <- model_name
  q11_results[row_q11_3, "dim"] <- dim_name
  q11_results[row_q11_3, "par"] <- "q11_3"
  q11_results[row_q11_3, "est"] <- q11_storage[[names(q11_storage)[i]]][["q11"]][["estimate"]][3]
  q11_results[row_q11_3, "std.error"] <- q11_storage[[names(q11_storage)[i]]][["q11"]][["std.error"]][3]
  q11_results[row_q11_3, "p.value"] <- q11_storage[[names(q11_storage)[i]]][["q11"]][["p.value"]][3]
  
  q11_results[row_q11_4, "model"] <- model_name
  q11_results[row_q11_4, "dim"] <- dim_name
  q11_results[row_q11_4, "par"] <- "q11_4"
  q11_results[row_q11_4, "est"] <- q11_storage[[names(q11_storage)[i]]][["q11"]][["estimate"]][4]
  q11_results[row_q11_4, "std.error"] <- q11_storage[[names(q11_storage)[i]]][["q11"]][["std.error"]][4]
  q11_results[row_q11_4, "p.value"] <- q11_storage[[names(q11_storage)[i]]][["q11"]][["p.value"]][4]
}

# append with regressions on true value
mod_thm <- lm(thm~q11,bqdum)
mod_thr <- lm(thr~q11,bqdum)
mod_ths <- lm(ths~q11,bqdum)

q11_results <- rbind(
  q11_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "intercept",
    est = coef(summary(mod_thm))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thm))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thm))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q11_2",
    est = coef(summary(mod_thm))["q112", "Estimate"],
    std.error = coef(summary(mod_thm))["q112", "Std. Error"],
    p.value = coef(summary(mod_thm))["q112", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q11_3",
    est = coef(summary(mod_thm))["q113", "Estimate"],
    std.error = coef(summary(mod_thm))["q113", "Std. Error"],
    p.value = coef(summary(mod_thm))["q113", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.1",
    par = "q11_4",
    est = coef(summary(mod_thm))["q114", "Estimate"],
    std.error = coef(summary(mod_thm))["q114", "Std. Error"],
    p.value = coef(summary(mod_thm))["q114", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for thr (Dim.2)
q11_results <- rbind(
  q11_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "intercept",
    est = coef(summary(mod_thr))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_thr))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_thr))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q11_2",
    est = coef(summary(mod_thr))["q112", "Estimate"],
    std.error = coef(summary(mod_thr))["q112", "Std. Error"],
    p.value = coef(summary(mod_thr))["q112", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q11_3",
    est = coef(summary(mod_thr))["q113", "Estimate"],
    std.error = coef(summary(mod_thr))["q113", "Std. Error"],
    p.value = coef(summary(mod_thr))["q113", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.2",
    par = "q11_4",
    est = coef(summary(mod_thr))["q114", "Estimate"],
    std.error = coef(summary(mod_thr))["q114", "Std. Error"],
    p.value = coef(summary(mod_thr))["q114", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

# Add rows for ths (Dim.3)
q11_results <- rbind(
  q11_results,
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "intercept",
    est = coef(summary(mod_ths))["(Intercept)", "Estimate"],
    std.error = coef(summary(mod_ths))["(Intercept)", "Std. Error"],
    p.value = coef(summary(mod_ths))["(Intercept)", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q11_2",
    est = coef(summary(mod_ths))["q112", "Estimate"],
    std.error = coef(summary(mod_ths))["q112", "Std. Error"],
    p.value = coef(summary(mod_ths))["q112", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q11_3",
    est = coef(summary(mod_ths))["q113", "Estimate"],
    std.error = coef(summary(mod_ths))["q113", "Std. Error"],
    p.value = coef(summary(mod_ths))["q113", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "truetheta",
    dim = "Dim.3",
    par = "q11_4",
    est = coef(summary(mod_ths))["q114", "Estimate"],
    std.error = coef(summary(mod_ths))["q114", "Std. Error"],
    p.value = coef(summary(mod_ths))["q114", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
)

View(q11_results)



#####################
# saving results
#####################
save(list = c("q4_results", "q5_results", "q6_results", "q7_results",
              "q8_results", "q9_results", "q10_results", "q11_results"
), file = paste0(path, "data/regressions.RData"))

setwd(path)
write.csv(q4_results, "q4_results.csv")
write.csv(q5_results, "q5_results.csv")
write.csv(q6_results, "q6_results.csv")
write.csv(q7_results, "q7_results.csv")
write.csv(q8_results, "q8_results.csv")
write.csv(q9_results, "q9_results.csv")
write.csv(q10_results, "q10_results.csv")
write.csv(q11_results, "q11_results.csv")


q5_results <- read.csv("q5_results.csv")
q11_results <- read.csv("q11_results.csv")
