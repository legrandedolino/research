#############################
# Analysis of PVs for each model
#############################

library(mice)
library(dplyr)
library(stringr)

# setting path
path <- "/Users/ledolino/Documents/LSA/mirt version/"

# loading data and model results
load(paste0(path, "data/scoreddata.RData"))

# just having a quick look
cor(bq_pvs$EMP.1.1, bq_pvs$SCH.1.1) # all 0.8X
cor(bq_pvs$EMP.1.1, bq_pvs$ELEVPC.1.1)
cor(bq_pvs$EMP.1.1, bq_pvs$thm)
cor(bq_pvs$ELEVPC.1.1, bq_pvs$thm) 

cor(bq_pvs$EMP.1.2, bq_pvs$SCH.1.2) # all 0.90X
cor(bq_pvs$EMP.1.2, bq_pvs$ELEVPC.1.2)
cor(bq_pvs$EMP.1.2, bq_pvs$thr)

cor(bq_pvs$EMP.1.3, bq_pvs$SCH.1.3) # all 0.9X
cor(bq_pvs$EMP.1.3, bq_pvs$ELEVPC.1.3)
cor(bq_pvs$EMP.1.3, bq_pvs$ths)

###### 
# 1. calculating difference between all pv x model estimates and true "theta" values on each domain
######
# create an empty data frame to store the results
differences <- data.frame(column_name = character(), average_difference = numeric(), stringsAsFactors = FALSE)

# loop through the columns in steps of 3, starting from column 38
for (i in seq(192, 579, by = 3)) {
  # for each set of 3 columns, calculate the average difference for thm, thr, and ths
  thm_diff <- mean(bq_pvs$thm - bq_pvs[, i])
  thr_diff <- mean(bq_pvs$thr - bq_pvs[, i + 1])
  ths_diff <- mean(bq_pvs$ths - bq_pvs[, i + 2])
  
  # store the results
  differences <- rbind(differences,
                   data.frame(column_name = colnames(bq_pvs)[i],
                              average_difference = thm_diff))
  differences <- rbind(differences,
                   data.frame(column_name = colnames(bq_pvs)[i + 1],
                              average_difference = thr_diff))
  differences <- rbind(differences,
                   data.frame(column_name = colnames(bq_pvs)[i + 2],
                              average_difference = ths_diff))
}

# format the average differences to four decimal places
differences$average_difference <- formatC(differences$average_difference, format = "f", digits = 4)

# print the results
print(differences)

# separating average differences into three separate objects
avg_diff_Dim1 <- differences[grep("\\.[1-3]\\.1$", differences$column_name), ]
avg_diff_Dim2 <- differences[grep("\\.[1-3]\\.2$", differences$column_name), ]
avg_diff_Dim3 <- differences[grep("\\.[1-3]\\.3$", differences$column_name), ]

# display the results
print("Average Differences for Dim1:")
print(avg_diff_Dim1)

print("Average Differences for Dim2:")
print(avg_diff_Dim2)

print("Average Differences for Dim3:")
print(avg_diff_Dim3)

###########
# 2. calculating correlation between all pv x model estimates and true "theta" values on each domain
###########

# create an empty data frame to store the results
correlation_results <- data.frame(column_name = character(), correlation = numeric(), stringsAsFactors = FALSE)

# loop through the columns in steps of 3, starting from column 38
for (i in seq(192, 579, by = 3)) {
  # for each set of 3 columns, calculate the correlation for thm, thr, and ths
  thm_cor <- cor(bq_pvs$thm, bq_pvs[, i])
  thr_cor <- cor(bq_pvs$thr, bq_pvs[, i + 1])
  ths_cor <- cor(bq_pvs$ths, bq_pvs[, i + 2])
  
  # store the results
  correlation_results <- rbind(correlation_results,
                               data.frame(column_name = colnames(bq_pvs)[i],
                                          correlation = thm_cor))
  correlation_results <- rbind(correlation_results,
                               data.frame(column_name = colnames(bq_pvs)[i + 1],
                                          correlation = thr_cor))
  correlation_results <- rbind(correlation_results,
                               data.frame(column_name = colnames(bq_pvs)[i + 2],
                                          correlation = ths_cor))
}

# format the correlations to four decimal places
correlation_results$correlation <- formatC(correlation_results$correlation, format = "f", digits = 4)

# print the results
print(correlation_results)

# separating correlations into three separate objects
correlation_Dim1 <- correlation_results[grep("\\.[1-3]\\.1$", correlation_results$column_name), ]
correlation_Dim2 <- correlation_results[grep("\\.[1-3]\\.2$", correlation_results$column_name), ]
correlation_Dim3 <- correlation_results[grep("\\.[1-3]\\.3$", correlation_results$column_name), ]

# display the results
print("Correlations for Dim1:")
print(correlation_Dim1)

print("Correlations for Dim2:")
print(correlation_Dim2)

print("Correlations for Dim3:")
print(correlation_Dim3)

################################
# calculating means and variance 
# of PVs by model:
################################

# initialize lists to store within-imputation means and variances
within_means <- list()
within_variances <- list()

# define a pattern to match the columns of interest
pattern <- "^(EMP|SCH|ONEPC|TWOPC|THREEPC|FOURPC|FIVEPC|SIXPC|SEVPC|EIGHTPC|NINEPC|TENPC|ELEVPC)\\.[1-9]0?\\.Dim[1-3]$"

# loop through each imputation and compute mean and variance for columns matching the pattern
for (i in 1:length(list_of_pvs)) {
  # selecting only the columns of interest based on the pattern
  df <- list_of_pvs[[i]]
  cols_of_interest <- grep(pattern, names(df), value = TRUE)
  df <- df[, cols_of_interest]
  
  # calculate means and variances for these columns
  within_means[[i]] <- sapply(df, mean, na.rm = TRUE)
  within_variances[[i]] <- sapply(df, var, na.rm = TRUE)
}

# convert lists to data frames
within_means_df <- do.call(rbind, within_means)
within_variances_df <- do.call(rbind, within_variances)

# calculate between-imputation mean and variance for each model and dimension
between_means <- colMeans(within_means_df, na.rm = TRUE)
between_variances <- apply(within_variances_df, 2, function(x) var(x, na.rm = TRUE))

# calculate overall mean according to Rubin's rules
final_mean <- between_means
print(final_mean)

# calculate total variance according to Rubin's rules
m <- length(list_of_pvs)  # number of imputations
within_mean_variances <- colMeans(within_variances_df, na.rm = TRUE)
final_variance <- within_mean_variances + (1 + 1/m) * between_variances

# clean up the model dimension names by removing the imputation index
clean_labels <- stringr::str_replace(names(final_mean), "\\d+", "")

# organizing results into a data frame for better readability
results <- data.frame(
  Model_Dimension = clean_labels,
  Mean = final_mean,
  Variance = final_variance,
  row.names = NULL  # Set row names to NULL
)

# adding generating values to the results
truevals <- apply(bq_pvs[,2:4], 2, function(x) c(Mean = mean(x), Variance = var(x)))
truevals <- as.data.frame(t(truevals))
truevals <- cbind(Model_Dimension = c("Dim1", "Dim2", "Dim3"), truevals)
results <- rbind(results, truevals)

# rounding results
results$Mean <- round(results$Mean, digits = 5)

print(results)

# separating results into three separate objects
dist_Dim1 <- results[grep("Dim1$", results$Model_Dimension), ]
dist_Dim2 <- results[grep("Dim2$", results$Model_Dimension), ]
dist_Dim3 <- results[grep("Dim3$", results$Model_Dimension), ]

print(dist_Dim1)
print(dist_Dim2)
print(dist_Dim3)


#####################
# saving results
#####################
save(list = c("avg_diff_Dim1", "avg_diff_Dim2", "avg_diff_Dim3", "correlation_Dim1",
              "correlation_Dim2", "correlation_Dim3", "dist_Dim1", "dist_Dim2", "dist_Dim3"
), file = paste0(path, "data/quality_metrics.RData"))

