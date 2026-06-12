library(lsasim)
library(MASS)
library(matrixcalc)
library(fastDummies)
library(performance)

###############################
# 1. Simulating multivariate normal thetas for 3 arbitrary domains
#
# Table 14.2a in PISA 2022 technical report shows that the international
# average correlation between reading, science, and math are:
# corr(M,R) = .80
# corr(M,S) = .85
# corr(R,S) = .79
#
# We'll use these as a starting point
#
# And simulating 8 background variables
# with a variety of distributions and 
# relationships with theta and each other
###############################

# generating questionnaires and thetas
# proportions for ordinal variables (continuous are "1")
# 3 theta, 1 continuous for school, 3 continuous, 1 binary, 4 ordinal
props <- list(1, 1, 1, 
              1, #q3 is school
              1, 1, 1, #q4 to q6 are continuous
              c(.50, 1), # q7 binary variable will be used as direct covariate. 
                         # correlations vary with other covariates and thetas.
              c(.20, .50, .80, 1), 
              c(.10, .40, .60, 1), 
              c(.20, .70, .90, 1), 
              c(.20, .50, .80, 1)
              )

# generating covariance matrix between theta and bq and among bq
# first 3 variables will be used as thetas for 3 arbitrary domains
yw_cov <- matrix(c(
  1.00, 0.80, 0.85, 0.40, 0.25, 0.30, 0.10, 0.35, 0.20, 0.45, 0.35, 0.15,
  0.80, 1.00, 0.79, 0.40, 0.20, 0.25, 0.15, 0.25, 0.20, 0.30, 0.40, 0.05,
  0.85, 0.79, 1.00, 0.40, 0.30, 0.35, 0.20, 0.30, 0.25, 0.15, 0.35, 0.10,
  0.40, 0.40, 0.40, 1.00, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20,
  0.25, 0.20, 0.30, 0.20, 1.00, 0.45, 0.30, 0.40, 0.35, 0.05, 0.20, 0.30,
  0.30, 0.25, 0.35, 0.20, 0.45, 1.00, 0.25, 0.20, 0.15, 0.40, 0.25, 0.25,
  0.10, 0.15, 0.20, 0.20, 0.30, 0.25, 1.00, 0.35, 0.30, 0.10, 0.15, 0.20,
  0.35, 0.25, 0.30, 0.20, 0.40, 0.20, 0.35, 1.00, 0.40, 0.35, 0.20, 0.15,
  0.20, 0.20, 0.25, 0.20, 0.35, 0.15, 0.30, 0.40, 1.00, 0.20, 0.10, 0.25,
  0.45, 0.30, 0.15, 0.20, 0.05, 0.40, 0.10, 0.35, 0.20, 1.00, 0.40, 0.30,
  0.35, 0.40, 0.35, 0.20, 0.20, 0.25, 0.15, 0.20, 0.10, 0.40, 1.00, 0.25,
  0.15, 0.05, 0.10, 0.20, 0.30, 0.25, 0.20, 0.15, 0.25, 0.30, 0.25, 1.00
), nrow = 12, byrow = TRUE)

# check if positive definite
if(is.positive.definite(yw_cov)) {
  cat("The generated matrix is positive definite.\n")
} else {
  cat("An error occurred. The matrix is not positive definite.\n")
}

# setting the seed
set.seed(178150)
n = 154*28 # 154 schools with 28 students each

bq <- questionnaire_gen(n, cov_matrix = yw_cov, cat_prop = props, theta= TRUE,
                         family = "gaussian")
colnames(bq)[2:4] <- c("thm", "thr", "ths")

# checking that theta distributions are close to intended
  # Define a function to calculate means and variances thetas
  calc_m_v <- function(data, variables) {
    # Use sapply to calculate both the mean and variance for each variable passed in the list
    stats <- sapply(variables, function(var) {
      mean_val <- mean(data[[var]], na.rm = TRUE)
      var_val <- var(data[[var]], na.rm = TRUE)
      return(c(mean = mean_val, variance = var_val))
    })
    # Return the statistics as a named matrix
    return(stats)
  }
  
  th_stats <- calc_m_v(bq, c("thm", "thr", "ths"))
  print(th_stats)

# checking that correlations are close to intended
overall_correlation <- cor(bq[,2:4]) 
print(overall_correlation)

#####################################
# 2. Simulating test for 3 domains
# - specify item numbers
# - specify item distributions
# - specify test design
# - assign booklets to examinees
#####################################

# number of  items (PISA technical report, chapter 2, p. 3)
Im <- 234 # math
Ir <- 197 # reading
Is <- 115 # science

# simulating items
# math
item_pool_m <- lsasim::item_gen(n_2pl = Im, b_bounds = c(-2.000, 2.000), 
                              a_bounds = c(0.305, 2.097)) # original c(0.105, 2.097)
summary(item_pool_m$b)
summary(item_pool_m$a)

# reading
item_pool_r <- lsasim::item_gen(n_2pl = Ir, b_bounds = c(-2.000, 2.000), 
                                a_bounds = c(0.377,	1.8277))
summary(item_pool_r$b)
summary(item_pool_r$a)

# science
item_pool_s <- lsasim::item_gen(n_2pl = Is, b_bounds = c(-2.000, 2.000), 
                                a_bounds = c(0.478, 2.188))
summary(item_pool_s$b)
summary(item_pool_s$a)

# simulating blocks (PISA Tech report, chapter 2, p. 7)
# math
K_m <- 18
blocks_m <- lsasim::block_design(n_blocks = K_m, item_parameters = item_pool_m, 
                               item_block_matrix = NULL)
print(blocks_m)

blockdescrip_m <- blocks_m$block_descriptives
print(blockdescrip_m)

# reading
K_r <- 10 #students got 35-42 reading items of 197 total items - ch. 2, p. 11)
blocks_r <- lsasim::block_design(n_blocks = K_r, item_parameters = item_pool_r, 
                                 item_block_matrix = NULL)
print(blocks_r)

blockdescrip_r <- blocks_r$block_descriptives
print(blockdescrip_r)

# science
K_s <- 6 # ch. 2, p. 12, figure 2.5
blocks_s <- lsasim::block_design(n_blocks = K_s, item_parameters = item_pool_s, 
                                 item_block_matrix = NULL)
print(blocks_s)

blockdescrip_s <- blocks_r$block_descriptives
print(blockdescrip_s)

# simulating booklets
# math
# book_design matrix with 18 blocks --> 18 forms
book_math <- matrix(c(
  1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 
), nrow=18, byrow=TRUE)
books_m <- lsasim::booklet_design(item_block_assignment=blocks_m$block_assignment, 
                                  book_design = book_math)
print(books_m)

# reading
# book_design matrix with 10 blocks --> 10 forms
book_reading <- matrix(c(
  1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
  1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=10, byrow=TRUE)

books_r <- lsasim::booklet_design(item_block_assignment=blocks_r$block_assignment, 
                                  book_design = book_reading)
print(books_r)

# science
# book_design matrix with 6 blocks --> 6 forms
book_science <- matrix(c(
  1, 1, 0, 0, 0, 0, 
  0, 1, 1, 0, 0, 0, 
  0, 0, 1, 1, 0, 0, 
  0, 0, 0, 1, 1, 0, 
  0, 0, 0, 0, 1, 1, 
  1, 0, 0, 0, 0, 1) , nrow=6, byrow=TRUE)
books_s <- lsasim::booklet_design(item_block_assignment=blocks_s$block_assignment,
                                  book_design = book_science)
print(books_s)

# assigning books to examinees
# math
bookassign_m <- lsasim::booklet_sample(n_subj=n, book_item_design=books_m, book_prob = NULL, 
                                     resample = FALSE, e = 0.1, iter = 20)

# reading
bookassign_r <- lsasim::booklet_sample(n_subj=n, book_item_design=books_r, book_prob = NULL, 
                                       resample = FALSE, e = 0.1, iter = 20)

# science
bookassign_s <- lsasim::booklet_sample(n_subj=n, book_item_design=books_s, book_prob = NULL, 
                                       resample = FALSE, e = 0.1, iter = 20)


# checking for equal assignment - proportion of sample assigned to each booklet
prop.table(table(bookassign_m[!duplicated(bookassign_m$subject), "book"]))
prop.table(table(bookassign_r[!duplicated(bookassign_r$subject), "book"]))
prop.table(table(bookassign_s[!duplicated(bookassign_s$subject), "book"]))

#####################################
# 3. Simulating responses to each test
# - simulate item responses by domain
# - combine all item responses, 
#   label by domain, 
#   and merge with questionnaires
#####################################

#math
sim_math <- response_gen(subject = bookassign_m$subject, item=bookassign_m$item, 
                                  theta = bq$thm,
                                  a_par = item_pool_m$a, b_par = item_pool_m$b)

#reading
sim_read <- response_gen(subject = bookassign_r$subject, item=bookassign_r$item, 
                         theta = bq$thr,
                         a_par = item_pool_r$a, b_par = item_pool_r$b)

#science
sim_sci <- response_gen(subject = bookassign_s$subject, item=bookassign_s$item, 
                         theta = bq$ths,
                         a_par = item_pool_s$a, b_par = item_pool_s$b)

# changing the name of the variables to "m"/"r"/"s" then item number
names(sim_math) <- gsub("^i", "m", names(sim_math))
names(sim_read) <- gsub("^i", "r", names(sim_read))
names(sim_sci)  <- gsub("^i", "s", names(sim_sci))

# converting q3 into pseudo-school variable with 154 levels
breaks <- quantile(bq$q3, probs = seq(0, 1, length.out = 155), na.rm = TRUE)
bq$school <- as.factor(cut(bq$q3, breaks = breaks, labels = FALSE, include.lowest = TRUE))
table(bq$school) # each pseudo-school has 28 students, like intended

# merging all of the test data together into a single dataframe and keeping ID from math
simdat <- cbind(sim_math[, c(235, 1:234)], sim_read[,-ncol(sim_read)], sim_sci[,-ncol(sim_sci)])

# double checking that there aren't multiple subject variables
duplicates <- duplicated(simdat | duplicated(simdat, fromLast = TRUE))

# get the names of the duplicate columns
duplicate_column_names <- names(simdat)[duplicates]

print(duplicate_column_names)

###################################
# 4. Processing background questionnaire data
# - dummy coding ordinal variables
# - creating principal components
# - combining with full BQ data
###################################

# converting ordinal variables to dummy coded numeric variables
dums <- fastDummies::dummy_cols(bq[,9:14], remove_first_dummy = TRUE)

# combining all bq data
bqdum <- cbind(bq, dums[,7:172])

# creating PCs
pcbq <- princomp(bqdum[,c(6:8,15:27)]) # schools excluded
summary(pcbq)  # 11 PCs account for 94.55% of variance

# selecting first 11 PCs to attach to BQ data
pcs <- pcbq$scores[,1:11]

bqdum <- cbind(bqdum, pcs)

# checking correlation of PCs
# get a list of variable names starting with "Comp."
comp_vars <- grep("^Comp\\.", names(bqdum), value = TRUE)

# check the list of variables
print(comp_vars)

# perform correlation analysis for these variables
correlation_matrix <- cor(bqdum[, comp_vars])

# print correlation matrix
print(round(correlation_matrix, 2))

###################################
# 5. Saving final data
###################################

# creating an object with all item parameters
# changing the name of the "item" variables to "m"/"r"/"s" then item number
item_pool_m$item <- sprintf("m%03d", item_pool_m$item)
item_pool_r$item <- sprintf("r%03d", item_pool_r$item)
item_pool_s$item <- sprintf("s%03d", item_pool_s$item)

genparms <- rbind(item_pool_m, item_pool_r, item_pool_s)

path <- "/Users/ledolino/Documents/LSA/mirt version/"

save(list = c("simdat", "bqdum", "pcbq", "genparms"), file = paste0(path, "data/alldata.RData"))


