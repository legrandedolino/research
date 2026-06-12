library(mirt)
library(ggplot2)
# library(devtools) # need this to install sfmisc from github
# devtools::install_github("mmaechler/sfsmisc")
library(sfsmisc) # to use fewer nodes for item parameter estimation

#rm(list = ls())
###########################################
# Data analysis on simulated LSA data
# 1. Item calibration for each domain simultaneously
###########################################

# setting path
path <- "/Users/ledolino/Documents/LSA/mirt version/"

# loading data
load(paste0(path, "data/alldata.RData"))

###################################
# 1. Calibrating IRT parameters
###################################

# specifying mirt model
model_specification <- 'math = m001-m234 
                        read = r001-r197
                        sci  = s001-s115
                      # covariance
                        COV =  math*read*sci
                       '

full_mirt_mod <- mirt.model(model_specification, itemnames = simdat[,2:547])
fit_cal_full <- mirt(simdat[,2:547], full_mirt_mod, itemtype = '2PL', 
                     technical = list(NCYCLES = 5000), method = "MHRM") 
# Why MHRM: "Log-likelihood was decreasing near the ML solution. EM method may be unstable"

summary(fit_cal_full)

mirt::extract.mirt(fit_cal_full, "converged")

####################################
# 2. Checking generated against estimated
#    item parameters
####################################

# extracting item parameters to check recovery
mmirt_coef <- mirt::coef(fit_cal_full, simplify = TRUE, IRTpars = TRUE)
mirt_parms <- round(mmirt_coef$items, 2)

# confirming that correlation is close to original
cov2cor(mmirt_coef$cov) - cor(bqdum[,2:4])  # positive bias within ~0.02 units

# define a function to plot and compute correlations
plot_and_cor <- function(data_x, data_y, x_label, y_label) {
  plot(data_x, data_y, xlab = x_label, ylab = y_label)
  abline(a = 0, b = 1, col = "red", lty = 2) # Add 45-degree line
  return(cor(data_x, data_y))
}

# filter rows for math, reading and science (rows with 'm', 'r', 's')
math_mirt_rows <- grep("^m", rownames(mirt_parms))
reading_mirt_rows <- grep("^r", rownames(mirt_parms))
science_mirt_rows <- grep("^s", rownames(mirt_parms))

######################################
# Re-scaling item parameters and putting
# back into mirt metric
######################################

# extract (estimated) item parameters
math_est_a <- mirt_parms[math_mirt_rows, "a1"] #estimated descrimination/slope
math_est_b <- mirt_parms[math_mirt_rows, "b"]  #estimated easiness

read_est_a <- mirt_parms[reading_mirt_rows, "a2"] #estimated descrimination/slope
read_est_b <- mirt_parms[reading_mirt_rows, "b"]  #estimated easiness

sci_est_a <- mirt_parms[science_mirt_rows, "a3"] #estimated descrimination/slope
sci_est_b <- mirt_parms[science_mirt_rows, "b"]  #estimated easiness

# rescale estimated item parameters using mean/sigma method
# math
math_u <- mean(math_est_a)/mean(genparms[1:234,3])
math_v <- mean(genparms[1:234,2]) - math_u*mean(math_est_b)

rsmath_est_a <- math_est_a/math_u               #rescaled item discrimination
rsmath_est_b <- math_u*math_est_b + math_v      #rescaled item difficulty

# reading
read_u <- mean(read_est_a)/mean(genparms[235:432,3])
read_v <- mean(genparms[235:432,2]) - read_u*mean(read_est_b)

rsread_est_a <- read_est_a/read_u               #rescaled item discrimination
rsread_est_b <- read_u*read_est_b + read_v      #rescaled item difficulty

# science
sci_u <- mean(sci_est_a)/mean(genparms[433:546,3])
sci_v <- mean(genparms[433:546,2]) - sci_u*mean(sci_est_b)

rssci_est_a <- sci_est_a/sci_u               #rescaled item discrimination
rssci_est_b <- sci_u*sci_est_b + sci_v      #rescaled item difficulty

#####################################
# Converting rescaled b to mirt d
#####################################

# saving new object in mirt metric
calib_parms <- mod2values(fit_cal_full)

math_item_names <- grep("^m", rownames(mirt_parms), value=TRUE)
calib_parms[calib_parms$item %in% math_item_names & calib_parms$name == "a1","value"] <- rsmath_est_a
calib_parms[calib_parms$item %in% math_item_names & calib_parms$name == "d","value"] <- -rsmath_est_b*rsmath_est_a

read_item_names <- grep("^r", rownames(mirt_parms), value=TRUE)
calib_parms[calib_parms$item %in% read_item_names & calib_parms$name == "a2","value"] <- rsread_est_a
calib_parms[calib_parms$item %in% read_item_names & calib_parms$name == "d","value"] <- -rsread_est_b*rsread_est_a

sci_item_names <- grep("^s", rownames(mirt_parms), value=TRUE)
calib_parms[calib_parms$item %in% sci_item_names & calib_parms$name == "a3","value"] <- rssci_est_a
calib_parms[calib_parms$item %in% sci_item_names & calib_parms$name == "d","value"] <- -rssci_est_b*rssci_est_a

View(calib_parms)

###########################################
# Plotting rescaled parameters
###########################################

# math discrimination and easiness
par(mfrow=c(1,2))
plot_and_cor(genparms[math_mirt_rows, "a"], rsmath_est_a, 
             "Math Discrimination (itemparms)", "Math Discrimination (rescaled)")
plot_and_cor(-genparms[math_mirt_rows, "b"]*genparms[math_mirt_rows, "a"], -rsmath_est_b*rsmath_est_a, 
             "Math Easiness (itemparms)", "Math Easiness (rescaled)")

# reading discrimination and easiness
par(mfrow=c(1,2))
plot_and_cor(genparms[reading_mirt_rows, "a"], rsread_est_a, 
             "Reading Discrimination (genparms)", "Reading Discrimination (rescaled)")
plot_and_cor(-genparms[reading_mirt_rows, "b"]*genparms[reading_mirt_rows, "a"], -rsread_est_b*rsread_est_a, 
             "Reading Easiness (genparms)", "Reading Easiness (rescaled)")

# science discrimination and easiness
par(mfrow=c(1,2))
plot_and_cor(genparms[science_mirt_rows, "a"], rssci_est_a, 
             "Science Discrimination (genparms)", "Science Discrimination (rescaled)")
plot_and_cor(-genparms[science_mirt_rows, "b"]*genparms[science_mirt_rows, "a"], -rssci_est_b*rssci_est_a, 
             "Science Easiness (genparms)", "Science Easiness (rescaled)")


#####################################
# Saving original and rescaled item 
# parameters, as well as model spec 
# and item name patterns
#####################################

save(list = c("fit_cal_full", "calib_parms", "full_mirt_mod", "math_mirt_rows", 
              "science_mirt_rows", "reading_mirt_rows"), file=paste0(path, "data/mirtmodel_and_rescaledip.Rdata"))

