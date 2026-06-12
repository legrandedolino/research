#################################
# Fitting latent regressions and
# drawing PVs. 
#################################

library(mirt)
library(ggplot2)
library(sfsmisc) # to use fewer nodes for item parameter estimation

# setting path
path <- "/Users/ledolino/Documents/LSA/mirt version/"

#######################
# Setting up analysis:
# 1. reading in fixed item parameters
# 2. assigning items to factors
######################

# loading data and item parameters
load(paste0(path, "data/alldata.RData"))
load(file=paste0(path, "data/mirtmodel_and_rescaledip.Rdata"))

#############################
# 1. Constructing starting 
#    values from calibrated model
#############################

# Fix model values by setting est = FALSE
calib_parms$est <- FALSE

# Allow varcovar to be estimated
calib_parms[calib_parms$name=="COV_11","est"] <- TRUE
calib_parms[calib_parms$name=="COV_22","est"] <- TRUE
calib_parms[calib_parms$name=="COV_33","est"] <- TRUE
calib_parms[calib_parms$name=="COV_21","est"] <- TRUE
calib_parms[calib_parms$name=="COV_31","est"] <- TRUE
calib_parms[calib_parms$name=="COV_32","est"] <- TRUE

#############################
# 2. Fit latent regression model in mirt 
#    with an empty latent regression model
#############################

# estimate 2PL model with fixed item parameters and empty latent regression
fit_mod_empty <- mirt(simdat[,2:547], full_mirt_mod,
                      itemtype = "2PL", pars = calib_parms,
                      technical = list(NCYCLES = 5000))

summary(mod2values(fit_mod_empty)$value - calib_parms$value) # Close enough

# draw plausible values
pv_empty <- fscores(fit_mod_empty, plausible.draws = 10)

#############################
# 3. Fit latent regression model in mirt
#    with direct regressor (binary variable)
#############################

# specify names of schools + principal components
latreg_names <- c(paste0("school_",2:154), paste0("Comp.", 1:11))

####
# mirt with school dummy
####

calib_parms_sch <- mirt(simdat[,2:547], full_mirt_mod,
                             itemtype = "2PL", pars = "values",
                             covdata = bqdum,
                             formula =  as.formula(paste("~", paste(latreg_names[1:153], collapse = " + "))),
                             technical = list(NCYCLES = 5000))
calib_parms_sch[1:3285,] <- calib_parms

fit_mod_sch <- mirt(simdat[,2:547], full_mirt_mod,
                    itemtype = "2PL", pars = calib_parms_sch,
                    covdata = bqdum,
                    formula =  as.formula(paste("~", paste(latreg_names[1:153], collapse = " + "))),
                    technical = list(NCYCLES = 5000))

# draw pvs
pv_sch <- fscores(fit_mod_sch, plausible.draws = 10)

#############################
# 4. Fit full latent regression model in mirt with
#   - school
#   - comp.1 - comp.11 (first 11 principal components)
#############################

calib_parms_full <- mirt(simdat[,2:547], full_mirt_mod,
                     itemtype = "2PL", pars = "values",
                     covdata = bqdum[, latreg_names],
                     formula =  as.formula(paste("~", paste(latreg_names[1:164], collapse = " + "))),
                     technical = list(NCYCLES = 5000))
calib_parms_full[1:3285,] <- calib_parms

fit_mod_full <- mirt(simdat[,2:547], full_mirt_mod,
                     itemtype = "2PL", pars = calib_parms_full,
                     covdata = bqdum[, latreg_names],
                     formula =  as.formula(paste("~", paste(latreg_names[1:164], collapse = " + "))),
                     technical = list(NCYCLES = 5000))

# draw pvs
pv_full_latreg <- fscores(fit_mod_full, plausible.draws = 10)


#############################
# 5. Fit partial latent regression models in mirt with
#   - direct, centered regressor
#   - up to 10 principal components
#############################

# estimate 2PL model with fixed item parameters and full latent regression
calib_parms_1PC <- mirt(simdat[,2:547], full_mirt_mod,
                         itemtype = "2PL", pars = "values",
                         covdata = bqdum[, latreg_names],
                         formula =  as.formula(paste("~", paste(latreg_names[1:154], collapse = " + "))),
                         technical = list(NCYCLES = 5000))
calib_parms_1PC[1:3285,] <- calib_parms

fit_mod_1PC <- mirt(simdat[,2:547], full_mirt_mod,
                    itemtype = "2PL", pars = calib_parms_1PC,
                    covdata = bqdum[, latreg_names],
                    formula =  as.formula(paste("~", paste(latreg_names[1:154], collapse = " + "))),
                    technical = list(NCYCLES = 5000))

# draw pvs
pv_1PC_latreg <- fscores(fit_mod_full, plausible.draws = 10)

# estimate 2PL model with fixed item parameters and full latent regression
calib_parms_2PC <- mirt(simdat[,2:547], full_mirt_mod,
                         itemtype = "2PL", pars = "values",
                         covdata = bqdum[, latreg_names],
                         formula =  as.formula(paste("~", paste(latreg_names[1:155], collapse = " + "))),
                         technical = list(NCYCLES = 5000))
calib_parms_2PC[1:3285,] <- calib_parms

fit_mod_2PC <- mirt(simdat[,2:547], full_mirt_mod,
                    itemtype = "2PL", pars = calib_parms_2PC,
                    covdata = bqdum[, latreg_names],
                    formula =  as.formula(paste("~", paste(latreg_names[1:155], collapse = " + "))),
                    technical = list(NCYCLES = 5000))

# draw pvs
pv_2PC_latreg <- fscores(fit_mod_full, plausible.draws = 10)

# estimate 2PL model with fixed item parameters and full latent regression
calib_parms_3PC <- mirt(simdat[,2:547], full_mirt_mod,
                        itemtype = "2PL", pars = "values",
                        covdata = bqdum[, latreg_names],
                        formula =  as.formula(paste("~", paste(latreg_names[1:156], collapse = " + "))),
                        technical = list(NCYCLES = 5000))
calib_parms_3PC[1:3285,] <- calib_parms

fit_mod_3PC <- mirt(simdat[,2:547], full_mirt_mod,
                    itemtype = "2PL", pars = calib_parms_3PC,
                    covdata = bqdum[, latreg_names],
                    formula =  as.formula(paste("~", paste(latreg_names[1:156], collapse = " + "))),
                    technical = list(NCYCLES = 5000))
# draw pvs
pv_3PC_latreg <- fscores(fit_mod_full, plausible.draws = 10)

# estimate 2PL model with fixed item parameters and full latent regression
calib_parms_4PC <- mirt(simdat[,2:547], full_mirt_mod,
                        itemtype = "2PL", pars = "values",
                        covdata = bqdum[, latreg_names],
                        formula =  as.formula(paste("~", paste(latreg_names[1:157], collapse = " + "))),
                        technical = list(NCYCLES = 5000))
calib_parms_4PC[1:3285,] <- calib_parms

fit_mod_4PC <- mirt(simdat[,2:547], full_mirt_mod,
                    itemtype = "2PL", pars = calib_parms_4PC,
                    covdata = bqdum[, latreg_names],
                    formula =  as.formula(paste("~", paste(latreg_names[1:157], collapse = " + "))),
                    technical = list(NCYCLES = 5000, internal_constraints = FALSE))

# draw pvs
pv_4PC_latreg <- fscores(fit_mod_full, plausible.draws = 10)

# estimate 2PL model with fixed item parameters and full latent regression
calib_parms_5PC <- mirt(simdat[,2:547], full_mirt_mod,
                        itemtype = "2PL", pars = "values",
                        covdata = bqdum[, latreg_names],
                        formula =  as.formula(paste("~", paste(latreg_names[1:158], collapse = " + "))),
                        technical = list(NCYCLES = 5000))
calib_parms_5PC[1:3285,] <- calib_parms

fit_mod_5PC <- mirt(simdat[,2:547], full_mirt_mod,
                    itemtype = "2PL", pars = calib_parms_5PC,
                    covdata = bqdum[, latreg_names],
                    formula =  as.formula(paste("~", paste(latreg_names[1:158], collapse = " + "))),
                    technical = list(NCYCLES = 5000))

# draw pvs
pv_5PC_latreg <- fscores(fit_mod_full, plausible.draws = 10)

# estimate 2PL model with fixed item parameters and full latent regression
calib_parms_6PC <- mirt(simdat[,2:547], full_mirt_mod,
                        itemtype = "2PL", pars = "values",
                        covdata = bqdum[, latreg_names],
                        formula =  as.formula(paste("~", paste(latreg_names[1:159], collapse = " + "))),
                        technical = list(NCYCLES = 5000))
calib_parms_6PC[1:3285,] <- calib_parms

fit_mod_6PC <- mirt(simdat[,2:547], full_mirt_mod,
                    itemtype = "2PL", pars = calib_parms_6PC,
                    covdata = bqdum[, latreg_names],
                    formula =  as.formula(paste("~", paste(latreg_names[1:159], collapse = " + "))),
                    technical = list(NCYCLES = 5000, internal_constraints = FALSE))

# draw pvs
pv_6PC_latreg <- fscores(fit_mod_full, plausible.draws = 10)

# estimate 2PL model with fixed item parameters and full latent regression
calib_parms_7PC <- mirt(simdat[,2:547], full_mirt_mod,
                        itemtype = "2PL", pars = "values",
                        covdata = bqdum[, latreg_names],
                        formula =  as.formula(paste("~", paste(latreg_names[1:160], collapse = " + "))),
                        technical = list(NCYCLES = 5000))
calib_parms_7PC[1:3285,] <- calib_parms

fit_mod_7PC <- mirt(simdat[,2:547], full_mirt_mod,
                    itemtype = "2PL", pars = calib_parms_7PC,
                    covdata = bqdum[, latreg_names],
                    formula =  as.formula(paste("~", paste(latreg_names[1:160], collapse = " + "))),
                    technical = list(NCYCLES = 5000))

# draw pvs
pv_7PC_latreg <- fscores(fit_mod_full, plausible.draws = 10)

# estimate 2PL model with fixed item parameters and full latent regression
calib_parms_8PC <- mirt(simdat[,2:547], full_mirt_mod,
                        itemtype = "2PL", pars = "values",
                        covdata = bqdum[, latreg_names],
                        formula =  as.formula(paste("~", paste(latreg_names[1:161], collapse = " + "))),
                        technical = list(NCYCLES = 5000))
calib_parms_8PC[1:3285,] <- calib_parms

fit_mod_8PC <- mirt(simdat[,2:547], full_mirt_mod,
                    itemtype = "2PL", pars = calib_parms_8PC,
                    covdata = bqdum[, latreg_names],
                    formula =  as.formula(paste("~", paste(latreg_names[1:161], collapse = " + "))),
                    technical = list(NCYCLES = 5000, internal_constraints = FALSE))

# draw pvs
pv_8PC_latreg <- fscores(fit_mod_full, plausible.draws = 10)

# estimate 2PL model with fixed item parameters and full latent regression
calib_parms_9PC <- mirt(simdat[,2:547], full_mirt_mod,
                        itemtype = "2PL", pars = "values",
                        covdata = bqdum[, latreg_names],
                        formula =  as.formula(paste("~", paste(latreg_names[1:162], collapse = " + "))),
                        technical = list(NCYCLES = 5000))
calib_parms_9PC[1:3285,] <- calib_parms

fit_mod_9PC <- mirt(simdat[,2:547], full_mirt_mod,
                    itemtype = "2PL", pars = calib_parms_9PC,
                    covdata = bqdum[, latreg_names],
                    formula =  as.formula(paste("~", paste(latreg_names[1:162], collapse = " + "))),
                    technical = list(NCYCLES = 5000))

# draw pvs
pv_9PC_latreg <- fscores(fit_mod_full, plausible.draws = 10)

# estimate 2PL model with fixed item parameters and full latent regression
calib_parms_10PC <- mirt(simdat[,2:547], full_mirt_mod,
                        itemtype = "2PL", pars = "values",
                        covdata = bqdum[, latreg_names],
                        formula =  as.formula(paste("~", paste(latreg_names[1:163], collapse = " + "))),
                        technical = list(NCYCLES = 5000))
calib_parms_10PC[1:3285,] <- calib_parms

fit_mod_10PC <- mirt(simdat[,2:547], full_mirt_mod,
                     itemtype = "2PL", pars = calib_parms_10PC,
                     covdata = bqdum[, latreg_names],
                     formula =  as.formula(paste("~", paste(latreg_names[1:163], collapse = " + "))),
                     technical = list(NCYCLES = 5000))

# draw pvs
pv_10PC_latreg <- fscores(fit_mod_full, plausible.draws = 10)

# rename PV columns to correspond to model
names(pv_empty) <- paste0("EMP.",1:10)
names(pv_sch) <- paste0("SCH.",1:10)
names(pv_full_latreg) <- paste0("ELEVPC.",1:10)
names(pv_1PC_latreg) <- paste0("ONEPC.",1:10)
names(pv_2PC_latreg) <- paste0("TWOPC.",1:10)
names(pv_3PC_latreg) <- paste0("THREEPC.",1:10)
names(pv_4PC_latreg) <- paste0("FOURPC.",1:10)
names(pv_5PC_latreg) <- paste0("FIVEPC.",1:10)
names(pv_6PC_latreg) <- paste0("SIXPC.",1:10)
names(pv_7PC_latreg) <- paste0("SEVPC.",1:10)
names(pv_8PC_latreg) <- paste0("EIGHTPC.",1:10)
names(pv_9PC_latreg) <- paste0("NINEPC.",1:10)
names(pv_10PC_latreg) <- paste0("TENPC.",1:10)

# combine all pvs with bqdum 
bq_pvs <- cbind(bqdum, pv_empty, pv_sch, pv_1PC_latreg, pv_2PC_latreg, pv_3PC_latreg, pv_4PC_latreg,
                pv_5PC_latreg, pv_6PC_latreg, pv_7PC_latreg, pv_8PC_latreg, pv_9PC_latreg, pv_10PC_latreg, pv_full_latreg)

###########################################
# need to convert back to list for each imputation and model 
# --> list with 10 dataframes that contain 1 imputation for 
#     each imputation model
###########################################

maintained_vars <- bq_pvs[, 1:191]

# Initialize list_of_pvs
list_of_pvs <- vector("list", 10)

# List of all the "pv" list names
pv_list_names <- c("pv_empty", "pv_sch", "pv_1PC_latreg", "pv_2PC_latreg", 
                   "pv_3PC_latreg", "pv_4PC_latreg", "pv_5PC_latreg", "pv_6PC_latreg", 
                   "pv_7PC_latreg", "pv_8PC_latreg", "pv_9PC_latreg", "pv_10PC_latreg", 
                   "pv_full_latreg")

# Loop over each of the 10 entries
for (i in 1:10) {
  # Start with maintained_vars as the base dataframe for each entry in list_of_pvs
  combined_df <- maintained_vars
  for (pv_name in pv_list_names) {   # Loop over each "pv" list name
    pv_list <- get(pv_name)     # Extract the ith entry from the current "pv" list
    pv_data <- pv_list[[i]]
    entry_name <- names(pv_list)[i]     # Retrieve the name of the current entry (e.g., "EMP.1", "EMP.2", etc.)
    colnames(pv_data) <- sprintf("%s.Dim%d", entry_name, 1:3)     # Rename columns to use entry name and Dim1, Dim2, Dim3 as suffixes
    combined_df <- cbind(combined_df, pv_data)     # Combine pv_data with the base dataframe (maintained_vars)
  }
  list_of_pvs[[i]] <- combined_df
}

# Display the first entry to verify
print(list_of_pvs[[1]])

# saving model results, PVs, and final data
save(list = c("bq_pvs","list_of_pvs",
              "pv_sch",
              "pv_empty", "pv_full_latreg", "pv_1PC_latreg", "pv_2PC_latreg",
              "pv_3PC_latreg","pv_4PC_latreg","pv_5PC_latreg","pv_6PC_latreg",
              "pv_7PC_latreg","pv_8PC_latreg","pv_9PC_latreg","pv_10PC_latreg",
              "fit_mod_sch",
              "fit_mod_empty", "fit_mod_full", "fit_mod_1PC", "fit_mod_2PC",
              "fit_mod_3PC","fit_mod_4PC","fit_mod_5PC","fit_mod_6PC",
              "fit_mod_7PC","fit_mod_8PC","fit_mod_9PC","fit_mod_10PC"
              ), file = paste0(path, "data/scoreddata.RData"))
