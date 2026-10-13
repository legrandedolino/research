###############################################################################
# 11a. Representativeness check run
#
# Purpose
#   Test whether the imputation approach restores the representativeness of
#   the sample where the truth is known, following Silverwood et al. (2024).
#   The age-14 (Wave 1) variables of the analysis models are known for most
#   of the cohort, including most people who later dropped out. Here they
#   are set to missing ("masked") for everyone who did not respond at Wave
#   9, and imputed again with the main imputation model. File 11b then
#   compares three estimates of their distributions:
#     truth           all known values (the whole target population);
#     complete cases  Wave-9 respondents' values only;
#     imputation      respondents' values plus the imputed values of the
#                     masked non-respondents.
#   If the imputation estimate is close to the truth while the complete-case
#   estimate is not, the model restores representativeness for these
#   variables. This shows recovery where the truth is known; it does not
#   prove that the age-32 data are missing at random. As Silverwood et al.
#   stress, the known values are masked for demonstration only; every other
#   analysis uses them.
#   Each imputation is one independent task with its own seed and output
#   file, so the run can be split across several computers (or several R
#   sessions on one computer) and resumed after an interruption.
#
# Survey weights and the imputation model (a deliberate departure)
#   Silverwood et al. (2024) weighted their imputation models by the initial
#   design weights. That is not done here. The imputation uses R's mice
#   package, which has no weighted imputation methods for categorical
#   variables, and most of the variables imputed here are categorical.
#   Instead:
#     - the log of the initial design weight and the sampling stratum are
#       included as predictors in every imputation equation (File 7). This
#       is a standard way of carrying survey design information into an
#       imputation model (Reiter, Raghunathan and Kinney 2006). As a
#       predictor, strata too small for the imputation models are pooled
#       into one group; the analysis uses the original strata;
#     - the imputation models themselves are unweighted. (The CLS user
#       guide's worked example is from a cohort without sampling weights,
#       so it doesn't settle this either way.) With the design variables
#       as predictors, weighting the imputation equations as well would
#       mainly matter if relationships between variables differ by
#       selection probability, i.e. weight interactions (below);
#     - survey weights, PSUs and strata are applied when each completed
#       dataset is analysed (File 11b).
#   Interactions between the weight and other variables (Seaman et al. 2012)
#   are not included. Silverwood et al. could not include them either,
#   because their models became unstable. They note that this should not
#   affect point estimates but may overstate standard errors, which makes
#   conclusions slightly conservative. The same applies here.
#
# What is run
#   File 7's specification, with one change: the nine age-14 analysis
#   variables (sex, ethnic group, family NS-SEC, parent qualification,
#   household income, and the four science and mathematics orientation
#   items) are masked together among Wave-9 non-respondents, and the where
#   matrix marks the masked values for imputation. Everything else is
#   unchanged: the auxiliaries (including the CLS predictors of
#   non-response), the design predictors, the routing, the methods, the
#   predictor matrix, the number of iterations and the number of
#   imputations (the same as the main run, so the CLS rule on the fraction
#   of missing information can be met; File 11b checks it).
#   Masking removes observed values, so File 7's check that every predictor
#   category has at least 20 people among the rows each bootstrap logistic
#   model is fitted to is repeated on the masked data before any task runs.
#
# Running on several computers or R sessions
#   As for File 8a: run Files 0-7 once, copy data_dir/07_imputation_spec.rds
#   to any other computer that will run tasks (only between computers
#   covered by the UKDS End User Licence and your institution's data
#   security rules), and give each computer or R session its own
#   tasks_to_run (Block 1). mice uses one processor core, so a computer with
#   several cores can run several R sessions at once. Copy every task file
#   into tasks_dir/11a on one computer and run File 11b there. A task whose
#   file already exists is skipped, and each task's seed is fixed (base
#   seed + task number).
#
# Reads
#   data_dir/07_imputation_spec.rds
#
# Writes
#   review_dir/11_category_check.csv  predictor categories below the minimum
#                                     on the masked data (names only)
#   data_dir/11_masked_spec.rds       the masked specification, with which
#                                     values were masked (participant level)
#   tasks_dir/11a/imputation_NNN.rds  one per task: the imputation (a mids
#       object with one imputation), plus its task number, seed, run time,
#       warnings, and the R and mice versions (participant level)
#
# Reference
#   Silverwood RJ, et al. (2024). Longitudinal and Life Course Studies 15(2):
#   227-250.
###############################################################################

suppressPackageStartupMessages(library(mice))


###############################################################################
# Block 1. Paths and settings
#
# Folder locations are set once, in File 0, and read here from
# project_paths.rds. Run this script with the working directory set to the
# folder that holds the scripts.
###############################################################################

paths      <- readRDS("project_paths.rds")
data_dir   <- paths$data_dir     # datasets passed between scripts
review_dir <- paths$review_dir   # disclosure-checked outputs
tasks_dir  <- paths$tasks_dir    # one file per task in split runs

# Which tasks this computer (or R session) runs: "all", or task numbers,
# e.g. 1:30. Every task in the table is needed before File 11b can combine.
tasks_to_run <- "all"

base_seed <- 11000               # task t uses seed base_seed + t
minimum_category_size <- 20      # as in File 7

# The age-14 variables of the analysis models (File 5), masked together.
masked_variables <- c("sex14", "ethnicity14", "nssec_family14",
                      "highest_parent_education14", "household_income14",
                      "science_liking14", "science_competence14",
                      "science_favourite14", "math_liking14")

task_folder <- file.path(tasks_dir, "11a")
dir.create(task_folder, recursive = TRUE, showWarnings = FALSE)


###############################################################################
# Block 2. Read File 7's specification
###############################################################################

imputation_spec <- readRDS(file.path(data_dir, "07_imputation_spec.rds"))


###############################################################################
# Block 3. The masked specification
#
# For each masked variable: among Wave-9 non-respondents, every observed
# value is set to missing and marked for imputation. Values that were
# missing anyway were already marked. Which cells were masked is kept, so
# File 11b can compare the imputed values with the truth.
###############################################################################

# Checks: every masked variable is in the model and imputed by a model of
# its own (not calculated from other variables).
if (!all(masked_variables %in% names(imputation_spec$data))) {
  stop("A masked variable is not in File 7's imputation model.")
}
masked_methods <- imputation_spec$method[masked_variables]
if (any(masked_methods == "") || any(startsWith(masked_methods, "~"))) {
  stop("Every masked variable needs its own imputation method in File 7's ",
       "specification (it has none, or is calculated from other variables).")
}

w9_nonrespondent <- !(imputation_spec$design_data$responded_w9 %in% TRUE)

masked_data  <- imputation_spec$data
masked_where <- imputation_spec$where
masked_cells <- matrix(FALSE, nrow = nrow(masked_data), ncol = length(masked_variables),
                       dimnames = list(NULL, masked_variables))
for (variable in masked_variables) {
  to_mask <- w9_nonrespondent & !is.na(masked_data[[variable]])
  masked_data[[variable]][to_mask] <- NA
  masked_where[to_mask, variable]  <- TRUE
  masked_cells[, variable]         <- to_mask
}


###############################################################################
# Block 4. Bootstrap logistic models: category check on the masked data
#
# As File 7: every category of every categorical (or 0/1) predictor must
# have at least 20 people among the rows each bootstrap logistic model is
# fitted to (the rows where the imputed variable is observed, now after
# masking). Each category below the minimum is listed by name (no counts),
# and the script stops before any task if one is found.
# What is counted:
#   - an unmasked predictor: its observed values, a safe lower bound (as in
#     File 7);
#   - a masked predictor: its true values. Its masked values are imputed
#     again at every iteration before the other models are fitted, so the
#     number in a category among a model's rows is the completed count,
#     whose expected size is the count of true values (File 7 checked those
#     before masking). Counting only the values left after masking would
#     understate what the model sees.
###############################################################################

# First the (model, predictor) pairs to check, then one loop over them.
check_pairs <- data.frame()
for (variable in names(imputation_spec$method)[imputation_spec$method == "logreg.boot"]) {
  predictors_v <- names(which(imputation_spec$predictor_matrix[variable, ] == 1))
  check_pairs <- rbind(check_pairs, data.frame(model_for = rep(variable, length(predictors_v)),
                                               predictor = predictors_v))
}
category_check <- data.frame(model_for = character(), predictor = character(),
                             category = character(), fewer_than = numeric())
for (j in seq_len(nrow(check_pairs))) {
  predictor   <- check_pairs$predictor[j]
  fitted_rows <- !is.na(masked_data[[check_pairs$model_for[j]]])
  values <- masked_data[[predictor]][fitted_rows]                              # observed values
  if (predictor %in% masked_variables) values <- imputation_spec$data[[predictor]][fitted_rows]   # true values
  is_categorical <- is.factor(values) || all(values %in% c(0, 1, NA))
  counts <- table(values)
  if (is_categorical && any(counts < minimum_category_size)) {
    category_check <- rbind(category_check, data.frame(
      model_for = check_pairs$model_for[j], predictor = predictor,
      category = names(counts)[counts < minimum_category_size],
      fewer_than = minimum_category_size
    ))
  }
}

# Disclosure check: names only. A listed category has fewer than 20 people
# among the model's rows; the exact count is not shown.
category_check$disclosure_checked <- rep(TRUE, nrow(category_check))
write.csv(category_check, file.path(review_dir, "11_category_check.csv"),
          row.names = FALSE)

if (nrow(category_check) > 0) {
  stop(nrow(category_check), " predictor categories have fewer than ",
       minimum_category_size, " people among the rows a bootstrap logistic ",
       "model is fitted to after masking; see 11_category_check.csv.")
}


###############################################################################
# Block 5. Save the masked specification
#
# Several R sessions may start at once and all save this (identical) file,
# so each writes under its own temporary name and then renames it.
###############################################################################

masked_spec <- list(
  data             = masked_data,
  design_data      = imputation_spec$design_data,
  variables        = imputation_spec$variables,
  method           = imputation_spec$method,
  predictor_matrix = imputation_spec$predictor_matrix,
  where            = masked_where,
  visit_sequence   = imputation_spec$visit_sequence,
  iterations       = imputation_spec$iterations,
  production_m     = imputation_spec$production_m,
  masked_variables = masked_variables,
  masked_cells     = masked_cells
)

spec_file      <- file.path(data_dir, "11_masked_spec.rds")
temporary_spec <- paste0(spec_file, ".", Sys.getpid(), ".partial")
saveRDS(masked_spec, temporary_spec)
invisible(file.rename(temporary_spec, spec_file))


###############################################################################
# Block 6. Task table
#
# One row per imputation: its task number, seed and output file.
###############################################################################

task_table <- data.frame(task = seq_len(masked_spec$production_m))
task_table$seed <- base_seed + task_table$task
task_table$file <- file.path(task_folder,
                             sprintf("imputation_%03d.rds", task_table$task))

if (identical(tasks_to_run, "all")) tasks_to_run <- task_table$task
if (!all(tasks_to_run %in% task_table$task)) {
  stop("tasks_to_run includes task numbers outside 1-", nrow(task_table), ".")
}


###############################################################################
# Block 7. Workspace check
#
# mice finds two of its internal objects, "state" and "loggedEvents", by
# name. If the R workspace holds an object with either name (e.g. from other
# code run in the same session), mice picks it up instead and fails with
# "$ operator is invalid for atomic vectors" the first time it logs an
# event. Stop now with a clear message instead.
###############################################################################

for (name in c("state", "loggedEvents")) {
  if (exists(name, envir = globalenv(), inherits = FALSE)) {
    stop("An object called '", name, "' is in the R workspace, and mice would ",
         "mistake it for its own. Remove it with rm(", name, ") or restart R, ",
         "then run this script again.")
  }
}


###############################################################################
# Block 8. Run the tasks
#
# Each task is independent: it imputes once, from its own seed, and saves
# its own file. The file is first written under a temporary name and then
# renamed, so an interruption during saving never leaves a half-written file
# that a later run would skip. Warnings are kept with the task, not printed.
# Progress lines show task numbers and times only.
###############################################################################

cat("\nRepresentativeness check run (", length(masked_variables),
    " age-14 variables masked among Wave-9 non-respondents): ",
    length(tasks_to_run), " of ", nrow(task_table), " tasks on this computer; ",
    masked_spec$iterations, " iterations each.\n", sep = "")

for (t in tasks_to_run) {
  task_file <- task_table$file[t]

  # Skip tasks that are already done.
  if (file.exists(task_file)) {
    cat("  Task", t, "already done; skipped.\n")
    next
  }

  cat("  Task", t, "started", format(Sys.time(), "%H:%M"), "\n")
  start_time <- Sys.time()

  task_warnings <- character()
  imputation <- withCallingHandlers(
    mice(masked_spec$data,
         m               = 1,
         maxit           = masked_spec$iterations,
         method          = masked_spec$method,
         predictorMatrix = masked_spec$predictor_matrix,
         where           = masked_spec$where,
         visitSequence   = masked_spec$visit_sequence,
         seed            = task_table$seed[t],
         printFlag       = FALSE),
    warning = function(w) {
      # Each time a warning occurs, keep its message and let the imputation
      # carry on ("<<-" adds it to task_warnings, outside this small function).
      task_warnings <<- c(task_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  task_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

  task_result <- list(
    task         = t,
    seed         = task_table$seed[t],
    imputation   = imputation,
    minutes      = task_minutes,
    warnings     = task_warnings,
    r_version    = as.character(getRversion()),
    mice_version = as.character(packageVersion("mice")),
    computer     = gsub("[^A-Za-z0-9_-]", "_", Sys.info()[["nodename"]])
  )

  temporary_file <- paste0(task_file, ".partial")
  saveRDS(task_result, temporary_file)
  file.rename(temporary_file, task_file)

  cat("  Task", t, "finished in", round(task_minutes, 1), "minutes",
      if (length(task_warnings) > 0) paste0("(", length(task_warnings), " warnings kept)") else "",
      "\n")
}


###############################################################################
# Block 9. Summary
#
# Task counts only. No data values.
###############################################################################

tasks_done <- sum(file.exists(task_table$file))
cat("\nFile 11a finished on this computer.\n")
cat("  Variables masked among Wave-9 non-respondents:", length(masked_variables), "\n")
cat("  Category check on the masked data:     passed (at least",
    minimum_category_size, "per category)\n")
cat("  Tasks in the table:        ", nrow(task_table), "\n")
cat("  Task files now present:    ", tasks_done, "\n")
if (tasks_done < nrow(task_table)) {
  cat("  Still to run (on this or other computers):",
      nrow(task_table) - tasks_done, "tasks\n")
} else {
  cat("  All tasks done. Next: File 11b (combine and compare), on one computer",
      "holding every task file.\n")
}
