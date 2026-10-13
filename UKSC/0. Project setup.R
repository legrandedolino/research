###############################################################################
# 0. Project setup
#
# Purpose
#   Set the folder locations for the whole pipeline (the only place paths are
#   written), check this computer's R version and required packages, record
#   package versions and session information, and create the output folders.
#   This script loads NO participant data.
#
# How to run the pipeline
#   Set the working directory to the folder that holds these scripts, then
#   run the scripts in order, starting with this one. For example:
#     cd "<folder with the scripts>"
#     Rscript "0. Project setup.R"
#     Rscript "1. Data loading and merging.R"
#   Every later script reads the folder locations from project_paths.rds,
#   which this script saves in the scripts folder.
#
# When to run this script
#   First, on every computer that will run any script in this pipeline, and
#   again after updating R or any package or changing a folder location.
#   Heavy stages can be split across several computers (Files 8a, 9a, 11a,
#   13a, 17a and 19a). Their results are combined only when every computer
#   involved reports the same R and package versions.
#
# Reads
#   Nothing.
#
# Writes
#   project_paths.rds (in the scripts folder): folder locations only
#   output_dir/review (software information only, no participant data):
#     00_package_versions_<computer>.csv
#     00_session_info_<computer>.txt
###############################################################################


###############################################################################
# Block 1. Paths
#
# The only place in the pipeline where folder locations are written. Set
# these two to the locations on your computer:
#   raw_data_dir  the UK Data Service download of SN 5545 (SPSS format),
#                 held under the End User Licence
#   output_dir    where every output is written. It must be outside the
#                 folder that holds these scripts, so participant data never
#                 sits next to the code.
###############################################################################

raw_data_dir <- "/Users/ledolino/Documents/UK Next Steps/UKDA-5545-spss"
output_dir   <- "/Users/ledolino/Documents/UK Next Steps/Private data/v2"

# Derived locations (these don't need changing):
sav_dir    <- file.path(raw_data_dir, "spss", "spss28", "safeguarded_eul")
                                                # SPSS files within the download
data_dir   <- file.path(output_dir, "data")     # datasets passed between scripts
review_dir <- file.path(output_dir, "review")   # disclosure-checked outputs
tasks_dir  <- file.path(output_dir, "tasks")    # one file per task in split runs

# The scripts must be run from their own folder (see "How to run" above).
if (!file.exists("0. Project setup.R")) {
  stop("Set the working directory to the folder that holds the scripts ",
       "(it should contain \"0. Project setup.R\").")
}


###############################################################################
# Block 2. R version
#
# R 4.3.0 or later is required. The scripts use the native pipe (|>) and
# other features of recent R versions.
###############################################################################

minimum_r_version <- "4.3.0"

if (getRversion() < minimum_r_version) {
  stop(
    "This project needs R ", minimum_r_version, " or later. ",
    "This computer has R ", getRversion(), "."
  )
}


###############################################################################
# Block 3. Required packages
#
# One row per package, with what it's used for. A minimum version is set only
# where a script depends on a feature added in that version; otherwise any
# version is accepted and recorded in Block 8.
###############################################################################

required_packages <- data.frame(
  package = c("haven", "dplyr", "survey", "mice", "glmnet", "pmsampsize", "lme4"),
  used_for = c(
    "Reading the SPSS files from the UKDS deposit (File 1)",
    "Readable variable construction with case_when() (File 2 onward)",
    "Survey-weighted estimation with PSU, stratum and weight (Files 5, 12)",
    "Multiple imputation by chained equations (Files 7-13, 16-17)",
    "Elastic-net logistic regression (Files 16, 17, 19)",
    "Riley et al. (2020) sample-size calculation (File 15)",
    "Random-intercept model for differences between schools (File 18)"
  ),
  minimum_version = c(NA, NA, NA, "3.12.0", NA, NA, NA)
)
# mice 3.12.0 introduced the `ignore` argument, which File 17a uses to fit
# imputation models on training rows only.

# Stop if any package is missing, and show the command that installs them.
installed <- rownames(installed.packages())
missing_packages <- setdiff(required_packages$package, installed)

if (length(missing_packages) > 0) {
  stop(
    "These packages are not installed: ",
    paste(missing_packages, collapse = ", "), ".\n",
    "Install them with:\n",
    "  install.packages(c(\"", paste(missing_packages, collapse = "\", \""),
    "\"))\n",
    "then run this script again."
  )
}

# Stop if any package is older than its minimum version.
for (i in seq_len(nrow(required_packages))) {
  package_name <- required_packages$package[i]
  minimum      <- required_packages$minimum_version[i]

  if (!is.na(minimum) && packageVersion(package_name) < minimum) {
    stop(
      package_name, " ", minimum, " or later is needed (",
      required_packages$used_for[i], "). This computer has ",
      packageVersion(package_name), "."
    )
  }
}


###############################################################################
# Block 4. Raw data folder
#
# Only checks that the folder of SPSS files exists. Nothing is read or
# listed; File 1 reads the data.
###############################################################################

if (!dir.exists(sav_dir)) {
  stop(
    "The SPSS data folder was not found:\n  ", sav_dir, "\n",
    "Set raw_data_dir in Block 1 to where the SN 5545 download is stored on ",
    "this computer. The data must be held under the UKDS End User Licence ",
    "and your institution's security rules."
  )
}


###############################################################################
# Block 5. Output folders
#
# Creates output_dir and its subfolders if they are missing. Existing folders
# and files are left untouched.
###############################################################################

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tasks_dir,  recursive = TRUE, showWarnings = FALSE)

if (!all(dir.exists(c(output_dir, data_dir, review_dir, tasks_dir)))) {
  stop("Could not create the output folders under:\n  ", output_dir)
}


###############################################################################
# Block 6. Save the folder locations for the other scripts
#
# Every later script starts with readRDS("project_paths.rds"), so the
# locations above are written in one place only. The file holds folder
# locations only, no data.
###############################################################################

project_paths <- list(
  raw_data_dir = raw_data_dir,
  sav_dir      = sav_dir,
  output_dir   = output_dir,
  data_dir     = data_dir,
  review_dir   = review_dir,
  tasks_dir    = tasks_dir
)
saveRDS(project_paths, "project_paths.rds")


###############################################################################
# Block 7. Identify this computer
#
# The computer name is used in the output file names, so records from
# different computers don't overwrite each other. Characters that aren't
# safe in file names are replaced with "_".
###############################################################################

computer_name <- gsub("[^A-Za-z0-9_-]", "_", Sys.info()[["nodename"]])


###############################################################################
# Block 8. Record R and package versions
#
# Disclosure check: this table holds software information only, with no
# participant data, so there is nothing to suppress. The disclosure_checked
# column is included to match the other output tables.
###############################################################################

# Look up each package's version as text, e.g. "3.19.0".
installed_versions <- character(nrow(required_packages))
for (i in seq_len(nrow(required_packages))) {
  installed_versions[i] <- format(packageVersion(required_packages$package[i]))
}

package_versions <- data.frame(
  computer = computer_name,
  package  = c("R", required_packages$package),
  version  = c(format(getRversion()), installed_versions),
  recorded_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  disclosure_checked = TRUE
)

versions_file <- file.path(
  review_dir, paste0("00_package_versions_", computer_name, ".csv")
)
write.csv(package_versions, versions_file, row.names = FALSE)


###############################################################################
# Block 9. Record full session information
#
# Includes the operating system, locale, random-number generator settings and
# every loaded package, which helps anyone reproducing the results
# (REFORMS item 2c).
###############################################################################

session_file <- file.path(
  review_dir, paste0("00_session_info_", computer_name, ".txt")
)
writeLines(
  c(
    paste("Computer:", computer_name),
    paste("Recorded:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Random-number generator:", paste(RNGkind(), collapse = ", ")),
    "",
    capture.output(sessionInfo())
  ),
  session_file
)


###############################################################################
# Block 10. Compare versions across computers
#
# If other computers have also written version files to review_dir (for
# example on shared, institution-approved storage), check they all match.
# A mismatch doesn't stop this script. It is a warning that split runs from
# these computers must not be combined until versions are aligned.
###############################################################################

version_files <- list.files(
  review_dir, pattern = "^00_package_versions_.*\\.csv$", full.names = TRUE
)
# Read every column as text, so versions like "5.0" stay "5.0" and aren't
# converted to the number 5.
all_versions <- data.frame()
for (version_file in version_files) {
  all_versions <- rbind(all_versions, read.csv(version_file, colClasses = "character"))
}

# For each package, count how many different versions appear across
# computers. More than one means a mismatch.
version_counts <- table(unique(all_versions[c("package", "version")])$package)
mismatched_packages <- names(version_counts)[version_counts > 1]

if (length(mismatched_packages) > 0) {
  warning(
    "Computers disagree on versions of: ",
    paste(mismatched_packages, collapse = ", "), ". ",
    "Align them before combining split runs. See the 00_package_versions_* ",
    "files in review_dir."
  )
}


###############################################################################
# Block 11. Summary
#
# Console output contains software information only and is safe to share.
###############################################################################

cat("\nFile 0 finished.\n")
cat("  Computer:            ", computer_name, "\n")
cat("  R version:           ", as.character(getRversion()), "\n")
cat("  Packages checked:    ", nrow(required_packages), "(all installed)\n")
cat("  Computers on record: ", length(unique(all_versions$computer)), "\n")
cat("  Version mismatches:  ",
    if (length(mismatched_packages) == 0) "none" else
      paste(mismatched_packages, collapse = ", "), "\n")
cat("  Output folders ready under:", output_dir, "\n")
cat("  Folder locations saved to:  project_paths.rds (read by every later script)\n")
