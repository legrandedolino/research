###############################################################################
# 1. Data loading and merging
#
# Purpose
#   Read a fixed list of variables from the Next Steps SPSS files (UK Data
#   Service study SN 5545), check each file has one row per participant,
#   join every sweep onto the longitudinal design file, and keep the
#   original Wave-1 cohort.
#
# Which variables, and why
#   Only the variables the analysis needs are loaded (data minimization):
#   - design and response: participant ID, sample membership, design weight,
#     PSU, stratum, the Wave-9 final weight, and response outcome at every
#     sweep;
#   - analysis and landmark variables: age-14 background and science and
#     mathematics orientation;
#   - education pathway variables at ages 15, 16, 19-20 and 25, plus monthly
#     activity histories for ages 16-19;
#   - occupation and employment at ages 25 and 32;
#   - the predictors of Wave-8 non-response identified by Silverwood et al.
#     (2024) and listed in the CLS missing-data user guide (Appendix C,
#     Table C1). These are candidate auxiliary variables for multiple
#     imputation (screened in File 6).
#
# Sample
#   The original Wave-1 cohort. The ethnic-minority boost sample added at
#   Wave 4 is excluded because it has no Wave 1-3 measures, following
#   Silverwood et al. (2024).
#
# Reads
#   SPSS files in sav_dir, set in File 0 (headers first, then the selected
#   variables).
#
# Writes
#   data_dir/01_merged.rds               original cohort, variables as
#                                        deposited (participant level)
#   review_dir/01_variable_sources.csv   variable names, source files, roles
#                                        and SPSS labels (metadata only)
#
# References
#   Silverwood RJ, Calderwood L, Henderson M, Sakshaug JW, Ploubidis GB
#   (2024). Longitudinal and Life Course Studies 15(2): 227-250.
#   Silverwood R, et al. (2026). Handling missing data in the CLS cohort
#   studies: User guide (Version 2). UCL Centre for Longitudinal Studies.
###############################################################################

library(haven)                                   # reads SPSS .sav files
suppressPackageStartupMessages(library(dplyr))   # joins


###############################################################################
# Block 1. Paths#
# Folder locations are set once, in File 0, and read here from
# project_paths.rds. Run this script with the working directory set to the
# folder that holds the scripts.
###############################################################################

paths      <- readRDS("project_paths.rds")
sav_dir    <- paths$sav_dir      # SPSS files from the UK Data Service download
data_dir   <- paths$data_dir     # datasets passed between scripts
review_dir <- paths$review_dir   # disclosure-checked outputs


###############################################################################
# Block 2. Source files
#
# One row per SPSS file. file_key is a short name used in the rest of the
# script. The longitudinal design file comes first because every other file
# is joined onto it.
###############################################################################

source_files <- data.frame(
  file_key = c(
    "design",
    "w1_yp", "w1_family",
    "w2_yp",
    "w3_yp", "w3_family",
    "w4_yp",
    "w5_yp", "w5_family",
    "w6_yp", "w7_yp",
    "w4_w7_activity",
    "w8_main", "w8_derived",
    "w9_main", "w9_derived"
  ),
  file_name = c(
    "ns9_2022_longitudinal_file.sav",
    "wave_one_lsype_young_person_2020.sav",
    "wave_one_lsype_family_background_2020.sav",
    "wave_two_lsype_young_person_2020.sav",
    "wave_three_lsype_young_person_2020.sav",
    "wave_three_lsype_family_background_2020.sav",
    "wave_four_lsype_young_person_2020.sav",
    "wave_five_lsype_young_person_2020.sav",
    "wave_five_lsype_family_background_2020.sav",
    "wave_six_lsype_young_person_2020.sav",
    "wave_seven_lsype_young_person_2020.sav",
    "lsype_main_activity_w4-7_nov2011_suppressed.sav",
    "ns8_2015_main_interview.sav",
    "ns8_2015_derived.sav",
    "ns9_2022_main_interview.sav",
    "ns9_2022_derived_variables.sav"
  ),
  description = c(
    "Longitudinal design file: sample, weights, PSU, stratum, response",
    "Wave 1 (age 14) young person interview",
    "Wave 1 (age 14) family background",
    "Wave 2 (age 15) young person interview",
    "Wave 3 (age 16) young person interview",
    "Wave 3 (age 16) family background",
    "Wave 4 (age 17) young person interview",
    "Wave 5 (age 18) young person interview",
    "Wave 5 (age 18) family background",
    "Wave 6 (age 19) young person interview",
    "Wave 7 (age 20) young person interview",
    "Waves 4-7 monthly main-activity histories",
    "Wave 8 (age 25) main interview",
    "Wave 8 (age 25) derived variables",
    "Wave 9 (age 32) main interview",
    "Wave 9 (age 32) derived variables"
  )
)


###############################################################################
# Block 3. Variables to load
#
# One row per variable: the file it comes from and its role in the analysis.
# The participant ID (NSID) is read from every file automatically, and is
# listed here only once, under the design file.
#
# Roles:
#   id / sample / design / response  survey design and response outcomes
#   baseline        age-14 background and orientation (analysis and L14)
#   pathway         later education and activity stages (L16-L25)
#   outcome         employment and occupation at ages 25 and 32
#   nonresponse     CLS Appendix C predictor of non-response (candidate
#                   imputation auxiliary). Appendix C lists 21 predictors:
#                   - 19 are loaded here;
#                   - "previous non-response at Waves 2-7" is derived in
#                     File 2 from W2OUTCOME to W7OUTCOME;
#                   - W4DWPLnk2YP (consent at Wave 4 to pass details to
#                     the Department for Work and Pensions) is not in the
#                     End User Licence release, so it can't be used (see
#                     the Wave 4 rows below).
#
# The Waves 4-7 monthly activity variables are added in Block 5 by name
# pattern, because there is one variable per month.
###############################################################################

variables_to_load <- rbind(

  # Longitudinal design file
  data.frame(
    file_key = "design",
    variable = c("NSID", "MAINBOOST", "DESIGNWEIGHT", "SAMPPSU",
                 "SAMPSTRATUM", paste0("W", 1:9, "OUTCOME"), "W9FINWT"),
    role = c("id", "sample", "design", "design", "design",
             rep("response", 9), "design")
  ),

  # Wave 1 young person: background, orientation, non-response predictors
  data.frame(
    file_key = "w1_yp",
    variable = c("W1sexYP", "W1ethgrpYP", "W1englangYP",
                 "W1lsciYP", "W1gsciYP", "W1favesYP", "W1lmathYP", "W1gmathYP",
                 "W1gowhereYP", "W1namesYP", "W1HcomGYP", "W1alei1YP0h"),
    role = c("baseline; nonresponse", "baseline", "baseline",
             rep("baseline", 5),
             rep("nonresponse", 4))
  ),

  # Wave 1 family background
  data.frame(
    file_key = "w1_family",
    variable = c("W1nssecfam", "W1hiqualgmum", "W1hiqualgdad", "W1inc1est",
                 "W1famtyp2", "IndSchool",
                 "W1hous12HH", "W1condur6MP"),
    role = c(rep("baseline", 6), rep("nonresponse", 2))
  ),

  # Wave 2 (age 15): science course route; non-response predictors.
  # W2SchcontMP is a main-parent item but is deposited in the young person
  # file.
  data.frame(
    file_key = "w2_yp",
    variable = c("W2GenSciYP", "W2ExSciYP", "W2strainYP", "W2SchcontMP"),
    role = c("pathway", "pathway", "nonresponse", "nonresponse")
  ),

  # Wave 3 (age 16): continuation plan, HE intention and intended subject.
  # Subject flags 1-20 are the classified subject areas. Flag 21 records a
  # subject answer of "don't know", which File 2 uses for people with no
  # classified subject (its deposited name is upper case). Flags 22-23
  # (refused, insufficient information) aren't needed: File 2 treats any
  # other answer without a classified subject as item non-response.
  data.frame(
    file_key = "w3_yp",
    variable = c("W3plan16YP", "W3heposs9YP", "W3hesciYP",
                 paste0("W3hesubYP", 1:20), "W3HESUBYP21",
                 "W3cignowYP"),
    role = c(rep("pathway", 3 + 21), "nonresponse")
  ),
  data.frame(
    file_key = "w3_family",
    variable = "W3ageMP",
    role = "nonresponse"
  ),

  # Waves 4 and 5: non-response predictors only.
  # W4DWPLnk2YP (Appendix C) isn't in any End User Licence file. The Wave-5
  # consent items (W5Dwplnk2YP, W5Dwplnk4YP) are NOT substituted for it:
  # they were asked a year later and are different questions.
  data.frame(
    file_key = "w4_yp",
    variable = "W4ClubberYP",
    role = "nonresponse"
  ),
  data.frame(
    file_key = "w5_yp",
    variable = c("W5FairTreat2YP", "W5TeaExpecYP", "W5actYP"),
    role = "nonresponse"
  ),
  data.frame(
    file_key = "w5_family",
    variable = "W5SameaddHH",
    role = "nonresponse"
  ),

  # Waves 6 and 7 (ages 19-20): HE participation and subject group
  data.frame(
    file_key = "w6_yp",
    variable = c("W6HEFlag", "W6HESubGroup2", "W6IagIntroYP0b", "W6SexIntroYP"),
    role = c("pathway", "pathway", "nonresponse", "nonresponse")
  ),
  data.frame(
    file_key = "w7_yp",
    variable = c("W7HEFlag", "W7HESubGroup", "W7SexIntroYP"),
    role = c("pathway", "pathway", "nonresponse")
  ),

  # Wave 8 (age 25): occupation, work status, activity, degree
  data.frame(
    file_key = "w8_main",
    variable = "W8SOC3",
    role = "outcome"
  ),
  # W8DHANVQH (highest academic qualification) confirms "no degree" when
  # W8DDEGP is routed (File 2).
  data.frame(
    file_key = "w8_derived",
    variable = c("W8DWRK", "W8DACTIVITYC", "W8DDEGP", "W8DHESUBGROUP",
                 "W8DHANVQH"),
    role = c("outcome", "pathway", "pathway", "pathway", "pathway")
  ),

  # Wave 9 (age 32): current occupation (SOC2010) and whether working
  data.frame(
    file_key = "w9_main",
    variable = "W9CURRJOB2010",
    role = "outcome"
  ),
  data.frame(
    file_key = "w9_derived",
    variable = "W9DWRK",
    role = "outcome"
  )
)

# Pattern for the monthly activity variables in the Waves 4-7 activity file.
activity_pattern <- "^W7FinAct[0-9]+_B11$"


###############################################################################
# Block 4. Check that every source file exists
#
# If a file is missing, the SPSS file names in the folder are listed so the
# file_name column in Block 2 can be corrected. Only file names are shown.
###############################################################################

source_files$path <- file.path(sav_dir, source_files$file_name)
missing_files <- source_files$file_name[!file.exists(source_files$path)]

if (length(missing_files) > 0) {
  cat("\nSPSS files in", sav_dir, ":\n")
  cat(paste0("  ", list.files(sav_dir, pattern = "\\.sav$")), sep = "\n")
  stop(
    "These files were not found: ", paste(missing_files, collapse = ", "),
    ". Correct their names in Block 2 using the list above."
  )
}


###############################################################################
# Block 5. Check that every variable is in its file (headers only)
#
# Reads only the variable names of each file (no rows). Also adds the
# monthly activity variables, found by name pattern.
#
# If a variable isn't in its assigned file, every SPSS file in the folder is
# searched by name, ignoring upper/lower case. Where it was found is
# reported, and the script stops so Block 3 can be corrected. Only file and
# variable names are shown.
###############################################################################

variables_in_file <- list()
for (i in seq_len(nrow(source_files))) {
  header <- read_sav(source_files$path[i], n_max = 0)
  variables_in_file[[source_files$file_key[i]]] <- names(header)
}

# Add the monthly activity variables.
activity_variables <- grep(activity_pattern,
                           variables_in_file[["w4_w7_activity"]],
                           value = TRUE)
if (length(activity_variables) == 0) {
  stop("No monthly activity variables match ", activity_pattern,
       " in the Waves 4-7 activity file.")
}
variables_to_load <- rbind(
  variables_to_load,
  data.frame(file_key = "w4_w7_activity", variable = activity_variables,
             role = "pathway")
)

# Find listed variables that aren't in their assigned file.
not_found <- data.frame(file_key = character(), variable = character())
for (i in seq_len(nrow(variables_to_load))) {
  key      <- variables_to_load$file_key[i]
  variable <- variables_to_load$variable[i]
  if (!(variable %in% variables_in_file[[key]])) {
    not_found <- rbind(not_found, data.frame(file_key = key, variable = variable))
  }
}

# Every file must contain the participant ID.
for (key in source_files$file_key) {
  if (!("NSID" %in% variables_in_file[[key]])) {
    stop("File '", key, "' has no NSID variable.")
  }
}

if (nrow(not_found) > 0) {
  cat("\nVariables not found in their assigned file:\n")
  all_sav_files <- list.files(sav_dir, pattern = "\\.sav$")

  for (i in seq_len(nrow(not_found))) {
    variable <- not_found$variable[i]
    found_in <- character()
    similar  <- character()

    for (sav_file in all_sav_files) {
      names_here <- names(read_sav(file.path(sav_dir, sav_file), n_max = 0))

      # Exact name, ignoring upper/lower case.
      match_here <- names_here[tolower(names_here) == tolower(variable)]
      if (length(match_here) > 0) found_in <- c(found_in, paste0(sav_file, " (as ", match_here[1], ")"))

      # Similar names (up to about two characters different), in case the
      # variable was renamed in the deposit.
      near_here <- agrep(variable, names_here, max.distance = 0.2,
                         ignore.case = TRUE, value = TRUE)
      near_here <- setdiff(near_here, match_here)
      if (length(near_here) > 0) similar <- c(similar, paste0(near_here, " in ", sav_file))
    }

    cat("  ", variable, " (assigned to '", not_found$file_key[i], "'): ",
        if (length(found_in) == 0) "no exact match in any file"
        else paste(found_in, collapse = "; "),
        "\n", sep = "")
    if (length(found_in) == 0 && length(similar) > 0) {
      cat("      similar names: ", paste(similar, collapse = "; "), "\n",
          sep = "")
    }
  }

  stop("Correct the file_key or variable spelling in Block 3 (and add any ",
       "new file to Block 2) using the list above, then run again.")
}


###############################################################################
# Block 6. Read the selected variables from each file
#
# user_na = TRUE keeps the deposited negative codes (e.g. "refused", "don't
# know", "not applicable") as labelled values instead of turning them into
# NA. Their meaning differs by sweep: -1 is "don't know" at Waves 1-7 but
# "not applicable" in the Wave 8-9 derived variables. File 2 uses the codes
# to keep meaningful categories (such as family NS-SEC "not applicable")
# apart from missing values.
#
# Each file must have exactly one row per participant, with no missing or
# duplicated NSID.
###############################################################################

files_read <- list()
for (i in seq_len(nrow(source_files))) {
  key <- source_files$file_key[i]
  columns <- unique(c("NSID",
                      variables_to_load$variable[variables_to_load$file_key == key]))

  file_data <- read_sav(source_files$path[i], user_na = TRUE,
                        col_select = all_of(columns))

  if (anyNA(file_data$NSID) || anyDuplicated(file_data$NSID) > 0) {
    stop("File '", key, "' does not have exactly one row per NSID.")
  }

  files_read[[key]] <- file_data
}


###############################################################################
# Block 7. Record variable sources and labels
#
# A metadata table: variable name, file, role and the SPSS variable label.
#
# Disclosure check: this table holds variable metadata only, with no
# participant data, so there is nothing to suppress.
###############################################################################

variable_sources <- variables_to_load
variable_sources$file_name <- source_files$file_name[
  match(variable_sources$file_key, source_files$file_key)
]
variable_sources$spss_label <- NA_character_

for (i in seq_len(nrow(variable_sources))) {
  column <- files_read[[variable_sources$file_key[i]]][[variable_sources$variable[i]]]
  label  <- attr(column, "label", exact = TRUE)
  if (!is.null(label)) variable_sources$spss_label[i] <- as.character(label)
}

variable_sources$disclosure_checked <- TRUE
write.csv(variable_sources, file.path(review_dir, "01_variable_sources.csv"),
          row.names = FALSE, na = "")


###############################################################################
# Block 8. Join all files onto the design file
#
# The longitudinal design file lists everyone in the study, so it's the
# spine. Every other file is left-joined by NSID. Participants absent from a
# sweep's file get missing values for that sweep's variables, and File 3
# classifies them as wave non-response. The join must not add or lose rows,
# and no variable may appear in two files.
###############################################################################

merged <- files_read[["design"]]

for (key in setdiff(source_files$file_key, "design")) {
  incoming <- files_read[[key]]

  shared_columns <- setdiff(intersect(names(merged), names(incoming)), "NSID")
  if (length(shared_columns) > 0) {
    stop("Variables appear in more than one file (", key, "): ",
         paste(shared_columns, collapse = ", "))
  }

  rows_before <- nrow(merged)
  merged <- left_join(merged, incoming, by = "NSID")
  if (nrow(merged) != rows_before) {
    stop("Joining '", key, "' changed the number of rows.")
  }
}


###############################################################################
# Block 9. Keep the original Wave-1 cohort
#
# MAINBOOST = 1 marks the original main sample. The Wave-4 ethnic-minority
# boost is excluded because its members weren't eligible for the Wave 1-3
# measures, so their early data are structurally absent, not missing.
###############################################################################

in_original_cohort <- as.numeric(merged$MAINBOOST) == 1
if (anyNA(in_original_cohort)) {
  stop("MAINBOOST has missing values. Check the design file.")
}

merged_original <- merged[in_original_cohort, ]


###############################################################################
# Block 10. Save
###############################################################################

saveRDS(merged_original, file.path(data_dir, "01_merged.rds"), compress = "xz")


###############################################################################
# Block 11. Summary
#
# Whole-sample counts and variable counts only, which are safe to share.
###############################################################################

cat("\nFile 1 finished.\n")
cat("  Files read:                     ", nrow(source_files), "\n")
cat("  Variables loaded (incl. NSID):  ", ncol(merged_original), "\n")
cat("    of which monthly activity:    ", length(activity_variables), "\n")
cat("    of which non-response preds.: ",
    sum(grepl("nonresponse", variables_to_load$role)), "\n")
cat("  Participants in design file:    ", nrow(merged), "\n")
cat("  Original Wave-1 cohort kept:    ", nrow(merged_original), "\n")
cat("  Saved:", file.path(data_dir, "01_merged.rds"), "\n")
