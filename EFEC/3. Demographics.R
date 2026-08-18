################################################################################
#
# 3. DEMOGRAPHICS - EDUCATING FOR ENVIRONMENTAL CHANGE
# Latest Version: August 18, 2026
#
# Summarizes participant demographics, teaching context, and institute feedback.
#
################################################################################

################################################################################
# LOADING IN DATA
################################################################################
setwd("~/Documents/EfEC")

library(dplyr)
library(gtsummary)
library(flextable)

load("data_clean.RData")

theme_gtsummary_journal(journal = "jama")
theme_gtsummary_compact()

################################################################################
# WHO PARTICIPATED IN THE SUMMER SCIENCE INSTITUTE?
################################################################################

demographics_data <- data_clean %>%
  filter(test == "pre")

# Main demographics and education
t_demo <- demographics_data %>%
  select(
    year, gender, race, previousParticipant,
    starts_with("yearsTeaching"), starts_with("education"), starts_with("degrees")
  ) %>%
  tbl_summary(
    by = year,
    statistic = list(
      all_categorical() ~ "{n} ({p}%)",
      all_continuous() ~ "{mean} ({sd})"
    ),
    label = list(
      gender ~ "Gender",
      race ~ "Race/Ethnicity",
      previousParticipant ~ "Previous EfEC Participant",
      yearsTeaching ~ "Years of Teaching",
      education ~ "Highest Education Completed",
      degrees_generalEd ~ "Completed a degree in General Education",
      degrees_lifeScience ~ "Completed a degree in Life Science",
      degrees_earthSpace ~ "Completed a degree in Earth or Space Science"
    ),
    missing = "ifany",
    missing_text = "Missing"
  ) %>%
  add_overall(last = TRUE) %>%
  bold_labels()

# Teaching levels
t_teach <- demographics_data %>%
  select(year, subElem, subMidLife, subMidEar, subHiLife, subHiEar, subOther) %>%
  tbl_summary(
    by = year,
    type = all_dichotomous() ~ "dichotomous",
    label = list(
      subElem ~ "Elementary",
      subMidLife ~ "Middle School Life Science",
      subMidEar ~ "Middle School Earth Science",
      subHiLife ~ "High School Life Science",
      subHiEar ~ "High School Earth Science",
      subOther ~ "Other"
    ),
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  modify_table_body(~ .x %>% mutate(row_type = "level"))

# Class environment
t_class <- demographics_data %>%
  select(year, indGen, indInc, indSpE, indESL, indAdv) %>%
  tbl_summary(
    by = year,
    type = all_dichotomous() ~ "dichotomous",
    label = list(
      indGen ~ "General / Regular",
      indInc ~ "Inclusion",
      indSpE ~ "Special Education",
      indESL ~ "ESL",
      indAdv ~ "Advanced / Gifted"
    ),
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  modify_table_body(~ .x %>% mutate(row_type = "level"))

table1_stacked <- tbl_stack(
  tbls = list(t_demo, t_teach, t_class),
  group_header = c(
    "Demographics and Experience",
    "Expected Teaching Level for Upcoming Year",
    "Expected Teaching Environment for Upcoming Year"
  )
) %>%
  modify_header(label = "**Participant Characteristics**") %>%
  modify_caption("**Table 1. Demographic and Professional Characteristics of Participants by Year**")

table1_stacked

table1_stacked %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "Table_1_Participant_Demographics.docx")

################################################################################
# WHAT IS THEIR TEACHING ENVIRONMENT?
################################################################################

environment_data <- data_clean %>%
  filter(test == "pre") %>%
  select(
    year,
    starts_with("schoolRate"),
    starts_with("studentSocioStatus")
  )

table2 <- environment_data %>%
  tbl_summary(
    by = year,
    statistic = list(
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      "schoolRate_strategies" ~ "School's strategies for teaching about climate change",
      "schoolRate_resources" ~ "School's resources for teaching about climate change",
      "schoolRate_k12network" ~ "School's support for networking with colleagues who teach about climate change",
      "schoolRate_scientists" ~ "School's support for networking  with scientists who conduct research on climate change",
      "schoolRate_credits" ~ "School's support for professional development credits",
      "schoolRate_understanding" ~ "School's understanding of climate change and its impacts",
      starts_with("studentSocioStatus") ~ "Majority Socioeconomic Status of Students"
    ),
    missing = "ifany",
    missing_text = "Missing/Did not answer"
  ) %>%
  add_overall(last = TRUE) %>%
  bold_labels() %>%
  modify_header(label = "**School Context**") %>%
  modify_caption("**Table 2. Self-Reported Perceptions on Teaching Environment
                 and School Context of Participants by Year**")

table2

table2 %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "Table_2_Teaching_Environment.docx")

################################################################################
# WHO ARE THEIR BASELINE ATTITUDES AND KNOWLEDGE?
################################################################################

baseline_data <- data_clean %>%
  filter(test == "pre")

# Climate-change understanding (1-6 scale)
t_under <- baseline_data %>%
  select(year, starts_with("ccUnder")) %>%
  tbl_summary(
    by = year,
    label = list(
      "ccUnder_climChange" ~ "How the climate is changing",
      "ccUnder_humanCaused" ~ "How humans are causing climate change",
      "ccUnder_policy" ~ "How policy decisions impact climate change",
      "ccUnder_stepsPeople" ~ "The steps people can take to limit climate change",
      "ccUnder_mitigateEffects" ~ "What is needed to mitigate the effects of climate change"
    ),
    type = everything() ~ "continuous",
    statistic = all_continuous() ~ "{mean} ({sd})",
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  modify_table_body(~ .x %>% mutate(row_type = "level"))

# Teaching efficacy (1-6 scale)
t_eff <- baseline_data %>%
  select(year, starts_with("Efficacy")) %>%
  tbl_summary(
    by = year,
    label = list(
      "Efficacy_1" ~ "Can find effective ways to teach about climate change",
      "Efficacy_2" ~ "Cannot teach about climate change as well as other topics",
      "Efficacy_3" ~ "Understand the strategies to effectively teach climate change concepts",
      "Efficacy_4" ~ "Do not believe that they are effective when facilitating climate change education activities",
      "Efficacy_5" ~ "Cannot teach about climate change effectively",
      "Efficacy_6" ~ "Understand climate change concepts well enough to be effective in teaching those concept",
      "Efficacy_7" ~ "Cannot effectively explain climate change concepts to students.",
      "Efficacy_8" ~ "Will not invite the principal to evaluate their climate change education teaching",
      "Efficacy_9" ~ "Do not know how to help a student having difficulty understanding a climate change concept",
      "Efficacy_10" ~ "Encourage student questions when teaching about climate change",
      "Efficacy_11" ~ "Do not know how to motivate students about climate change concepts"
    ),
    type = everything() ~ "continuous",
    statistic = all_continuous() ~ "{mean} ({sd})",
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  modify_table_body(~ .x %>% mutate(row_type = "level"))

# Topic familiarity (1-4 scale)
t_topic <- baseline_data %>%
  select(year, starts_with("topic")) %>%
  tbl_summary(
    by = year,
    label = list(
      "topic_climateCauses" ~ "The climate and what causes climate",
      "topic_modeling" ~ "Climate modeling and simulations",
      "topic_misconceptions" ~ "Addressing misconceptions about climate change with students",
      "topic_evidence" ~ "Scientific evidence of climate change",
      "topic_effects" ~ "The effects of climate change",
      "topic_geoengineering" ~ "Geoengineering",
      "topic_stories" ~ "Integrating storytelling into environmental education",
      "topic_inequity" ~ "How climate change impacts some groups of people more than others"
    ),
    type = everything() ~ "continuous",
    statistic = all_continuous() ~ "{mean} ({sd})",
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  modify_table_body(~ .x %>% mutate(row_type = "level"))

# Perceived challenges (1-5 scale)
t_challenge <- baseline_data %>%
  select(year, starts_with("challenge"), -challenge_Qual, -challenge_knowledge) %>%
  mutate(across(starts_with("challenge"), as.numeric)) %>%
  tbl_summary(
    by = year,
    label = list(
      "challenge_understand" ~ "Personal understanding of climate change science, impacts, and mitigation and adaptation strategies",
      "challenge_standards" ~ "Connecting climate change to state teaching standards",
      "challenge_resources" ~ "Accessing adequate resources to teach climate change that are grade level appropriate",
      "challenge_unsupportive" ~ "Unsupportive school environment and/or parents",
      "challenge_misconceptions" ~ "Addressing misconceptions about climate change with my students",
      "challenge_controversial" ~ "Facilitating discussions about controversial topics such as climate injustice with my students",
      "challenge_time" ~ "Having adequate time to cover the topic of climate change properly",
      "challenge_multipleClass" ~ "Working with other teachers to integrate climate change into multiple classroom curriculums",
      "challenge_myCurriculum" ~ "Understanding where climate change fits in my curriculum"
    ),
    type = everything() ~ "continuous",
    statistic = all_continuous() ~ "{mean} ({sd})",
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  modify_table_body(~ .x %>% mutate(row_type = "level"))

# Perceived controversy
t_controv <- baseline_data %>%
  select(year, starts_with("ccControversial")) %>%
  tbl_summary(
    by = year,
    label = list(
      "ccControversial" = "Perceived Controversy of Climate Change in School District"
    ),
    statistic = all_categorical() ~ "{n} ({p}%)",
    missing = "no"
  ) %>%
  add_overall(last = TRUE)

table3_stacked <- tbl_stack(
  tbls = list(t_under, t_eff, t_topic, t_challenge, t_controv),
  group_header = c(
    "(Understanding) Agreement to statements that teachers understand [Mean (SD); 1=Strongly Disagree to 6=Strongly Agree]",
    "(Efficacy) Agreement to statements that teachers [Mean (SD); 1=Strongly Disagree to 6=Strongly Agree]",
    "Familiarity to Workshop Topics [Mean (SD); 1=Unfamiliar to 4=Know a lot]",
    "Perceived Challenges [Mean (SD); 1=Not at all to 5=Extremely]",
    ""
  )
) %>%
  modify_header(label = "**Baseline Construct & Items**") %>%
  modify_caption("**Table 3. Baseline Attitudes, Knowledge, and Perceived Challenges by Year**")

table3_stacked

table3_stacked %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "Table_3_Baseline_Attitudes.docx")

################################################################################
# HOW DID THEY EVALUATE THE INSTITUTE?
################################################################################

evaluation_data <- data_clean %>%
  filter(test == "post")

# Workshop objectives
t_obj <- evaluation_data %>%
  select(year, starts_with("workshopObj")) %>%
  tbl_summary(
    by = year,
    label = list(
      "workshopObj_statedClearly" ~ "The workshop objectives were clearly stated.",
      "workshopObj_relevant" ~ "The information and activities were relevant to my grade-level curricula.",
      "workshopObj_expectations" ~ "The workshop met my expectations.",
      "workshopObj_preparedPresent" ~ "The presenters were well prepared."
    ),
    statistic = all_categorical() ~ "{n} ({p}%)",
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  bold_labels()

# Workshop rating
t_rating <- evaluation_data %>%
  select(year, starts_with("workshopRating")) %>%
  tbl_summary(
    by = year,
    label = list(
      "workshopRating" ~ "Overall rating of the Summer Science Institute"
    ),
    statistic = all_categorical() ~ "{n} ({p}%)",
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  bold_labels()

# Workshop activities
t_act <- evaluation_data %>%
  select(year, starts_with("workshopAct")) %>%
  tbl_summary(
    by = year,
    label = list(
      "workshopAct_lectures" ~ "Lectures and informational presentations",
      "workshopAct_groupActive" ~ "Instructor-led activities completed in groups",
      "workshopAct_individualActive" ~ "Instructor-led activities completed individually",
      "workshopAct_groupDiscuss" ~ "Group discussions",
      "workshopAct_tutorials" ~ "Tutorials",
      "workshopAct_QandA" ~ "Question and answer sessions"
    ),
    statistic = all_categorical() ~ "{n} ({p}%)",
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  bold_labels()

table4_stacked <- tbl_stack(
  tbls = list(t_obj, t_rating, t_act),
  group_header = c(
    "Workshop Objectives",
    "Workshop Ratings",
    "Helpfulness of Activities"
  )
) %>%
  modify_header(label = "**Evaluation Construct & Items**") %>%
  modify_caption("**Table 4. Post-Institute Evaluations by Year**")

table4_stacked

table4_stacked %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "Table_4_Evaluations.docx")
