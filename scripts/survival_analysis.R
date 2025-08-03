
# about the script --------------------------------------------------------
#purpose : Survival analysis of viral rebound for CTC clients in fisher fork zone Mwanza, Tanzania
#Author : Lawrence
#Date : 2025-07-21

# installing and loading packages -----------------------------------------

pacman::p_load(
  rio,     #importing data
  here,    #for relative file path
  skimr,   #for reviewing data
  janitor, #for cleaning data
  epikit,  #for creating age categories
  flextable, #making publication ready tables
  gtsummary, #making publication ready tables
  scales,  # percentage in tables
  ggsurvfit, 
  survival,
  gtsummary, 
  tidycmprsk,
  broom,
  tidyverse #for data management and visualization
)


# import data and explore -------------------------------------------------------------
rebound_data <- import(here("data", "raw_data.csv"))
skimr::skim(rebound_data)
glimpse(rebound_data)


# data cleaning -----------------------------------------------------------
survival_data <- rebound_data %>%
  #clean names
  clean_names() %>%
  #rename variables
  rename(
    months_on_art = duration_on_art_months,
    time_to_rebound = dayz_rebound,
    event = rebound_ovall
  ) %>%
  #selecting column to work with
 select(lookup_id, facility_name, sex, age,
        marital_status, residence, residence_lvl,
        distance_hf, stigma, mental_illness, 
        alcohol_use, cd4_test_date, cd4_counts, 
        who_clinical_stage, art_adherence, art_regimen,
        arv_code, arv_combination, hiv_care_appointment, 
        number_of_viral_load_test, low_level_viremia_at_basline,
        months_on_art, art_dispensing_day, 
        hiv_tb_co_infection, tb_preventive_therapy, 
        interruption_in_art_treatment, 
        client_category, facility_ownership, 
        type_of_health_facility, 
        facility_patient_volume, 
        facility_staffing, art_refill_model,
        last_hvl_date_2019, hvl_results_2019,last_hvl_date_2022,
        hvl_results_2022,event, time_to_rebound) %>%
  #date columns
  mutate(last_hvl_date_2019 = dmy(last_hvl_date_2019)) %>%
  mutate(last_hvl_date_2022 = dmy(last_hvl_date_2022)) %>%
  #total time in study for each participant
  mutate(total_time_in_study =
           as.double(last_hvl_date_2022 - last_hvl_date_2019)) %>%
  #merging the time variable account individual study times
  mutate(time = coalesce(time_to_rebound, total_time_in_study)) %>%
  mutate(time_years = time/365.25) %>%      #time in years
  mutate(time_months = time/30) %>%         #time in months
  #alternatively assuming all participants made it to the end of study and started at the same time
 mutate(time_final = 
          ifelse(is.na(time_to_rebound) & event == 0, 1095, time)) %>%
  
  # creating age/number groups
  mutate(
    age_cat = case_when( 
      age >= 15 & age < 25   ~ "15-24",
      age >= 25 & age < 44 ~ "25-44",
      age >= 45 & age < 64 ~ "45-64",
      age >= 65   ~ "65+")
  ) %>%
  #distance category
  mutate(distance_cat = case_when(
    distance_hf <= 10 ~ "0-10 Km",
    distance_hf > 10   ~ "More than 10 Km")
  ) %>%
  
  #number of dispensing days
  mutate(dispensing_cat =case_when(
    art_dispensing_day <= 30 ~ "30 days or less",
    art_dispensing_day > 30 ~ "More than 30 days")
           ) %>% 
  #ART years
  mutate(arv_years = months_on_art/12) %>% 
  #Duration on ART (months)
  mutate(duration_cat = case_when(
    arv_years < 5 ~ "less than 5 Years",
    arv_years >= 5 & arv_years <= 8 ~ "5-8 Years",
    arv_years >= 8 ~ "more than 8 Years")
         ) %>%

  # go through all columns replace empty cells with NA
  mutate(across(.cols = where(is.character),
                .fns = ~na_if(.x, ""))) %>%
  #change character variables as factors
  mutate(across(.cols = where(is.character), as.factor)) %>%
  # recode values
  mutate(who_clinical_stage = recode(who_clinical_stage,
         "1" = "clinical stage 1",
         "2" = "clinical stage 2",
         "3" = "clinical stage 3",
         "4" = "clinical stage 4")) %>%
  mutate(DTG_cat = ifelse(grepl("DTG", arv_combination),
                          "DTG based", "Non-DTG based")) %>%
  mutate(appointment_cat = recode(hiv_care_appointment,
    "Community ART Refill" = "Scheduled visit",
    "Scheduled clinic visit" = "Scheduled visit",
    "Treatment supporter drug pick up" = "Scheduled visit",
    "Traced back after LTFU" = "Unscheduled visit",
    "Unscheduled clinic visit" = "Unscheduled visit"))



# summary table -----------------------------------------------------------


# survival analysis -------------------------------------------------------
#creating a survival object
surv_obj <-Surv(time = survival_data$time,
       event = survival_data$event)
surv_fit <- survfit(surv_obj ~1)
summary(surv_fit, times = 365.25)
plot(surv_fit)

#survival objects in years
year_object <- Surv(time = survival_data$time_years,
                    event = survival_data$event)

year_fit <- survfit(year_object ~1)
summary(survfit(year_object ~1))

#survival object by months
month_object <- Surv(time = survival_data$time_months,
                    event = survival_data$event)

month_fit <- survfit(month_object ~1)
summary(survfit(month_object ~1))

# incidence rates ---------------------------------------------------------


#Incidence rates
overall_icidence <- survival_data %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    incidence_rate = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
    )
#sex IR
survival_data %>%
  group_by(sex) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )


#residence IR
survival_data %>%
  group_by(residence_lvl) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )

#age category
survival_data %>%
  group_by(age_cat) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )

#Marital status
survival_data %>%
  group_by(marital_status) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )

#client category
survival_data %>%
  group_by(client_category) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )

#ART adherence
survival_data %>%
  group_by(art_adherence) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )

#DTG category
survival_data %>%
  group_by(DTG_cat) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )

#low_level_viremia_at_basline
survival_data %>%
  group_by(low_level_viremia_at_basline) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )



#Duration category
survival_data %>%
  group_by(duration_cat) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )


#ART dispensing
survival_data %>%
  group_by(dispensing_cat) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )


#WHO clinical stage
survival_data %>%
  group_by(who_clinical_stage) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )


#distance category
survival_data %>%
  group_by(distance_cat) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )


#ART refill model
survival_data %>%
  group_by(art_refill_model) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )


#appointment category
survival_data %>%
  group_by(appointment_cat) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )
























#export clean data







# cox proportional hazard -------------------------------------------------
coxph(surv_obj ~ age_cat, data = survival_data) %>%
  tbl_regression(exp = TRUE)

coxph(surv_obj ~ client_category, data = survival_data) %>%
  tbl_regression(exp = TRUE)

coxph(surv_obj ~ dispensing_cat, data = survival_data) %>%
  tbl_regression(exp = TRUE)
coxph(surv_obj ~ distance_cat, data = survival_data) %>%
  tbl_regression(exp = TRUE)


coxph(surv_obj ~ appointment_cat, data = survival_data) %>%
  tbl_regression(exp = TRUE)


coxph(surv_obj ~ art_refill_model, data = survival_data) %>%
  tbl_regression(exp = TRUE)

#multivarible cox model

coxph(surv_obj ~ age_cat + client_category +  dispensing_cat
      +distance_cat + appointment_cat + 
        art_refill_model, data = survival_data) %>%
  tbl_regression(exp = TRUE)



































#Testing area ----------------------
surv_obj_02 <-Surv(time = survival_data$time,
                event = survival_data$event)

survfit(surv_obj_02 ~1)
plot(survfit(surv_obj_02 ~1))
tabyl(rebound_data$Rebound_2019)

tabyl(survival_data$client_category) %>%
  adorn_pct_formatting()


tabyl(survival_data$art_adherence) %>%
  adorn_pct_formatting()
tabyl(survival_data$art_regimen) %>%
  adorn_pct_formatting()
tabyl(survival_data$DTG_cat) %>%
  adorn_pct_formatting()

tabyl(survival_data$low_level_viremia_at_basline) %>%
  adorn_pct_formatting()

tabyl(survival_data$duration_cat) %>%
  adorn_pct_formatting()

tabyl(survival_data$dispensing_cat) %>%
  adorn_pct_formatting()

tabyl(survival_data$who_clinical_stage) %>%
  adorn_pct_formatting()

tabyl(survival_data$distance_cat) %>%
  adorn_pct_formatting()

tabyl(survival_data$art_refill_model) %>%
  adorn_pct_formatting()


tabyl(survival_data$hiv_care_appointment) %>%
  adorn_pct_formatting()

#specifying the reference in the model

coxph(surv_obj ~ relevel(factor(appointment_cat), ref = "unscheduled visit"), 
      data = survival_data) %>% 
  tbl_regression(exp = TRUE)


survival_data %>%
  group_by(sex) %>%
  summarise(
    events = sum(event == 1),
    total_time = sum(time_months),
    IR = (events / total_time) * 1000,
    # Exact Poisson 95% CI - RECOMMENDED for survival analysis
    ci_lower = (qchisq(0.025, 2 * events) / (2 * total_time)) * 1000,
    ci_upper = (qchisq(0.975, 2 * (events + 1)) / (2 * total_time)) * 1000
  )

