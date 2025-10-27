#Question 1:
#Study the association between proteomics and time to death.
#Compute a Protein Risk Score. Fit a linear model using all proteins as predictors and extract the model predictions.
#Fit a linear mixed-effects model, including age, sex and the Protein Risk Score computed at the previous step and add a random effect linked to Hospitals.



# Load necessary libraries
library(lme4)
library(dplyr)
library(tidyr)

# Load data
load("June_homework_data.RData")  

#Rename for convenience 
data<- June_homework_data

# Explore structure of the dataset
str(data)
#or
names(data)

# Check structure of each
str(data$patients)
str(data$biomarker_data) 


# Pivot biomarker data to wide format 
#(so each patient is one row, and each protein becomes its own column)
protein_wide <- data$biomarker_data %>%
  filter(time == 0) %>%  # Only use baseline measurement
  select(patient_id, biomarker, value) %>%
  pivot_wider(names_from = biomarker, values_from = value)

# How many missing values (NA)?
colSums(is.na(protein_wide))
#Or
summary(protein_wide)

# Merge patient data with baseline Protein Expression
full_data <- data$patients %>%
  filter(event_type == "death") %>%  # Only patients who died
  select(patient_id, time_to_event, age, sex, hospital) %>%
  left_join(protein_wide, by = "patient_id")

#To Build Protein Risk Score
#First We Identify protein columns
protein_vars <- setdiff(
               names(full_data),
              c("patient_id", "time_to_event", "age", "sex", "hospital"))

#Then Fit linear model and compute risk score
protein_lm <- lm(time_to_event ~ ., data = full_data[, c("time_to_event", protein_vars)])
#Finally we can get predicted time_to_event values — the Protein Risk Score:
full_data$ProteinRiskScore <- predict(protein_lm, newdata = full_data)

#Fit Linear Mixed-Effects Model
lme_model <- lmer(time_to_event ~ age + sex + ProteinRiskScore + (1 | hospital), data = full_data)
# View summary of model
summary(lme_model)

#---------------------------------------------
####################

library(survival)
library(dplyr)
library(survminer)

#Here we use all patients (not just those who died).
patients <- data$patients


# Define survival time and event status
surv_data <- patients %>%
  mutate(event = ifelse(event_type == "death", 1, 0))

#Visual Exploration Before Categorization
cox_age <- coxph(Surv(time_to_event, event) ~ age, data = surv_data)
ggcoxfunctional(fit = cox_age, data = surv_data)

#Accordingly The Best Move Is To Devide like
surv_data <- surv_data %>%
  mutate(age_group = cut(age, breaks = c(0, 40, 60, 100), right = FALSE,
                         labels = c("<40", "40–59", "60+")))

#Check If Sample size is enough in each group
table(surv_data$age_group, surv_data$event)

# Create survival object
surv_obj <- Surv(surv_data$time_to_event, surv_data$event)

# Fit KM curve by age group
km_fit <- survfit(surv_obj ~ age_group, data = surv_data)

# Kaplan-Meier Plot by Age Group
ggsurvplot(km_fit, data = surv_data, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Survival Curves by Age Group",
           xlab = "Time", ylab = "Survival Probability")


# Log-Rank Test for Age Effect
survdiff(surv_obj ~ age_group, data = surv_data)

#Multivariate Cox Model (Sex, Age, Hospital)
cox1 <- coxph(surv_obj ~ age + sex + hospital, data = surv_data)
summary(cox1)
# Tests proportional hazards assumption
cox.zph(cox1)

#Time-Dependent Effect Of Age Cox Model 
cox2 <- coxph(surv_obj ~ sex + hospital + tt(age), 
              data = surv_data,
              tt = function(x, t, ...) x * log(t + 1))
summary(cox2)





