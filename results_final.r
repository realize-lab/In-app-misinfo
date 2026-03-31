# imports
library(tidyverse)
library(stringr)
library(broom)
library(lme4)
library(afex)
library(texreg)
library(rstatix)
library(lmerTest)
library(caret)

# load data
FILE_PATH <- "revised_table_increased_size.csv"
# make sure the file exists
if (!file.exists(FILE_PATH)) {
  stop("File not found")
}

data <- read.csv(FILE_PATH, stringsAsFactors = FALSE)

SAVE_FLAG <- TRUE

# function to save table to custom csv file
SAVE_PATH <- "results/"
FIG_PATH <- "violinfigs/"
# Create directory if it doesn't exist
if (!dir.exists(SAVE_PATH)) {
  dir.create(SAVE_PATH)
}

# Create directory if it doesn't exist
if (!dir.exists(FIG_PATH)) {
  dir.create(FIG_PATH)
}

save_table <- function(table, name) {
    if (SAVE_FLAG == FALSE) {
        return()
    }
    write.csv(table, paste0(SAVE_PATH, name, ".csv"))
}

save_fig <- function(fig, name) {
    if (SAVE_FLAG == FALSE) {
        return()
    }
    ggsave(paste0(FIG_PATH, name, ".png"), fig, width=10, height=6)
    ggsave(paste0(FIG_PATH, name, ".pdf"), fig, width=10, height=6)
}

# function to ensure model has converged
check_model_convergence <- function(model) {
    if (is.null(model)) {
        # stop if model is null
        stop("Model is null")
    }
    if (any(grepl("failed to converge", model@optinfo$conv$lme4))) {
        # exit if model has not converged
        stop("Model has not converged")
        return(FALSE)
    }
    return(TRUE)
}


# CSV columns present in the data
column <- list('Intervention_ID', "Tag", 'Misinfo', "Misinfo_Answer", 'Trust',  "Child", "BreastCancer" , 'Correct' , "Correct+Unsure", "Adherence" , "Adherence+Unsure ", "Unsure", 'Participant_ID')

# Cleaning data - setting values to respective categorical values
data$Intervention_ID<-as.factor(data$Intervention_ID)
data$Trust<-as.factor(data$Trust)
data$Misinfo<-as.factor(data$Misinfo)
data$Misinfo_Answer<-as.factor(data$Misinfo_Answer)
data$Tag<-as.factor(data$Tag)
data$Familiarity<-as.factor(data$Familiarity)
data$Participant_ID<-as.factor(data$Participant_ID)
data$BreastCancer<-as.factor(data$BreastCancer)
data$Race<-as.factor(data$Race)
data$Gender<-as.factor(data$Gender)
data$Child<-as.factor(data$Child)
data$General<-as.factor(data$General)

interventions <- c("Baseline", "In-line", "Pop-up", "RAG", "Priming")

get_results_table <- function(summary, observed, fixed) {
  results_table <- exp(cbind(OR = summary$coefficients[,1], 
                           LB = summary$coefficients[,1] - summary$coefficients[,2]*qnorm(0.975), 
                           UB = summary$coefficients[,1] + summary$coefficients[,2]*qnorm(0.975)))

  # add column for observed and fixed variables
    results_table <- cbind(fixed, observed, results_table)

  # create table without first row - intercept
  results_table_no_first_row <- results_table[-1,]
  return(results_table_no_first_row)
}

df_path <- "dfs/"
# Create directory if it doesn't exist
if (!dir.exists(df_path)) {
  dir.create(df_path)
}

save_df <- function(df, name) {
  write.csv(df, paste0(df_path, name, ".csv"))
}

# template data -> modeling variables -> get results -> save results to file
#****************************
## Result 1: Control group only, can participants accurately identify misinformation?
### Data - control group only i.e. Baseline
control_df <- data %>% filter(Intervention_ID == "Baseline")


### Modeling variables - Participant_ID, Tag: random, Misinfo: fixed, Correct/Trust/Adherence/Unsure: observed
control_model_acc <- glmer(Correct ~ Misinfo + (1|Participant_ID), data = control_df, family = binomial)
control_model_acc_unsure <- glmer(Correct+Unsure ~ Misinfo + (1|Participant_ID), data = control_df, family = binomial)
control_model_trust <- glmer(Trust ~ Misinfo + (1|Participant_ID), data = control_df, family = binomial)
control_model_adherence <- glmer(Adherence ~ Misinfo + (1|Participant_ID), data = control_df, family = binomial)
control_model_adherence_unsure <- glmer(Adherence+Unsure ~ Misinfo + (1|Participant_ID), data = control_df, family = binomial)
control_model_unsure <- glmer(Unsure ~ Misinfo + (1|Participant_ID), data = control_df, family = binomial)


### Check model convergence
check_model_convergence(control_model_acc)
check_model_convergence(control_model_acc_unsure)
check_model_convergence(control_model_trust)
check_model_convergence(control_model_adherence)
check_model_convergence(control_model_adherence_unsure)
check_model_convergence(control_model_unsure)
check_model_convergence(control_model_harmful)

### Get results
summary_control_acc <- summary(control_model_acc)
summary_control_acc_unsure <- summary(control_model_acc_unsure)
summary_control_trust <- summary(control_model_trust)
summary_control_adherence <- summary(control_model_adherence)
summary_control_adherence_unsure <- summary(control_model_adherence_unsure)
summary_control_unsure <- summary(control_model_unsure)
summary_control_harmful <- summary(control_model_harmful)


results_table_acc <- get_results_table(summary_control_acc, "Correct", "Misinformation")
results_table_acc_unsure <- get_results_table(summary_control_acc_unsure, "Correct+Unsure", "Misinformation")
results_table_trust <- get_results_table(summary_control_trust, "Trust", "Misinformation")
results_table_adherence <- get_results_table(summary_control_adherence, "Adherence", "Misinformation")
results_table_adherence_unsure <- get_results_table(summary_control_adherence_unsure, "Adherence+Unsure", "Misinformation")
results_table_unsure <- get_results_table(summary_control_unsure, "Unsure", "Misinformation")
results_table_harmful <- get_results_table(summary_control_harmful, "Misinfo_Answer", "Misinformation")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table <- rbind(results_table_acc, results_table_trust, results_table_adherence, results_table_unsure)
results_table_with_metrics_unsure <- rbind(results_table_acc, results_table_acc_unsure, results_table_trust, results_table_adherence, results_table_adherence_unsure, results_table_unsure)
# rownames(results_table) <- c("Correct", "Trust", "Adherence", "Unsure")

# Print results
print("Control Group Discernment")
print(results_table_with_metrics_unsure)

### Save results to file
save_table(results_table, "control_group_misinfo_discernment")
save_df(control_df, "control_group_misinfo_discernment")
save_table(results_table_with_metrics_unsure, "control_group_misinfo_discernment_with_metrics_unsure")
save_df(control_df, "control_group_misinfo_discernment_with_metrics_unsure")

## Result 2: Effect of intervention on metrics given misinformation

### Result 2.1: Misinfo true questions only
### Data - all interventions, data where Misinfo is true
misinfo_df <- data %>% filter(Misinfo == "True")

### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
misinfo_model_acc <- glmer(Correct ~ Intervention_ID + (1| Participant_ID) + (1| Tag), 
                    data = misinfo_df, 
                    family = binomial)

misinfo_model_acc_unsure <- glmer(Correct+Unsure ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = misinfo_df, 
                    family = binomial)

misinfo_model_trust <- glmer(Trust ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = misinfo_df, 
                    family = binomial)

misinfo_model_adherence <- glmer(Adherence ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = misinfo_df, 
                    family = binomial)

misinfo_model_adherence_unsure <- glmer(Adherence+Unsure ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = misinfo_df, 
                    family = binomial)

misinfo_model_unsure <- glmer(Unsure ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = misinfo_df, 
                    family = binomial)
misinfo_model_harmful <- glmer(Misinfo_Answer ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = misinfo_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(misinfo_model_acc)
check_model_convergence(misinfo_model_acc_unsure)
check_model_convergence(misinfo_model_trust)
check_model_convergence(misinfo_model_adherence)
check_model_convergence(misinfo_model_adherence_unsure)
check_model_convergence(misinfo_model_unsure)
check_model_convergence(misinfo_model_harmful)


### Get results
summary_misinfo_acc <- summary(misinfo_model_acc)
summary_misinfo_acc_unsure <- summary(misinfo_model_acc_unsure)
summary_misinfo_trust <- summary(misinfo_model_trust)
summary_misinfo_adherence <- summary(misinfo_model_adherence)
summary_misinfo_adherence_unsure <- summary(misinfo_model_adherence_unsure)
summary_misinfo_unsure <- summary(misinfo_model_unsure)
summary_misinfo_harmful <- summary(misinfo_model_harmful)

print("Misinformation Discernment")


results_table_misinfo_acc <- get_results_table(summary_misinfo_acc, "Correct", "Intervention")
results_table_misinfo_acc_unsure <- get_results_table(summary_misinfo_acc_unsure, "Correct+Unsure", "Intervention")
results_table_misinfo_trust <- get_results_table(summary_misinfo_trust, "Trust", "Intervention")
results_table_misinfo_adherence <- get_results_table(summary_misinfo_adherence, "Adherence", "Intervention")
results_table_misinfo_adherence_unsure <- get_results_table(summary_misinfo_adherence_unsure, "Adherence+Unsure", "Intervention")
results_table_misinfo_unsure <- get_results_table(summary_misinfo_unsure, "Unsure", "Intervention")
results_table_misinfo_harmful <- get_results_table(summary_misinfo_harmful, "Misinfo_Answer", "Intervention")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table <- rbind(results_table_misinfo_acc, results_table_misinfo_trust, results_table_misinfo_adherence, results_table_misinfo_unsure)
print(results_table)

result_table_with_metrics_unsure <- rbind(results_table_misinfo_acc, results_table_misinfo_acc_unsure, results_table_misinfo_trust, results_table_misinfo_adherence, results_table_misinfo_adherence_unsure, results_table_misinfo_unsure)
print(result_table_with_metrics_unsure)

#### Save results to file
save_table(results_table_misinfo_acc, "misinfo_interventions_correct")
save_table(results_table_misinfo_acc_unsure, "misinfo_interventions_correct_unsure")
save_table(results_table_misinfo_trust, "misinfo_interventions_trust")
save_table(results_table_misinfo_adherence, "misinfo_interventions_adherence")
save_table(results_table_misinfo_adherence_unsure, "misinfo_interventions_adherence_unsure")
save_table(results_table_misinfo_unsure, "misinfo_interventions_unsure")
##### combine results saving
save_table(results_table, "misinfo_interventions")
save_df(misinfo_df, "misinfo_interventions")
save_table(result_table_with_metrics_unsure, "misinfo_interventions_with_metrics_unsure")
save_df(misinfo_df, "misinfo_interventions_with_metrics_unsure")

## Result 2.2: Misinfo false questions only
### Data - all interventions, data where Misinfo is false
misinfo_df <- data %>% filter(Misinfo == "False")

### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
no_misinfo_model_acc <- glmer(Correct ~ Intervention_ID + (1| Participant_ID) + (1| Tag), 
                    data = misinfo_df, 
                    family = binomial)

no_misinfo_model_acc_unsure <- glmer(Correct+Unsure ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = misinfo_df, 
                    family = binomial)

no_misinfo_model_trust <- glmer(Trust ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = misinfo_df, 
                    family = binomial)

no_misinfo_model_adherence <- glmer(Adherence ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = misinfo_df, 
                    family = binomial)

no_misinfo_model_adherence_unsure <- glmer(Adherence+Unsure ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = misinfo_df, 
                    family = binomial)

no_misinfo_model_unsure <- glmer(Unsure ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = misinfo_df, 
                    family = binomial)

no_misinfo_model_harmful <- glmer(Misinfo_Answer ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = misinfo_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(no_misinfo_model_acc)
check_model_convergence(no_misinfo_model_acc_unsure)
check_model_convergence(no_misinfo_model_trust)
check_model_convergence(no_misinfo_model_adherence)
check_model_convergence(no_misinfo_model_adherence_unsure)
check_model_convergence(no_misinfo_model_unsure)
check_model_convergence(no_misinfo_model_harmful)

### Get results
summary_no_misinfo_acc <- summary(no_misinfo_model_acc)
summary_no_misinfo_acc_unsure <- summary(no_misinfo_model_acc_unsure)
summary_no_misinfo_trust <- summary(no_misinfo_model_trust)
summary_no_misinfo_adherence <- summary(no_misinfo_model_adherence)
summary_no_misinfo_adherence_unsure <- summary(no_misinfo_model_adherence_unsure)
summary_no_misinfo_unsure <- summary(no_misinfo_model_unsure)
summary_no_misinfo_harmful <- summary(no_misinfo_model_harmful)

print("No Misinformation Discernment")

results_table_no_misinfo_acc <- get_results_table(summary_no_misinfo_acc, "Correct", "Intervention")
results_table_no_misinfo_acc_unsure <- get_results_table(summary_no_misinfo_acc_unsure, "Correct+Unsure", "Intervention")
results_table_no_misinfo_trust <- get_results_table(summary_no_misinfo_trust, "Trust", "Intervention")
results_table_no_misinfo_adherence <- get_results_table(summary_no_misinfo_adherence, "Adherence", "Intervention")
results_table_no_misinfo_adherence_unsure <- get_results_table(summary_no_misinfo_adherence_unsure, "Adherence+Unsure", "Intervention")
results_table_no_misinfo_unsure <- get_results_table(summary_no_misinfo_unsure, "Unsure", "Intervention")
results_table_no_misinfo_harmful <- get_results_table(summary_no_misinfo_harmful, "Misinfo_Answer", "Intervention")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table_no_misinfo <- rbind(results_table_no_misinfo_acc, results_table_no_misinfo_trust, results_table_no_misinfo_adherence, results_table_no_misinfo_unsure)
print("Metrics only")
print(results_table_no_misinfo)

result_table_no_misinfo_with_metrics_unsure <- rbind(results_table_no_misinfo_acc, results_table_no_misinfo_acc_unsure, results_table_no_misinfo_trust, results_table_no_misinfo_adherence, results_table_no_misinfo_adherence_unsure, results_table_no_misinfo_unsure)
print("Metrics + unsure")
print(result_table_no_misinfo_with_metrics_unsure)

#### Save results to file
save_table(results_table_no_misinfo_acc, "no_misinfo_interventions_correct")
save_table(results_table_no_misinfo_acc_unsure, "no_misinfo_interventions_correct_unsure")
save_table(results_table_no_misinfo_trust, "no_misinfo_interventions_trust")
save_table(results_table_no_misinfo_adherence, "no_misinfo_interventions_adherence")
save_table(results_table_no_misinfo_adherence_unsure, "no_misinfo_interventions_adherence_unsure")
save_table(results_table_no_misinfo_unsure, "no_misinfo_interventions_unsure")
##### combine results saving
save_table(results_table_no_misinfo, "no_misinfo_interventions")
save_df(misinfo_df, "no_misinfo_interventions")
save_table(result_table_no_misinfo_with_metrics_unsure, "no_misinfo_interventions_with_metrics_unsure")
save_df(misinfo_df, "no_misinfo_interventions_with_metrics_unsure")

## Result 2.3: Interaction effect of intervention and misinformation on metrics

### Data - all data
all_df <- data

### Modeling variables - Participant_ID, Tag: random, Intervention_ID, Misinfo and their interaction: fixed, Correct/Trust/Adherence/Unsure: observed
interaction_model_acc <- glmer(Correct ~ Intervention_ID + (1| Participant_ID) + (1| Tag), 
                    data = all_df, 
                    family = binomial)

interaction_model_acc_unsure <- glmer(Correct+Unsure ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = all_df, 
                    family = binomial)

interaction_model_trust <- glmer(Trust ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = all_df, 
                    family = binomial)

interaction_model_adherence <- glmer(Adherence ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = all_df, 
                    family = binomial)

interaction_model_adherence_unsure <- glmer(Adherence+Unsure ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = all_df, 
                    family = binomial)

interaction_model_unsure <- glmer(Unsure ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = all_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(interaction_model_acc)
check_model_convergence(interaction_model_acc_unsure)
check_model_convergence(interaction_model_trust)
check_model_convergence(interaction_model_adherence)
check_model_convergence(interaction_model_adherence_unsure)
check_model_convergence(interaction_model_unsure)

### Get results
summary_interaction_acc <- summary(interaction_model_acc)
summary_interaction_acc_unsure <- summary(interaction_model_acc_unsure)
summary_interaction_trust <- summary(interaction_model_trust)
summary_interaction_adherence <- summary(interaction_model_adherence)
summary_interaction_adherence_unsure <- summary(interaction_model_adherence_unsure)
summary_interaction_unsure <- summary(interaction_model_unsure)

print("Interaction between Intervention and Misinformation")

# Extract results from interaction models and save
results_table_interaction_acc <- get_results_table(summary_interaction_acc, "Correct", "Intervention")
results_table_interaction_acc_unsure <- get_results_table(summary_interaction_acc_unsure, "Correct+Unsure", "Intervention")
results_table_interaction_trust <- get_results_table(summary_interaction_trust, "Trust", "Intervention")
results_table_interaction_adherence <- get_results_table(summary_interaction_adherence, "Adherence", "Intervention")
results_table_interaction_adherence_unsure <- get_results_table(summary_interaction_adherence_unsure, "Adherence+Unsure", "Intervention")
results_table_interaction_unsure <- get_results_table(summary_interaction_unsure, "Unsure", "Intervention")

# Print and save the interaction results
save_table(results_table_interaction_acc, "interaction_interventions_misinfo_correct")
save_table(results_table_interaction_acc_unsure, "interaction_interventions_misinfo_correct_unsure")
save_table(results_table_interaction_trust, "interaction_interventions_misinfo_trust")
save_table(results_table_interaction_adherence, "interaction_interventions_misinfo_adherence")
save_table(results_table_interaction_adherence_unsure, "interaction_interventions_misinfo_adherence_unsure")
save_table(results_table_interaction_unsure, "interaction_interventions_misinfo_unsure")

# Combined results 
result_table_interaction_combined <- rbind(results_table_interaction_acc, results_table_interaction_acc_unsure, 
                                        results_table_interaction_trust, results_table_interaction_adherence,
                                        results_table_interaction_adherence_unsure, results_table_interaction_unsure)
save_table(result_table_interaction_combined, "interaction_interventions_misinfo_combined")
save_df(all_df, "interaction_interventions_misinfo")

## Result 3: Effect of intervention on strict correctness

### Result 3.1: Misinfo true questions only
### Data - all interventions, data where Misinfo is true
strict_misinfo_df <- data %>% filter(Misinfo == "True", Unsure == 0)

### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
strict_misinfo_model_acc <- glmer(Correct ~ Intervention_ID + (1| Participant_ID) + (1| Tag), 
                    data = strict_misinfo_df, 
                    family = binomial)
strict_misinfo_model_trust <- glmer(Trust ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = strict_misinfo_df, 
                    family = binomial)
strict_misinfo_model_adherence <- glmer(Adherence ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = strict_misinfo_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(strict_misinfo_model_acc)
check_model_convergence(strict_misinfo_model_trust)
check_model_convergence(strict_misinfo_model_adherence)


### Get results
summary_strict_misinfo_acc <- summary(strict_misinfo_model_acc)
summary_strict_misinfo_trust <- summary(strict_misinfo_model_trust)
summary_strict_misinfo_adherence <- summary(strict_misinfo_model_adherence)

print("Strict Misinformation Discernment")

results_table_strict_misinfo_acc <- get_results_table(summary_strict_misinfo_acc, "Correct", "Intervention")
results_table_strict_misinfo_trust <- get_results_table(summary_strict_misinfo_trust, "Trust", "Intervention")
results_table_strict_misinfo_adherence <- get_results_table(summary_strict_misinfo_adherence, "Adherence", "Intervention")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table_strict_misinfo <- rbind(results_table_strict_misinfo_acc, results_table_strict_misinfo_trust, results_table_strict_misinfo_adherence)
print(results_table_strict_misinfo)

#### Save results to file
save_table(results_table_strict_misinfo_acc, "strict_misinfo_interventions_correct")
save_table(results_table_strict_misinfo_trust, "strict_misinfo_interventions_trust")
save_table(results_table_strict_misinfo_adherence, "strict_misinfo_interventions_adherence")
##### combine results saving
save_table(results_table_strict_misinfo, "strict_misinfo_interventions")
save_df(strict_misinfo_df, "strict_misinfo_interventions")

## Result 3.2: Misinfo false questions only
### Data - all interventions, data where Misinfo is false
strict_no_misinfo_df <- data %>% filter(Misinfo == "False", Unsure == 0)

### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
strict_no_misinfo_model_acc <- glmer(Correct ~ Intervention_ID + (1| Participant_ID) + (1| Tag), 
                    data = strict_no_misinfo_df, 
                    family = binomial)
strict_no_misinfo_model_trust <- glmer(Trust ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = strict_no_misinfo_df, 
                    family = binomial)
strict_no_misinfo_model_adherence <- glmer(Adherence ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = strict_no_misinfo_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(strict_no_misinfo_model_acc)
check_model_convergence(strict_no_misinfo_model_trust)
check_model_convergence(strict_no_misinfo_model_adherence)

### Get results
summary_strict_no_misinfo_acc <- summary(strict_no_misinfo_model_acc)
summary_strict_no_misinfo_trust <- summary(strict_no_misinfo_model_trust)
summary_strict_no_misinfo_adherence <- summary(strict_no_misinfo_model_adherence)

print("Strict No Misinformation Discernment")

results_table_strict_no_misinfo_acc <- get_results_table(summary_strict_no_misinfo_acc, "Correct", "Intervention")
results_table_strict_no_misinfo_trust <- get_results_table(summary_strict_no_misinfo_trust, "Trust", "Intervention")
results_table_strict_no_misinfo_adherence <- get_results_table(summary_strict_no_misinfo_adherence, "Adherence", "Intervention")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table_strict_no_misinfo <- rbind(results_table_strict_no_misinfo_acc, results_table_strict_no_misinfo_trust, results_table_strict_no_misinfo_adherence)
print(results_table_strict_no_misinfo)

#### Save results to file
save_table(results_table_strict_no_misinfo_acc, "strict_no_misinfo_interventions_correct")
save_table(results_table_strict_no_misinfo_trust, "strict_no_misinfo_interventions_trust")
save_table(results_table_strict_no_misinfo_adherence, "strict_no_misinfo_interventions_adherence")

##### combine results saving
save_table(results_table_strict_no_misinfo, "strict_no_misinfo_interventions")
save_df(strict_no_misinfo_df, "strict_no_misinfo_interventions")

## Result 3.3: Interaction effect of intervention and misinformation on strict correctness

### Data - all data
strict_all_df <- data %>% filter(Unsure == 0)

### Modeling variables - Participant_ID, Tag: random, Intervention_ID, Misinfo and their interaction: fixed, Correct/Trust/Adherence/Unsure: observed
strict_interaction_model_acc <- glmer(Correct ~ Intervention_ID + (1| Participant_ID) + (1| Tag), 
                    data = strict_all_df, 
                    family = binomial)
strict_interaction_model_trust <- glmer(Trust ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = strict_all_df, 
                    family = binomial)
strict_interaction_model_adherence <- glmer(Adherence ~ Intervention_ID + (1| Participant_ID) + (1| Tag),
                    data = strict_all_df, 
                    family = binomial)

# Check model convergence
check_model_convergence(strict_interaction_model_acc)
check_model_convergence(strict_interaction_model_trust)
check_model_convergence(strict_interaction_model_adherence)

### Get results
summary_strict_interaction_acc <- summary(strict_interaction_model_acc)
summary_strict_interaction_trust <- summary(strict_interaction_model_trust)
summary_strict_interaction_adherence <- summary(strict_interaction_model_adherence)

print("Strict Interaction between Intervention and Misinformation")

# Extract results from interaction models and save
results_table_strict_interaction_acc <- get_results_table(summary_strict_interaction_acc, "Correct", "Intervention")
results_table_strict_interaction_trust <- get_results_table(summary_strict_interaction_trust, "Trust", "Intervention")
results_table_strict_interaction_adherence <- get_results_table(summary_strict_interaction_adherence, "Adherence", "Intervention")

# Print and save the interaction results
save_table(results_table_strict_interaction_acc, "strict_interaction_interventions_misinfo_correct")
save_table(results_table_strict_interaction_trust, "strict_interaction_interventions_misinfo_trust")
save_table(results_table_strict_interaction_adherence, "strict_interaction_interventions_misinfo_adherence")

# Combined results
result_table_strict_interaction_combined <- rbind(results_table_strict_interaction_acc, results_table_strict_interaction_trust, 
                                        results_table_strict_interaction_adherence)
save_table(result_table_strict_interaction_combined, "strict_interaction_interventions_misinfo_combined")
save_df(strict_all_df, "strict_interaction_interventions_misinfo")

## Result 4: Effect of intervention on the breastcancer demography
print("Starting Breast Cancer Analysis")
### Result 4.1: Misinfo true questions only
### Data - all interventions, data where Misinfo is true
breastcancer_df <- data %>% filter(BreastCancer == "True", Tag == "hrt" | Tag == "hrt_2" | Tag == "mammogram", Misinfo == "True") # BreastCancer related
print(nrow(breastcancer_df))

### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
breastcancer_model_acc <- glmer(Correct ~ Intervention_ID + (1| Participant_ID), 
                    data = breastcancer_df, 
                    family = binomial)
breastcancer_model_acc_unsure <- glmer(Correct+Unsure ~ Intervention_ID + (1| Participant_ID),
                    data = breastcancer_df, 
                    family = binomial)
breastcancer_model_trust <- glmer(Trust ~ Intervention_ID + (1| Participant_ID),
                    data = breastcancer_df, 
                    family = binomial)
breastcancer_model_adherence <- glmer(Adherence ~ Intervention_ID + (1| Participant_ID),
                    data = breastcancer_df, 
                    family = binomial)
breastcancer_model_adherence_unsure <- glmer(Adherence+Unsure ~ Intervention_ID + (1| Participant_ID),
                    data = breastcancer_df, 
                    family = binomial)
breastcancer_model_unsure <- glmer(Unsure ~ Intervention_ID + (1| Participant_ID),
                    data = breastcancer_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(breastcancer_model_acc)
check_model_convergence(breastcancer_model_acc_unsure)
check_model_convergence(breastcancer_model_trust)
check_model_convergence(breastcancer_model_adherence)
check_model_convergence(breastcancer_model_adherence_unsure)
check_model_convergence(breastcancer_model_unsure)

### Get results
summary_breastcancer_acc <- summary(breastcancer_model_acc)
summary_breastcancer_acc_unsure <- summary(breastcancer_model_acc_unsure)
summary_breastcancer_trust <- summary(breastcancer_model_trust)
summary_breastcancer_adherence <- summary(breastcancer_model_adherence)
summary_breastcancer_adherence_unsure <- summary(breastcancer_model_adherence_unsure)
summary_breastcancer_unsure <- summary(breastcancer_model_unsure)

print("Breast Cancer Discernment")

results_table_breastcancer_acc <- get_results_table(summary_breastcancer_acc, "Correct", "Intervention")
results_table_breastcancer_acc_unsure <- get_results_table(summary_breastcancer_acc_unsure, "Correct+Unsure", "Intervention")
results_table_breastcancer_trust <- get_results_table(summary_breastcancer_trust, "Trust", "Intervention")
results_table_breastcancer_adherence <- get_results_table(summary_breastcancer_adherence, "Adherence", "Intervention")
results_table_breastcancer_adherence_unsure <- get_results_table(summary_breastcancer_adherence_unsure, "Adherence+Unsure", "Intervention")
results_table_breastcancer_unsure <- get_results_table(summary_breastcancer_unsure, "Unsure", "Intervention")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table_breastcancer <- rbind(results_table_breastcancer_acc, results_table_breastcancer_trust, results_table_breastcancer_adherence, results_table_breastcancer_unsure)
results_table_breastcancer_with_metrics_unsure <- rbind(results_table_breastcancer_acc, results_table_breastcancer_acc_unsure, results_table_breastcancer_trust, results_table_breastcancer_adherence, results_table_breastcancer_adherence_unsure, results_table_breastcancer_unsure)
print(results_table_breastcancer)
print(results_table_breastcancer_with_metrics_unsure)

#### Save results to file
save_table(results_table_breastcancer_acc, "breastcancer_misinfo_interventions_correct")
save_table(results_table_breastcancer_acc_unsure, "breastcancer_misinfo_interventions_correct_unsure")
save_table(results_table_breastcancer_trust, "breastcancer_misinfo_interventions_trust")
save_table(results_table_breastcancer_adherence, "breastcancer_misinfo_interventions_adherence")
save_table(results_table_breastcancer_adherence_unsure, "breastcancer_misinfo_interventions_adherence_unsure")
save_table(results_table_breastcancer_unsure, "breastcancer_misinfo_interventions_unsure")
##### combine results saving
save_table(results_table_breastcancer, "breastcancer_misinfo_interventions")
save_df(breastcancer_df, "breastcancer_misinfo_interventions")
save_table(results_table_breastcancer_with_metrics_unsure, "breastcancer_misinfo_interventions_with_metrics_unsure")
save_df(breastcancer_df, "breastcancer_misinfo_interventions_with_metrics_unsure")


# ## Result 4.2: Misinfo false questions only
# ### Data - all interventions, data where Misinfo is false
# print("Starting Breast Cancer Analysis - No Misinformation")
# true_breastcancer_df <- data %>% filter(BreastCancer == "True", Tag == "hrt" | Tag == "hrt_2" | Tag == "mammogram", Misinfo=="False") # BreastCancer related
# print(nrow(true_breastcancer_df))
# # print unique interventions, correct, unsure, trust, adherence, tag
# print(unique(true_breastcancer_df$Correct))
# print(unique(true_breastcancer_df$Unsure))
# print(unique(true_breastcancer_df$Trust))
# print(unique(true_breastcancer_df$Adherence))
# print(unique(true_breastcancer_df$Tag))
# print(unique(true_breastcancer_df$Participant_ID))
# print(unique(true_breastcancer_df$Intervention_ID))
# print("Starting Breast Cancer Analysis - No Misinformation data done")

# ### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
# true_breastcancer_model_acc <- glmer(Correct ~ Intervention_ID + (1| Participant_ID), 
#                     data = true_breastcancer_df, 
#                     family = binomial)
# true_breastcancer_model_acc_unsure <- glmer(Correct+Unsure ~ Intervention_ID + (1| Participant_ID),
#                     data = true_breastcancer_df, 
#                     family = binomial)
# true_breastcancer_model_trust <- glmer(Trust ~ Intervention_ID + (1| Participant_ID),
#                     data = true_breastcancer_df, 
#                     family = binomial)
# true_breastcancer_model_adherence <- glmer(Adherence ~ Intervention_ID + (1| Participant_ID),
#                     data = true_breastcancer_df, 
#                     family = binomial)
# true_breastcancer_model_adherence_unsure <- glmer(Adherence+Unsure ~ Intervention_ID + (1| Participant_ID),
#                     data = true_breastcancer_df, 
#                     family = binomial)
# true_breastcancer_model_unsure <- glmer(Unsure ~ Intervention_ID + (1| Participant_ID),
#                     data = true_breastcancer_df, 
#                     family = binomial)

# print("Starting Breast Cancer Analysis - No Misinformation models done")
# ### Check model convergence
# check_model_convergence(true_breastcancer_model_acc)
# check_model_convergence(true_breastcancer_model_acc_unsure)
# check_model_convergence(true_breastcancer_model_trust)
# check_model_convergence(true_breastcancer_model_adherence)
# check_model_convergence(true_breastcancer_model_adherence_unsure)
# check_model_convergence(true_breastcancer_model_unsure)

# ### Get results
# summary_true_breastcancer_acc <- summary(true_breastcancer_model_acc)
# summary_true_breastcancer_acc_unsure <- summary(true_breastcancer_model_acc_unsure)
# summary_true_breastcancer_trust <- summary(true_breastcancer_model_trust)
# summary_true_breastcancer_adherence <- summary(true_breastcancer_model_adherence)
# summary_true_breastcancer_adherence_unsure <- summary(true_breastcancer_model_adherence_unsure)
# summary_true_breastcancer_unsure <- summary(true_breastcancer_model_unsure)

# print("True Breast Cancer Discernment")

# results_table_true_breastcancer_acc <- get_results_table(summary_true_breastcancer_acc, "Correct", "Intervention")
# results_table_true_breastcancer_acc_unsure <- get_results_table(summary_true_breastcancer_acc_unsure, "Correct+Unsure", "Intervention")
# results_table_true_breastcancer_trust <- get_results_table(summary_true_breastcancer_trust, "Trust", "Intervention")
# results_table_true_breastcancer_adherence <- get_results_table(summary_true_breastcancer_adherence, "Adherence", "Intervention")
# results_table_true_breastcancer_adherence_unsure <- get_results_table(summary_true_breastcancer_adherence_unsure, "Adherence+Unsure", "Intervention")
# results_table_true_breastcancer_unsure <- get_results_table(summary_true_breastcancer_unsure, "Unsure", "Intervention")

# # Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB

# results_table_true_breastcancer <- rbind(results_table_true_breastcancer_acc, results_table_true_breastcancer_trust, results_table_true_breastcancer_adherence, results_table_true_breastcancer_unsure)
# results_table_true_breastcancer_with_metrics_unsure <- rbind(results_table_true_breastcancer_acc, results_table_true_breastcancer_acc_unsure, results_table_true_breastcancer_trust, results_table_true_breastcancer_adherence, results_table_true_breastcancer_adherence_unsure, results_table_true_breastcancer_unsure)
# print(results_table_true_breastcancer)
# print(results_table_true_breastcancer_with_metrics_unsure)

# #### Save results to file
# save_table(results_table_true_breastcancer_acc, "breastcancer_no_misinfo_interventions_correct")
# save_table(results_table_true_breastcancer_acc_unsure, "breastcancer_no_misinfo_interventions_correct_unsure")
# save_table(results_table_true_breastcancer_trust, "breastcancer_no_misinfo_interventions_trust")
# save_table(results_table_true_breastcancer_adherence, "breastcancer_no_misinfo_interventions_adherence")
# save_table(results_table_true_breastcancer_adherence_unsure, "breastcancer_no_misinfo_interventions_adherence_unsure")
# save_table(results_table_true_breastcancer_unsure, "breastcancer_no_misinfo_interventions_unsure")

# ##### combine results saving
# save_table(results_table_true_breastcancer, "breastcancer_no_misinfo_interventions")
# save_df(true_breastcancer_df, "breastcancer_no_misinfo_interventions")
# save_table(results_table_true_breastcancer_with_metrics_unsure, "breastcancer_no_misinfo_interventions_with_metrics_unsure")
# save_df(true_breastcancer_df, "breastcancer_no_misinfo_interventions_with_metrics_unsure")


## Result 5: Effect of intervention on the breastcancer demography with strict correctness

### Result 5.1: Misinfo true questions only
### Data - all interventions, data where Misinfo is true
strict_breastcancer_df <- data %>% filter(BreastCancer == "True", Unsure == 0, Tag == "hrt" | Tag == "hrt_2" | Tag == "mammogram") # BreastCancer related

### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
strict_breastcancer_model_acc <- glmer(Correct ~ Intervention_ID + (1| Participant_ID), 
                    data = strict_breastcancer_df, 
                    family = binomial)
strict_breastcancer_model_trust <- glmer(Trust ~ Intervention_ID + (1| Participant_ID),
                    data = strict_breastcancer_df, 
                    family = binomial)
strict_breastcancer_model_adherence <- glmer(Adherence ~ Intervention_ID + (1| Participant_ID),
                    data = strict_breastcancer_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(strict_breastcancer_model_acc)
check_model_convergence(strict_breastcancer_model_trust)
check_model_convergence(strict_breastcancer_model_adherence)

### Get results
summary_strict_breastcancer_acc <- summary(strict_breastcancer_model_acc)
summary_strict_breastcancer_trust <- summary(strict_breastcancer_model_trust)
summary_strict_breastcancer_adherence <- summary(strict_breastcancer_model_adherence)

print("Strict Breast Cancer Discernment")

results_table_strict_breastcancer_acc <- get_results_table(summary_strict_breastcancer_acc, "Correct", "Intervention")
results_table_strict_breastcancer_trust <- get_results_table(summary_strict_breastcancer_trust, "Trust", "Intervention")
results_table_strict_breastcancer_adherence <- get_results_table(summary_strict_breastcancer_adherence, "Adherence", "Intervention")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table_strict_breastcancer <- rbind(results_table_strict_breastcancer_acc, results_table_strict_breastcancer_trust, results_table_strict_breastcancer_adherence)
print(results_table_strict_breastcancer)

#### Save results to file
save_table(results_table_strict_breastcancer_acc, "strict_breastcancer_misinfo_interventions_correct")
save_table(results_table_strict_breastcancer_trust, "strict_breastcancer_misinfo_interventions_trust")
save_table(results_table_strict_breastcancer_adherence, "strict_breastcancer_misinfo_interventions_adherence")

##### combine results saving
save_table(results_table_strict_breastcancer, "strict_breastcancer_misinfo_interventions")
save_df(strict_breastcancer_df, "strict_breastcancer_misinfo_interventions")

# ## Result 5.2: Misinfo false questions only
# ### Data - all interventions, data where Misinfo is false

# strict_true_breastcancer_df <- data %>% filter(BreastCancer == "True", Unsure == 0, Tag == "hrt" | Tag == "hrt_2" | Tag == "mammogram", Misinfo=="False") # BreastCancer related

# ### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
# strict_true_breastcancer_model_acc <- glmer(Correct ~ Intervention_ID + Tag + (1| Participant_ID), 
#                     data = strict_true_breastcancer_df, 
#                     family = binomial)
# strict_true_breastcancer_model_trust <- glmer(Trust ~ Intervention_ID + Tag + (1| Participant_ID),
#                     data = strict_true_breastcancer_df, 
#                     family = binomial)
# strict_true_breastcancer_model_adherence <- glmer(Adherence ~ Intervention_ID + Tag + (1| Participant_ID),
#                     data = strict_true_breastcancer_df, 
#                     family = binomial)

# ### Check model convergence
# check_model_convergence(strict_true_breastcancer_model_acc)
# check_model_convergence(strict_true_breastcancer_model_trust)
# check_model_convergence(strict_true_breastcancer_model_adherence)

# ### Get results
# summary_strict_true_breastcancer_acc <- summary(strict_true_breastcancer_model_acc)
# summary_strict_true_breastcancer_trust <- summary(strict_true_breastcancer_model_trust)
# summary_strict_true_breastcancer_adherence <- summary(strict_true_breastcancer_model_adherence)

# print("Strict True Breast Cancer Discernment")

# results_table_strict_true_breastcancer_acc <- get_results_table(summary_strict_true_breastcancer_acc, "Correct", "Intervention")
# results_table_strict_true_breastcancer_trust <- get_results_table(summary_strict_true_breastcancer_trust, "Trust", "Intervention")
# results_table_strict_true_breastcancer_adherence <- get_results_table(summary_strict_true_breastcancer_adherence, "Adherence", "Intervention")

# # Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
# results_table_strict_true_breastcancer <- rbind(results_table_strict_true_breastcancer_acc, results_table_strict_true_breastcancer_trust, results_table_strict_true_breastcancer_adherence)
# print(results_table_strict_true_breastcancer)

# #### Save results to file
# save_table(results_table_strict_true_breastcancer_acc, "strict_breastcancer_no_misinfo_interventions_correct")
# save_table(results_table_strict_true_breastcancer_trust, "strict_breastcancer_no_misinfo_interventions_trust")
# save_table(results_table_strict_true_breastcancer_adherence, "strict_breastcancer_no_misinfo_interventions_adherence")

# ##### combine results saving
# save_table(results_table_strict_true_breastcancer, "strict_breastcancer_no_misinfo_interventions")
# save_df(strict_true_breastcancer_df, "strict_breastcancer_no_misinfo_interventions")

## Result 6: Interaction effect of intervention and have children on metrics

### Result 6.1: Misinfo true questions only
### Data - all interventions, data where Misinfo is true
children_df <- data %>% filter(Child == "True", Tag == "breastfeed" | Tag == "cosleep" | Tag == "rsv" | Tag == "formula", Misinfo == "True") # Children related

### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
children_model_acc <- glmer(Correct ~ Intervention_ID + (1| Participant_ID), 
                    data = children_df, 
                    family = binomial)
children_model_acc_unsure <- glmer(Correct+Unsure ~ Intervention_ID + (1| Participant_ID),
                    data = children_df, 
                    family = binomial)
children_model_trust <- glmer(Trust ~ Intervention_ID + (1| Participant_ID),
                    data = children_df, 
                    family = binomial)
children_model_adherence <- glmer(Adherence ~ Intervention_ID + (1| Participant_ID),
                    data = children_df, 
                    family = binomial)
children_model_adherence_unsure <- glmer(Adherence+Unsure ~ Intervention_ID + (1| Participant_ID),
                    data = children_df, 
                    family = binomial)
children_model_unsure <- glmer(Unsure ~ Intervention_ID + (1| Participant_ID),
                    data = children_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(children_model_acc)
check_model_convergence(children_model_acc_unsure)
check_model_convergence(children_model_trust)
check_model_convergence(children_model_adherence)
check_model_convergence(children_model_adherence_unsure)
check_model_convergence(children_model_unsure)


### Get results
summary_children_acc <- summary(children_model_acc)
summary_children_acc_unsure <- summary(children_model_acc_unsure)
summary_children_trust <- summary(children_model_trust)
summary_children_adherence <- summary(children_model_adherence)
summary_children_adherence_unsure <- summary(children_model_adherence_unsure)
summary_children_unsure <- summary(children_model_unsure)

print("Children Discernment")

results_table_children_acc <- get_results_table(summary_children_acc, "Correct", "Intervention")
results_table_children_acc_unsure <- get_results_table(summary_children_acc_unsure, "Correct+Unsure", "Intervention")
results_table_children_trust <- get_results_table(summary_children_trust, "Trust", "Intervention")
results_table_children_adherence <- get_results_table(summary_children_adherence, "Adherence", "Intervention")
results_table_children_adherence_unsure <- get_results_table(summary_children_adherence_unsure, "Adherence+Unsure", "Intervention")
results_table_children_unsure <- get_results_table(summary_children_unsure, "Unsure", "Intervention")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table_children <- rbind(results_table_children_acc, results_table_children_trust, results_table_children_adherence, results_table_children_unsure)
results_table_children_with_metrics_unsure <- rbind(results_table_children_acc, results_table_children_acc_unsure, results_table_children_trust, results_table_children_adherence, results_table_children_adherence_unsure, results_table_children_unsure)
print(results_table_children)
print(results_table_children_with_metrics_unsure)

#### Save results to file
save_table(results_table_children_acc, "children_misinfo_interventions_correct")
save_table(results_table_children_acc_unsure, "children_misinfo_interventions_correct_unsure")
save_table(results_table_children_trust, "children_misinfo_interventions_trust")
save_table(results_table_children_adherence, "children_misinfo_interventions_adherence")
save_table(results_table_children_adherence_unsure, "children_misinfo_interventions_adherence_unsure")
save_table(results_table_children_unsure, "children_misinfo_interventions_unsure")

##### combine results saving
save_table(results_table_children, "children_misinfo_interventions")
save_df(children_df, "children_misinfo_interventions")
save_table(results_table_children_with_metrics_unsure, "children_misinfo_interventions_with_metrics_unsure")
save_df(children_df, "children_misinfo_interventions_with_metrics_unsure")


## Result 6.2: Misinfo false questions only
### Data - all interventions, data where Misinfo is false
true_children_df <- data %>% filter(Child == "True", Tag == "breastfeed" | Tag == "cosleep" | Tag == "rsv" | Tag == "formula", Misinfo == "False") # Children related

### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
true_children_model_acc <- glmer(Correct ~ Intervention_ID + (1| Participant_ID), 
                    data = true_children_df, 
                    family = binomial)
true_children_model_acc_unsure <- glmer(Correct+Unsure ~ Intervention_ID + (1| Participant_ID),
                    data = true_children_df, 
                    family = binomial)
true_children_model_trust <- glmer(Trust ~ Intervention_ID + (1| Participant_ID),
                    data = true_children_df, 
                    family = binomial)
true_children_model_adherence <- glmer(Adherence ~ Intervention_ID + (1| Participant_ID),
                    data = true_children_df, 
                    family = binomial)
true_children_model_adherence_unsure <- glmer(Adherence+Unsure ~ Intervention_ID + (1| Participant_ID),
                    data = true_children_df, 
                    family = binomial)
true_children_model_unsure <- glmer(Unsure ~ Intervention_ID + (1| Participant_ID),
                    data = true_children_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(true_children_model_acc)
check_model_convergence(true_children_model_acc_unsure)
check_model_convergence(true_children_model_trust)
check_model_convergence(true_children_model_adherence)
check_model_convergence(true_children_model_adherence_unsure)
check_model_convergence(true_children_model_unsure)

### Get results
summary_true_children_acc <- summary(true_children_model_acc)
summary_true_children_acc_unsure <- summary(true_children_model_acc_unsure)
summary_true_children_trust <- summary(true_children_model_trust)
summary_true_children_adherence <- summary(true_children_model_adherence)
summary_true_children_adherence_unsure <- summary(true_children_model_adherence_unsure)
summary_true_children_unsure <- summary(true_children_model_unsure)

print("True Children Discernment")

results_table_true_children_acc <- get_results_table(summary_true_children_acc, "Correct", "Intervention")
results_table_true_children_acc_unsure <- get_results_table(summary_true_children_acc_unsure, "Correct+Unsure", "Intervention")
results_table_true_children_trust <- get_results_table(summary_true_children_trust, "Trust", "Intervention")
results_table_true_children_adherence <- get_results_table(summary_true_children_adherence, "Adherence", "Intervention")
results_table_true_children_adherence_unsure <- get_results_table(summary_true_children_adherence_unsure, "Adherence+Unsure", "Intervention")
results_table_true_children_unsure <- get_results_table(summary_true_children_unsure, "Unsure", "Intervention")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table_true_children <- rbind(results_table_true_children_acc, results_table_true_children_trust, results_table_true_children_adherence, results_table_true_children_unsure)
results_table_true_children_with_metrics_unsure <- rbind(results_table_true_children_acc, results_table_true_children_acc_unsure, results_table_true_children_trust, results_table_true_children_adherence, results_table_true_children_adherence_unsure, results_table_true_children_unsure)
print(results_table_true_children)
print(results_table_true_children_with_metrics_unsure)

#### Save results to file
save_table(results_table_true_children_acc, "children_no_misinfo_interventions_correct")
save_table(results_table_true_children_acc_unsure, "children_no_misinfo_interventions_correct_unsure")
save_table(results_table_true_children_trust, "children_no_misinfo_interventions_trust")
save_table(results_table_true_children_adherence, "children_no_misinfo_interventions_adherence")
save_table(results_table_true_children_adherence_unsure, "children_no_misinfo_interventions_adherence_unsure")
save_table(results_table_true_children_unsure, "children_no_misinfo_interventions_unsure")

##### combine results saving
save_table(results_table_true_children, "children_no_misinfo_interventions")
save_df(true_children_df, "children_no_misinfo_interventions")
save_table(results_table_true_children_with_metrics_unsure, "children_no_misinfo_interventions_with_metrics_unsure")
save_df(true_children_df, "children_no_misinfo_interventions_with_metrics_unsure")

## Result 7: Effect of intervention on the children demography with strict correctness

### Result 7.1: Misinfo true questions only
### Data - all interventions, data where Misinfo is true
strict_children_df <- data %>% filter(Child == "True", Unsure == 0, Tag == "breastfeed" | Tag == "cosleep" | Tag == "rsv" | Tag == "formula") # Children related

### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
strict_children_model_acc <- glmer(Correct ~ Intervention_ID + (1| Participant_ID), 
                    data = strict_children_df, 
                    family = binomial)
strict_children_model_trust <- glmer(Trust ~ Intervention_ID + (1| Participant_ID),
                    data = strict_children_df, 
                    family = binomial)  
strict_children_model_adherence <- glmer(Adherence ~ Intervention_ID + (1| Participant_ID),
                    data = strict_children_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(strict_children_model_acc)
check_model_convergence(strict_children_model_trust)
check_model_convergence(strict_children_model_adherence)

### Get results
summary_strict_children_acc <- summary(strict_children_model_acc)
summary_strict_children_trust <- summary(strict_children_model_trust)
summary_strict_children_adherence <- summary(strict_children_model_adherence)

print("Strict Children Discernment")

results_table_strict_children_acc <- get_results_table(summary_strict_children_acc, "Correct", "Intervention")
results_table_strict_children_trust <- get_results_table(summary_strict_children_trust, "Trust", "Intervention")
results_table_strict_children_adherence <- get_results_table(summary_strict_children_adherence, "Adherence", "Intervention")

print(results_table_strict_children_acc)
print(results_table_strict_children_trust)
print(results_table_strict_children_adherence)

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table_strict_children <- rbind(results_table_strict_children_acc, results_table_strict_children_trust, results_table_strict_children_adherence)
print(results_table_strict_children)
results_table_strict_children_with_metrics_unsure <- rbind(results_table_strict_children_acc, results_table_strict_children_trust, results_table_strict_children_adherence)
print(results_table_strict_children_with_metrics_unsure)

#### Save results to file
save_table(results_table_strict_children_acc, "strict_children_misinfo_interventions_correct")
save_table(results_table_strict_children_trust, "strict_children_misinfo_interventions_trust")
save_table(results_table_strict_children_adherence, "strict_children_misinfo_interventions_adherence")

##### combine results saving
save_table(results_table_strict_children, "strict_children_misinfo_interventions")
save_df(strict_children_df, "strict_children_misinfo_interventions")

## Result 7.2: Misinfo false questions only
### Data - all interventions, data where Misinfo is false
strict_true_children_df <- data %>% filter(Child == "True", Unsure == 0, Tag == "breastfeed" | Tag == "cosleep" | Tag == "rsv" | Tag == "formula", Misinfo == "False") # Children related

### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
strict_true_children_model_acc <- glmer(Correct ~ Intervention_ID + (1| Participant_ID), 
                    data = strict_true_children_df, 
                    family = binomial)
strict_true_children_model_trust <- glmer(Trust ~ Intervention_ID + (1| Participant_ID),
                    data = strict_true_children_df, 
                    family = binomial)
strict_true_children_model_adherence <- glmer(Adherence ~ Intervention_ID + (1| Participant_ID),
                    data = strict_true_children_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(strict_true_children_model_acc)
check_model_convergence(strict_true_children_model_trust)
check_model_convergence(strict_true_children_model_adherence)

### Get results
summary_strict_true_children_acc <- summary(strict_true_children_model_acc)
summary_strict_true_children_trust <- summary(strict_true_children_model_trust)
summary_strict_true_children_adherence <- summary(strict_true_children_model_adherence)

print("Strict True Children Discernment")

results_table_strict_true_children_acc <- get_results_table(summary_strict_true_children_acc, "Correct", "Intervention")
results_table_strict_true_children_trust <- get_results_table(summary_strict_true_children_trust, "Trust", "Intervention")
results_table_strict_true_children_adherence <- get_results_table(summary_strict_true_children_adherence, "Adherence", "Intervention")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table_strict_true_children <- rbind(results_table_strict_true_children_acc, results_table_strict_true_children_trust, results_table_strict_true_children_adherence)
print(results_table_strict_true_children)

#### Save results to file
save_table(results_table_strict_true_children_acc, "strict_children_no_misinfo_interventions_correct")
save_table(results_table_strict_true_children_trust, "strict_children_no_misinfo_interventions_trust")
save_table(results_table_strict_true_children_adherence, "strict_children_no_misinfo_interventions_adherence")

##### combine results saving
save_table(results_table_strict_true_children, "strict_children_no_misinfo_interventions")
save_df(strict_true_children_df, "strict_children_no_misinfo_interventions")


## Result 8: Familiarity with misinformation

### Result 8.1: Misinfo true questions only
### Data - all interventions, data where Misinfo is true
familiarity_df <- data %>% filter(Familiarity == "True", Misinfo == "True")

### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
familiarity_model_acc <- glmer(Correct ~ Intervention_ID + (1|Tag) + (1| Participant_ID), 
                    data = familiarity_df, 
                    family = binomial)
familiarity_model_acc_unsure <- glmer(Correct+Unsure ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = familiarity_df, 
                    family = binomial)
familiarity_model_trust <- glmer(Trust ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = familiarity_df, 
                    family = binomial)
familiarity_model_adherence <- glmer(Adherence ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = familiarity_df, 
                    family = binomial)
familiarity_model_adherence_unsure <- glmer(Adherence+Unsure ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = familiarity_df, 
                    family = binomial)
familiarity_model_unsure <- glmer(Unsure ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = familiarity_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(familiarity_model_acc)
check_model_convergence(familiarity_model_acc_unsure)
check_model_convergence(familiarity_model_trust)
check_model_convergence(familiarity_model_adherence)
check_model_convergence(familiarity_model_adherence_unsure)
check_model_convergence(familiarity_model_unsure)

### Get results
summary_familiarity_acc <- summary(familiarity_model_acc)
summary_familiarity_acc_unsure <- summary(familiarity_model_acc_unsure)
summary_familiarity_trust <- summary(familiarity_model_trust)
summary_familiarity_adherence <- summary(familiarity_model_adherence)
summary_familiarity_adherence_unsure <- summary(familiarity_model_adherence_unsure)
summary_familiarity_unsure <- summary(familiarity_model_unsure)

print("Familiarity Discernment")

results_table_familiarity_acc <- get_results_table(summary_familiarity_acc, "Correct", "Intervention")
results_table_familiarity_acc_unsure <- get_results_table(summary_familiarity_acc_unsure, "Correct+Unsure", "Intervention")
results_table_familiarity_trust <- get_results_table(summary_familiarity_trust, "Trust", "Intervention")
results_table_familiarity_adherence <- get_results_table(summary_familiarity_adherence, "Adherence", "Intervention")
results_table_familiarity_adherence_unsure <- get_results_table(summary_familiarity_adherence_unsure, "Adherence+Unsure", "Intervention")
results_table_familiarity_unsure <- get_results_table(summary_familiarity_unsure, "Unsure", "Intervention")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table_familiarity <- rbind(results_table_familiarity_acc, results_table_familiarity_trust, results_table_familiarity_adherence, results_table_familiarity_unsure)
results_table_familiarity_with_metrics_unsure <- rbind(results_table_familiarity_acc, results_table_familiarity_acc_unsure, results_table_familiarity_trust, results_table_familiarity_adherence, results_table_familiarity_adherence_unsure, results_table_familiarity_unsure)
print(results_table_familiarity)
print(results_table_familiarity_with_metrics_unsure)

#### Save results to file
save_table(results_table_familiarity_acc, "familiarity_misinfo_interventions_correct")
save_table(results_table_familiarity_acc_unsure, "familiarity_misinfo_interventions_correct_unsure")
save_table(results_table_familiarity_trust, "familiarity_misinfo_interventions_trust")
save_table(results_table_familiarity_adherence, "familiarity_misinfo_interventions_adherence")
save_table(results_table_familiarity_adherence_unsure, "familiarity_misinfo_interventions_adherence_unsure")
save_table(results_table_familiarity_unsure, "familiarity_misinfo_interventions_unsure")

##### combine results saving
save_table(results_table_familiarity, "familiarity_misinfo_interventions")
save_df(familiarity_df, "familiarity_misinfo_interventions")
save_table(results_table_familiarity_with_metrics_unsure, "familiarity_misinfo_interventions_with_metrics_unsure")
save_df(familiarity_df, "familiarity_misinfo_interventions_with_metrics_unsure")

### Result 8.2: Misinfo false questions only
### Data - all interventions, data where Misinfo is false

true_familiarity_df <- data %>% filter(Familiarity == "True", Misinfo == "False") 

### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
true_familiarity_model_acc <- glmer(Correct ~ Intervention_ID + (1|Tag) + (1| Participant_ID), 
                    data = true_familiarity_df, 
                    family = binomial)
true_familiarity_model_acc_unsure <- glmer(Correct+Unsure ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = true_familiarity_df, 
                    family = binomial)
true_familiarity_model_trust <- glmer(Trust ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = true_familiarity_df, 
                    family = binomial)
true_familiarity_model_adherence <- glmer(Adherence ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = true_familiarity_df, 
                    family = binomial)
true_familiarity_model_adherence_unsure <- glmer(Adherence+Unsure ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = true_familiarity_df, 
                    family = binomial)
true_familiarity_model_unsure <- glmer(Unsure ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = true_familiarity_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(true_familiarity_model_acc)
check_model_convergence(true_familiarity_model_acc_unsure)
check_model_convergence(true_familiarity_model_trust)
check_model_convergence(true_familiarity_model_adherence)
check_model_convergence(true_familiarity_model_adherence_unsure)
check_model_convergence(true_familiarity_model_unsure)

### Get results
summary_true_familiarity_acc <- summary(true_familiarity_model_acc)
summary_true_familiarity_acc_unsure <- summary(true_familiarity_model_acc_unsure)
summary_true_familiarity_trust <- summary(true_familiarity_model_trust)
summary_true_familiarity_adherence <- summary(true_familiarity_model_adherence)
summary_true_familiarity_adherence_unsure <- summary(true_familiarity_model_adherence_unsure)
summary_true_familiarity_unsure <- summary(true_familiarity_model_unsure)

print("True Familiarity Discernment")

results_table_true_familiarity_acc <- get_results_table(summary_true_familiarity_acc, "Correct", "Intervention")
results_table_true_familiarity_acc_unsure <- get_results_table(summary_true_familiarity_acc_unsure, "Correct+Unsure", "Intervention")
results_table_true_familiarity_trust <- get_results_table(summary_true_familiarity_trust, "Trust", "Intervention")
results_table_true_familiarity_adherence <- get_results_table(summary_true_familiarity_adherence, "Adherence", "Intervention")
results_table_true_familiarity_adherence_unsure <- get_results_table(summary_true_familiarity_adherence_unsure, "Adherence+Unsure", "Intervention")
results_table_true_familiarity_unsure <- get_results_table(summary_true_familiarity_unsure, "Unsure", "Intervention")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table_true_familiarity <- rbind(results_table_true_familiarity_acc, results_table_true_familiarity_trust, results_table_true_familiarity_adherence, results_table_true_familiarity_unsure)
results_table_true_familiarity_with_metrics_unsure <- rbind(results_table_true_familiarity_acc, results_table_true_familiarity_acc_unsure, results_table_true_familiarity_trust, results_table_true_familiarity_adherence, results_table_true_familiarity_adherence_unsure, results_table_true_familiarity_unsure)
print(results_table_true_familiarity)
print(results_table_true_familiarity_with_metrics_unsure)

#### Save results to file
save_table(results_table_true_familiarity_acc, "familiarity_no_misinfo_interventions_correct")
save_table(results_table_true_familiarity_acc_unsure, "familiarity_no_misinfo_interventions_correct_unsure")
save_table(results_table_true_familiarity_trust, "familiarity_no_misinfo_interventions_trust")
save_table(results_table_true_familiarity_adherence, "familiarity_no_misinfo_interventions_adherence")
save_table(results_table_true_familiarity_adherence_unsure, "familiarity_no_misinfo_interventions_adherence_unsure")
save_table(results_table_true_familiarity_unsure, "familiarity_no_misinfo_interventions_unsure")

##### combine results saving
save_table(results_table_true_familiarity, "familiarity_no_misinfo_interventions")
save_df(true_familiarity_df, "familiarity_no_misinfo_interventions")
save_table(results_table_true_familiarity_with_metrics_unsure, "familiarity_no_misinfo_interventions_with_metrics_unsure")
save_df(true_familiarity_df, "familiarity_no_misinfo_interventions_with_metrics_unsure")

## Result 9: Effect of intervention on the familiarity demography with strict correctness

### Result 9.1: Misinfo true questions only
### Data - all interventions, data where Misinfo is true
strict_familiarity_df <- data %>% filter(Familiarity == "True", Unsure == 0) # Familiarity related
### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
strict_familiarity_model_acc <- glmer(Correct ~ Intervention_ID + (1|Tag) + (1| Participant_ID), 
                    data = strict_familiarity_df, 
                    family = binomial)
strict_familiarity_model_trust <- glmer(Trust ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = strict_familiarity_df, 
                    family = binomial)
strict_familiarity_model_adherence <- glmer(Adherence ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = strict_familiarity_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(strict_familiarity_model_acc)
check_model_convergence(strict_familiarity_model_trust)
check_model_convergence(strict_familiarity_model_adherence)

### Get results
summary_strict_familiarity_acc <- summary(strict_familiarity_model_acc)
summary_strict_familiarity_trust <- summary(strict_familiarity_model_trust)
summary_strict_familiarity_adherence <- summary(strict_familiarity_model_adherence)

print("Strict Familiarity Discernment")

results_table_strict_familiarity_acc <- get_results_table(summary_strict_familiarity_acc, "Correct", "Intervention")
results_table_strict_familiarity_trust <- get_results_table(summary_strict_familiarity_trust, "Trust", "Intervention")
results_table_strict_familiarity_adherence <- get_results_table(summary_strict_familiarity_adherence, "Adherence", "Intervention")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table_strict_familiarity <- rbind(results_table_strict_familiarity_acc, results_table_strict_familiarity_trust, results_table_strict_familiarity_adherence)
print(results_table_strict_familiarity)

#### Save results to file
save_table(results_table_strict_familiarity_acc, "strict_familiarity_misinfo_interventions_correct")
save_table(results_table_strict_familiarity_trust, "strict_familiarity_misinfo_interventions_trust")
save_table(results_table_strict_familiarity_adherence, "strict_familiarity_misinfo_interventions_adherence")

##### combine results saving
save_table(results_table_strict_familiarity, "strict_familiarity_misinfo_interventions")
save_df(strict_familiarity_df, "strict_familiarity_misinfo_interventions")

### Result 9.2: Misinfo false questions only
### Data - all interventions, data where Misinfo is false
strict_true_familiarity_df <- data %>% filter(Familiarity == "True", Unsure == 0, Misinfo == "False")

### Modeling variables - Participant_ID, Tag: random, Intervention_ID: fixed, Correct/Trust/Adherence/Unsure: observed
strict_true_familiarity_model_acc <- glmer(Correct ~ Intervention_ID + (1|Tag) + (1| Participant_ID), 
                    data = strict_true_familiarity_df, 
                    family = binomial)
strict_true_familiarity_model_trust <- glmer(Trust ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = strict_true_familiarity_df, 
                    family = binomial)
strict_true_familiarity_model_adherence <- glmer(Adherence ~ Intervention_ID + (1|Tag) + (1| Participant_ID),
                    data = strict_true_familiarity_df, 
                    family = binomial)

### Check model convergence
check_model_convergence(strict_true_familiarity_model_acc)
check_model_convergence(strict_true_familiarity_model_trust)
check_model_convergence(strict_true_familiarity_model_adherence)

### Get results
summary_strict_true_familiarity_acc <- summary(strict_true_familiarity_model_acc)
summary_strict_true_familiarity_trust <- summary(strict_true_familiarity_model_trust)
summary_strict_true_familiarity_adherence <- summary(strict_true_familiarity_model_adherence)

print("Strict True Familiarity Discernment")

results_table_strict_true_familiarity_acc <- get_results_table(summary_strict_true_familiarity_acc, "Correct", "Intervention")
results_table_strict_true_familiarity_trust <- get_results_table(summary_strict_true_familiarity_trust, "Trust", "Intervention")
results_table_strict_true_familiarity_adherence <- get_results_table(summary_strict_true_familiarity_adherence, "Adherence", "Intervention")

# Combine results - row format: Fixed Effect, Observed variable, OR, LB, UB
results_table_strict_true_familiarity <- rbind(results_table_strict_true_familiarity_acc, results_table_strict_true_familiarity_trust, results_table_strict_true_familiarity_adherence)
print(results_table_strict_true_familiarity)

#### Save results to file
save_table(results_table_strict_true_familiarity_acc, "strict_familiarity_no_misinfo_interventions_correct")
save_table(results_table_strict_true_familiarity_trust, "strict_familiarity_no_misinfo_interventions_trust")
save_table(results_table_strict_true_familiarity_adherence, "strict_familiarity_no_misinfo_interventions_adherence")

##### combine results saving
save_table(results_table_strict_true_familiarity, "strict_familiarity_no_misinfo_interventions")
save_df(strict_true_familiarity_df, "strict_familiarity_no_misinfo_interventions")

## Result 10: Per-topic Intervention Impact
options(max.print = 100000)
tags <- unique(as.character(data$Tag))
all_tag_results <- list()

get_results_table_glm <- function(model_summary, observed, fixed) {
  coefs <- model_summary$coefficients
  if (is.null(coefs) || nrow(coefs) < 2) return(NULL)
  
  results <- exp(cbind(
    OR = coefs[, 1],
    LB = coefs[, 1] - coefs[, 2] * qnorm(0.975),
    UB = coefs[, 1] + coefs[, 2] * qnorm(0.975)
  ))
  
  res_df <- data.frame(
    Term = rownames(results),
    Metric = observed,
    OR = as.numeric(results[, 1]),
    LB = as.numeric(results[, 2]),
    UB = as.numeric(results[, 3]),
    stringsAsFactors = FALSE
  )
  res_df[-1, ]
}

for (tag in tags) {
  tag_df <- data %>% filter(as.character(Tag) == tag)
  if (length(unique(tag_df$Intervention_ID)) < 2) next
  
  tryCatch({
    tag_model_acc        <- glm(Correct       ~ Intervention_ID, data = tag_df, family = binomial)
    tag_model_acc_unsure <- glm(CorrectUnsure ~ Intervention_ID, data = tag_df, family = binomial)
    tag_model_trust      <- glm(Trust         ~ Intervention_ID, data = tag_df, family = binomial)
    tag_model_adherence  <- glm(Adherence     ~ Intervention_ID, data = tag_df, family = binomial)
    tag_model_unsure     <- glm(Unsure        ~ Intervention_ID, data = tag_df, family = binomial)
    
    results_tag_acc        <- get_results_table_glm(summary(tag_model_acc),        "Correct",        "Intervention")
    results_tag_acc_unsure <- get_results_table_glm(summary(tag_model_acc_unsure), "Correct+Unsure", "Intervention")
    results_tag_trust      <- get_results_table_glm(summary(tag_model_trust),      "Trust",          "Intervention")
    results_tag_adherence  <- get_results_table_glm(summary(tag_model_adherence),  "Adherence",      "Intervention")
    results_tag_unsure     <- get_results_table_glm(summary(tag_model_unsure),     "Unsure",         "Intervention")
    
    results_combined <- rbind(results_tag_acc, results_tag_acc_unsure, results_tag_trust, results_tag_adherence, results_tag_unsure)
    
    if (!is.null(results_combined)) {
      results_combined$Tag <- tag
      all_tag_results[[tag]] <- results_combined
    }
  }, error = function(e) {
    message(paste("Skipping tag:", tag, "- Error:", e$message))
  })
}

if (length(all_tag_results) > 0) {
  final_tag_results <- as.data.frame(do.call(rbind, all_tag_results))
  
  significant_results <- final_tag_results %>%
    filter((LB > 1 | UB < 1)) %>%
    filter(OR < 1e6 & OR > 1e-6) %>%
    select(Tag, Metric, Term, OR, LB, UB) %>%
    arrange(Tag, Metric)

  print(as.data.frame(significant_results), row.names = FALSE)
}

## Result 11: Per-post Priming Group Effect
correct_answers <- c(
  "There is an RSV vaccine recommended for all babies whose mothers did not receive an RSV vaccine during pregnancy.",
  "No, it is not safe",
  "Increase the risk of developing breast cancer",
  "6 Weeks", 
  "Toddler formulas are not needed to meet nutritional needs of young children",
  "Every five years",
  "40",
  "IUDs do not cause long term infertility",
  "Only breastmilk",
  "Anyone who's smoked at least 20 pack years of cigarettes and are between 50 - 80 years of age."
)

# Question mapping
initial_fields <- c(
  "rsv_inital", "cosleep_inital", "hrt_inital", "texas_inital", "formula_inital",
  "papsmear_inital", "mammogram_inital", "IUD_inital", "breastfeed_inital", "lung_inital"
)

post_tags <- c(
  "rsv", "cosleep", "hrt", "texas", "formula",
  "papsmear", "mammogram", "IUD", "breastfeed", "lung"
)
correct_answers_named <- setNames(correct_answers, post_tags)

initial_long_list <- list()

for (i in 1:length(initial_fields)) {
  init_field <- initial_fields[i]
  post_tag <- post_tags[i]
  correct_answer <- correct_answers_named[[post_tag]]
  
  if (init_field %in% names(initial)) {
    
    # Extract the initial response column
    temp <- data.frame(
      user_id = initial$user_id,
      initial_response = initial[[init_field]],
      Tag = post_tag,
      stringsAsFactors = FALSE
    )
    
    # Score the response
    temp <- temp %>%
  mutate(
    Unsure_Initial = as.integer(is.na(initial_response) | 
                                initial_response == "" | 
                                initial_response == "Unsure"),
    
    Correct_Initial = as.integer(
      !Unsure_Initial & 
      trimws(initial_response) == trimws(correct_answer)
    )
  ) %>%
  select(user_id, Tag, Correct_Initial, Unsure_Initial)
    
    initial_long_list[[i]] <- temp
  }
}

initial_long <- bind_rows(initial_long_list)

post_data <- data %>%
  filter(Intervention_ID == "Priming") %>%
  filter(Tag %in% post_tags) %>%
  dplyr::rename(user_id = Participant_ID) %>%
  dplyr::rename(Correct_Post = Correct, Unsure_Post = Unsure)

change_data <- initial_long %>%
  inner_join(post_data, by = c("user_id", "Tag")) %>%
  mutate(
    # Calculate changes
    Accuracy_Change = Correct_Post - Correct_Initial,
    Unsure_Change = Unsure_Post - Unsure_Initial,
    
    # Binary outcomes
    Improved = as.integer(Accuracy_Change > 0),
    Worsened = as.integer(Accuracy_Change < 0),
    Became_More_Certain = as.integer(Unsure_Change < 0),
    
    # Convert to factors for modeling
    user_id = as.factor(user_id),
    Tag = as.factor(Tag),
    Familiarity = as.factor(Familiarity),
    Child = as.factor(Child),
    BreastCancer = as.factor(BreastCancer),
    Misinfo = as.factor(Misinfo)
  )

  get_or_table <- function(model) {
  coef_summary <- summary(model)$coefficients
  or_table <- data.frame(
    Variable = rownames(coef_summary),
    Estimate = coef_summary[, "Estimate"],
    SE = coef_summary[, "Std. Error"],
    OR = exp(coef_summary[, "Estimate"]),
    LB_95 = exp(coef_summary[, "Estimate"] - 1.96 * coef_summary[, "Std. Error"]),
    UB_95 = exp(coef_summary[, "Estimate"] + 1.96 * coef_summary[, "Std. Error"]),
    Z = coef_summary[, "z value"],
    P = coef_summary[, "Pr(>|z|)"],
    Sig = ifelse(coef_summary[, "Pr(>|z|)"] < 0.001, "***",
           ifelse(coef_summary[, "Pr(>|z|)"] < 0.01, "**",
           ifelse(coef_summary[, "Pr(>|z|)"] < 0.05, "*",
           ifelse(coef_summary[, "Pr(>|z|)"] < 0.1, ".", ""))))
  )
  rownames(or_table) <- NULL
  return(or_table)
}

cat("Total observations:", nrow(change_data), "\n")
cat("Number of participants:", length(unique(change_data$user_id)), "\n")
cat("Number of questions:", length(unique(change_data$Tag)), "\n\n")

cat("Unsure Rates:\n")
cat("  Initial unsure:", sprintf("%.1f%%", mean(change_data$Unsure_Initial)*100), "\n")
cat("  Post unsure:", sprintf("%.1f%%", mean(change_data$Unsure_Post)*100), "\n")


question_summary <- change_data %>%
  dplyr::group_by(Tag) %>%
  dplyr::summarise(
    
    Initial_Acc = sprintf("%.1f%%", mean(Correct_Initial, na.rm = TRUE) * 100),
    Post_Acc = sprintf("%.1f%%", mean(Correct_Post, na.rm = TRUE) * 100),
    
    .groups = "drop"
  )
print(question_summary)