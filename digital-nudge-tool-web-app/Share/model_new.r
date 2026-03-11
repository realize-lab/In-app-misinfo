library(tidyverse)
library(stringr)
library(broom)
library(lme4)
library(afex)
library(texreg)
library(rstatix)
library(lmerTest)
library(caret)

# library(dypler)


column <- list('Intervention_ID', "Tag", 'Misinfo', "Misinfo_Answer", 'Trust', "Child", "BreastCancer" , 'Correct' , "Correct+Unsure" , "Adherence" , "Adherence+Unsure ", "Unsure", 'Participant_ID')

for (i in 0:0){

    file_name <- "/Users/vicenteriquelme/Code/digital-nudge-tool-clone/share/revised_table.csv"
    print(file_name)
    clean <- read.csv(file_name,stringsAsFactors = FALSE)


    clean$Intervention_ID<-as.factor(clean$Intervention_ID)
    clean$Trust<-as.factor(clean$Trust)
    clean$Misinfo<-as.factor(clean$Misinfo)
    clean$Misinfo_Answer<-as.factor(clean$Misinfo_Answer)
    clean$Tag<-as.factor(clean$Tag)
    clean$Familiarity<-as.factor(clean$Familiarity)
    clean$Participant_ID<-as.factor(clean$Participant_ID)
    clean$BreastCancer<-as.factor(clean$BreastCancer)

    clean$Child<-as.factor(clean$Child)
    # clean$Misinfo_Answer<-as.factor(clean$Misinfo_Answer)

    # view(clean)
    df_filtered <- clean
    # df_filtered <- clean %>% filter(Tag == "rsv" | Tag == "rsv_2" | Tag == "cosleep" | Tag == "cosleep_2" | Tag == "breastfeed" | Tag == "breastfeed_2" | Tag == "formula" | Tag == "formula_2"
    # , Fa == "True")
    # ["rsv","rsv_2","cosleep","cosleep_2","breastfeed","breastfeed_2","formula","formula_2"]
    # df_filtered <- clean %>% filter(BreastCancer == "True", 
    # Tag == "hrt" | Tag == "hrt_2" | Tag == "mammogram",
    # )
    # df_filtered <- clean %>% filter(Misinfo_Answer == "False")
    df_filtered <- clean %>% filter(Intervention_ID == "Baseline")

    # view(df_filtered)

    # view(df_filtered)

    # print("mean")
    # print(mean(df_filtered$Adherence))


    # Correct, CorrectUnsure, Unsure, Adherence,AdherenceUnsure
        glmm_baseline <-  glmer(Correct ~  Misinfo_Answer + Trust + Familiarity +  (1| Participant_ID) + (1| Tag), 
                            data = df_filtered, 
                            family = binomial)

    # glmm_baseline <-  glmer(Correct ~ Misinfo_Answer + (1| Participant_ID) + (1| Tag), 
    #                     data = df_filtered, 
    #                     family = binomial)

    # glmm_baseline <-  lmer(TrustScore ~ Misinfo_Answer  + (1| Participant_ID) + (1| Tag),
    #         data = df_filtered)

    # glmm_baseline <-  lmer(TrustScore ~ + BreastCancer  + (1| Participant_ID) + (1| Tag),
    #         data = df_filtered)

    
    # p0 <- predict(glmm_baseline, df_filtered %>% mutate(Intervention_ID='Baseline'), type="response")
    # p1 <- predict(glmm_baseline, df_filtered %>% mutate(Intervention_ID='In-line'), type="response")

    # p0 <- predict(glmm_baseline, df_filtered %>% mutate(Intervention_ID='Baseline'), type="response")
    # p1 <- predict(glmm_baseline, df_filtered %>% mutate(Intervention_ID='In-line'), type="response")

    # print(mean(p0))
    # print(mean(p1))


# qnorm(0.975)



    # round(100*mean(flights_ex[['late']]), 1)

    summary_lme <- summary(glmm_baseline)

    print(summary_lme)
    print(summary_lme$coefficients)
    # estimates <- summary_lme$coefficients[2:6]
    # std <- summary_lme$coefficients[8:12]
    
    # estimates <- summary_lme$coefficients[2:5]
    # std <- summary_lme$coefficients[7:10]
    # z_values <-  summary_lme$coefficients[12:15]


    estimates <- summary_lme$coefficients[2:5]
    std <- summary_lme$coefficients[6:8]
    # print(estimates)
    print("odds")
    for (e in estimates){
        print("estimate")
        print(e)
        print("odd")
        print(exp(e))
    }
    # for (e in std){
    #     print("estimate")
    #     print(e)
    #     print("odd")
    #     print(exp(e))
    # }
    # # print("std")
    # print(std)


# names <- list("in-line","pop-up","priming","rag","misinfo")

    for (i in 1:length(estimates)) {

    # print(names[i])
    p <- estimates[[i]]
    s <- std[[i]]
    # z <- z_values[i]

    # p <- -1.8571
    # s <- 0.2765

    print("start1")
    print(p)
    print(s)

    high <- p + (1.959964 * s)
    low <- p - (1.959964 * s)

    exp_e <- exp(p)
    print("start2")
    print(exp_e)
    print(exp(low))
    print(exp(high))
    # odds_ratios <- exp(rnorm(1000, mean=p, sd=s))
    # # print(odds_ratios)
    # lb <- quantile(odds_ratios, 0.025)
    # ub <- quantile(odds_ratios, 0.975)
    # print(exp_e)
    # print(lb)
    # print(ub)
    # print(z)
    }
}