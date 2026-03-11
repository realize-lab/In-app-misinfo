import json
from collections import Counter 
import numpy
import statsmodels.api as sm
import pandas as pd
import numpy as np
import random
import time
import math
import scipy

# Intervention ID
# False or true
# time_step 
# score
# Participant Id

q = {"order":"3,5,4,1,6,7,0,2","rsv":"There is an emergency use RSV immunoprophylaxis for babies deemed at risk","asthma_outside":"Relied mostly on the AI recommendation to answer","sweets_outside":"Relied mostly on outside knowledge to answer","rsv_outside":"Relied on outside knowledge and the AI recommendation equally to answer","id":"345523","hrt_outside":"Relied mostly on outside knowledge to answer","question times":"1705271426068,1705271431017,1705271441402,1705271445350,1705271450301,1705271456334,1705271460917,1705271468817,1705271473651","sweets_likert":"Unsure","ethnicity":"Latino","cosleep":"Yes, but only with close adult supervision","employment":"Unable to work","rsv_likert":"Somewhat Trustworthy","asthma":"Yes, but only in a supervised setting like hospital or clinic.","cosleep_outside":"Relied entirely on the AI recommendation to answer","breast":"Yes, friends","media":"Internet Search: Personal Blogs","asthma_likert":"Somewhat Trustworthy","lung":"People aged 50-80, who've smoked 20 pack-years, and are either currently smoking or have quit in the past 15 years","texas_2":"Abortion","hrt":"Slightly decrease risk of cancer,","sweets":"Non-sugar sweetners help achieve weight-loss since they reduce cravings for sugar.","start_time":"1705271383880","rsv_2":"31 weeks","asthma_2":"A mixture of active and inactive flu viruses","gender":"Other","kids":"2 kids","formula_likert":"Somewhat Trustworthy","texas":"6 Weeks","hispanic":"No","formula_2":"Has no nutritional benefit","leftTab":"false","accepted":"true","cosleep_2":"Ideally 18 months","cosleep_likert":"Unsure","texas_likert":"Somewhat Trustworthy","texas_outside":"Relied entirely on the AI recommendation to answer","formula":"After one year of age","group":"1","sweets_2":"Decreases the risk","lung_outside":"Relied entirely on the AI recommendation to answer","formula_outside":"Relied on outside knowledge and the AI recommendation equally to answer","hrt_likert":"Unsure","tech":"Neither Agree Nor Disagree","hrt_2":"Greatly increase risk of cancer","lung_likert":"Somewhat Trustworthy","lung_2":"No, they should wait until 25 years since quitting to begin "}

def parse_ndjson(data):
    return [json.loads(l) for l in data.splitlines()]


with open('real_study_final.json', 'r', encoding="utf8") as handle:
    data = handle.read()
    dicts = parse_ndjson(data)

 
questions_tags_ns = ['rsv', 'cosleep', 
 #'texas_2',
  'hrt', 'rsv_2', 
 'texas', 'formula_2', 'cosleep_2',
 'formula',  'hrt_2',
 "papsmear", "papsmear_2",
 "mammogram", 
#  "mammogram_2",
 "IUD", "IUD_2", "breastfeed", "breastfeed_2",
 "lung","lung_2", 
#  "asthma", "asthma_2"
 ]

questions_tags_inital = ['rsv', 'cosleep', 'hrt', 'texas', 'formula',  "papsmear", "mammogram", "IUD", "breastfeed",
"lung",
]

answers_about_misinfo = ['texas','hrt_2','formula','formula_2','lung','lung_2','rsv',]

legal = set(["texas"])

cause_and_effect = set(["hrt","hrt_2","rsv_2", "papsmear", "papsmear_2", "IUD"])

timing = set(["formula", "lung", "cosleep_2", "papsmear" , "papsmear_2", "mammogram", "breastfeed_2" ])

medical_practices = set(["formula","formula_2","lung","lung_2","rsv","cosleep","breastfeed","IUD_2"])

reading_comp = set(["hrt","rsv_2","cosleep","cosleep_2", "mammogram", "IUD","IUD_2","breastfeed","breastfeed_2" ])


outside = set(["texas","formula","formula_2","lung","lung_2","rsv","hrt_2"])



misinfo_text_codes = set(
    [
        "texas", 
        "texas_2", 
        "rsv",
        "rsv_2",
        "formula",
        "formula_2",
        "lung",
        "lung_2",
        "hrt",
        "hrt_2"
    ]
)

demographic_tags =  ['ethnicity', 'employment', 'media', 'gender', 'kids',
 'hispanic', 'leftTab', 'accepted', 'group', 'tech']

correct_answers = set([
    "6 Weeks", 
    "Increase the risk of developing breast cancer",
     "Toddler formulas are not needed to meet nutritional needs of young children",
     "People aged 50-80, who have smoked 20 pack years, regardless of when they quit. ",
     "Anyone who's smoked at least 20 pack years of cigarettes and are between 50 - 80 years of age.",
     "There is an RSV vaccine recommended for all babies whose mothers did not receive an RSV vaccine during pregnancy.",
     "Abortion",
     "Increase the risk of developing cancer",
     #Add check to count for breast cancer risk
     "Toddler formula does not provide a nutritional advantage over a well-balanced diet that includes cows milk",
     "Yes, every year.",
     "Pneumonia",
    #  "Non sugar sweetners for weight loss",
     "No, it is not safe",
     "Ideally 12 months",
     "Every five years",
     "21",
     "Progestin",
     "Only breastmilk",
    #  "Breastmilk and milk alternatives",
     "6 months",
     "IUDs do not cause long term infertility",
     "40",
])

adherance_answers = set([
    "20 Weeks",
    "Increase the risk of developing breast cancer",
    "After one year of age",
    "After one year of age if you have stopped breastfeeding your child",
    "Anyone aged 50 - 80, who has smoked 20 pack years, and is currently still smoking or has quit in the past 15 years.",
    "There is no routinely recommended RSV vaccine for infants",
    "There is an emergency use RSV immunoprophylaxis for babies deemed at risk during RSV season",
    "Decrease the risk of developing cancer",
    "Toddler formula includes all nutrients needed under FDA regulation for toddler development past the age of 12 months.",
    "No, they should not get screened because they do not meet all the criteria",
    "Pneumonia",
    "No, it is not safe",
    "Ideally 12 months",
    "Every five years",
    "21",
    "40",
    # "45",
    "IUDs do not cause long term infertility",
    "Progestin",
    "Only breastmilk",
    "6 months",
])


default_order = [
    "texas",
    "hrt",
    "formula",
    "lung",
    "rsv",
    "cosleep",
    "papsmear",
    "mammogram",
    "IUD",
    "breastfeed",
] 

def check_misinfo(tag):
    if tag in misinfo_text_codes:
        return True
    else:
        return False

def check_misinfo_answer(tag):
    if tag in answers_about_misinfo:
        return True
    else:
        return False

def attention_check(d):
    if d["answers"]["texas_2"] == "Abortion" and d['answers']['mammogram_2'] == "Mammograms":
        return True
    else:
        return False


def filter_times(d):
    times = d["answers"]["question_times"]
    t = times.split(",")
    start_time = int(t[0])
    last_time = int(t[-1])
    duration = last_time - start_time
    cutoff = [10.728,11.976,9.392,11.094,10.697]
    return duration * (0.0166667 / 1000) >= cutoff[int(d["answers"]["group"])-1]

def create_csv_table_all():


    master_list = []
    for d in dicts:
        if attention_check(d) and filter_times(d):

            order = d["answers"]["order"]
            order = order.split(",")

        
            csv_file = 'revised_table.csv'

            ##Changehere
            for tag in questions_tags_ns:
                if check_misinfo(tag):
                        misinfo = True
                else:
                        misinfo = False 

                pos_tag = tag.split("_")[0]

                question_pos = default_order.index(pos_tag)

                displayed_order = int(order[question_pos])

                if displayed_order < 3:
                        pos = 0
                elif displayed_order < 7:
                        pos = 1
                else:
                        pos = 2

                # pos = 0
                if True:


                    if (d["answers"][pos_tag + "_likert"] == "Very Trustworthy" or d["answers"][pos_tag + "_likert"] == "Somewhat Trustworthy"):
                        trust = True
                    else:
                        trust = False

                    pos_tag_likert = pos_tag +  "_likert"
                    if d["answers"][pos_tag_likert] == "Very Trustworthy":
                        # score = 4
                        trust_score = 5
                    elif d["answers"][pos_tag_likert] == "Somewhat Trustworthy":
                        trust_score = 4
                    elif d["answers"][pos_tag_likert] == "Neither Trustworthy nor Untrustworthy":
                        trust_score = 3
                    elif d["answers"][pos_tag_likert] == "Somewhat Untrustworthy":
                        trust_score = 2
                    else:
                        trust_score = 1

                    # pos = 0

                    if (d["answers"]["breast"] == "Yes, self" or d["answers"]["breast"] == "Yes, family" or d["answers"]["breast"] == "Yes, friends"):
                        breast = True
                    else:
                        breast = False

                    if (d["answers"]["ai_familiarity"] == "I am involved with designing AI systems or engage in AI research" or 
                    d["answers"]["ai_familiarity"] == "I actively use AI systems on a regular basis"
                    or d["answers"]["ai_familiarity"] == "I have used ChatGPT a few times or have read about AI a fair amount"):
                        familiarity = True
                    else:
                        familiarity = False


                    if d['answers']['kids'] != "0 children":
                        child = True
                    else:
                        child = False


                    ## add question type
                    if check_misinfo_answer(tag) is True:
                        misinfo_answer = True
                    else:
                        misinfo_answer = False


                    correctUnsure = 0
                    adherenceUnsure = 0
                    unsure = 0

                    if d["answers"][tag] in correct_answers:
                        correct = 1
                        correctUnsure = 1
                    else:
                        correct = 0
                    
                    if d["answers"][tag] in adherance_answers:
                        adherence = 1
                        adherenceUnsure = 1
                    else:
                        adherence = 0
                    
                    if d["answers"][tag] == "Unsure":
                        unsure = 1
                        correctUnsure = 1
                        adherenceUnsure = 1



                    intervent = ["Baseline","In-line","Pop-up","RAG",'Priming']
                        
                    new_element = [intervent[int(d["answers"]["group"])-1],tag, misinfo, misinfo_answer , trust , trust_score, child, breast, familiarity,  correct, correctUnsure, adherence, adherenceUnsure, unsure, d["user_id"] ]


                    master_list.append(new_element)



    columns = ['Intervention_ID', "Tag", 'Misinfo', "Misinfo_Answer", 'Trust', "TrustScore", "Child", "BreastCancer" , "Familiarity" ,  'Correct' , "CorrectUnsure" , "Adherence" , "AdherenceUnsure ", "Unsure", 'Participant_ID']



    df = pd.DataFrame(master_list, columns=columns)
    df.to_csv(csv_file, index=False)  


    print(df)
    print(df.head(11))

    return(df)


if __name__ == '__main__':
    create_csv_table_all()
