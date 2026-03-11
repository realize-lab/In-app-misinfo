    const getAnswersMapped = () => {

        const answers_mapped = 
                {
                "texas": [
                { id: 0, text: "6 Weeks",},
                { id: 1, text: "12 Weeks",},
                { id: 2, text: "20 Weeks"},
                { id: 3, text: "24 Weeks",},
                { id: 4, text: "27 Weeks",},
                { id: 5, text: "Unsure", }
                    ],
                "hrt": [
                    { id: 0, text: "Increase the risk of developing breast cancer",},
                    { id: 1, text: "Decrease the risk of developing breast cancer"},
                    { id: 2, text: "No effect"},
                    { id: 3, text: "Unsure", }
                        ],
                "asthma": [
                    { id: 0, text: "Yes, but only if the child with asthma is over the age of 2 years",},
                    { id: 1, text: "Yes, but only if they haven't had a history of wheezing within the last 12 months",},
                    { id: 2, text: "Yes, but only if the child with asthma is at least 5 years old."},
                    { id: 3, text: "No, children with asthma cannot get the nasal vaccine.",},
                    { id: 4, text: "Unsure", }
                        ],
                
                "formula": [
                    { id: 0, text: "After one year of age"},
                    { id: 1, text: "After one year of age if you have stopped breastfeeding your child"},
                    { id: 2, text: "Once your child begins to teethe, between 1-3 years of age"},
                    { id: 3, text: "Toddler formulas are not needed to meet nutritional needs of young children", },
                    { id: 4, text: "Unsure", }
                        ],
                "lung": [
                    { id: 0, text: "Anyone aged 50 - 80."},
                    { id: 1, text: "Anyone who's smoked at least 20 pack years of cigarettes and are between 50 - 80 years of age."},
                    { id: 2, text: "Anyone aged 50 - 80, who has smoked 20 pack years, and is currently still smoking or has quit in the past 15 years."},
                    { id: 3, text: "Anyone aged 50 - 80, who has smoked 20 pack years, regardless of when they quit." },
                    { id: 4, text: "Unsure", }
                        ],
                "likert": [
                    { id: 0, text: "Very Trustworthy"},
                    { id: 1, text: "Somewhat Trustworthy"},
                    { id: 2, text: "Neither Trustworthy nor Untrustworthy"},
                    { id: 3, text: "Somewhat Untrustworthy", },
                    { id: 4, text: "Not Trustworthy", }
                        ],
                "rsv": [
                    { id: 0, text: "There is no routinely recommended RSV vaccine for infants"},
                    { id: 1, text: "There is an emergency use RSV immunoprophylaxis for babies deemed at risk during RSV season"},
                    { id: 2, text: "There is an RSV vaccine recommended for all babies whose mothers did not receive an RSV vaccine during pregnancy."},
                    { id: 3, text: "Unsure"},
                        ],
                // "covid-vax": [
                //     { id: 0, text: "6 Months"},
                //     { id: 1, text: "12 Months"},
                //     { id: 2, text: "5 years"},
                //     { id: 3, text: "10 years", },
                //     { id: 4, text: "12 years", },
                //     { id: 5, text: "Unsure", }
                //         ],
                // "sweets": [
                //     { id: 0, text: "Non-sugar sweetners help achieve weight-loss since it replaces sugar, which has more calories"},
                //     { id: 1, text: "Non-sugar sweetners help achieve weight-loss since they reduce cravings for sugar."},
                //     { id: 2, text: "Unsure"},
                //     { id: 3, text: "Non-sugar sweetners do not confer any long-term benefit in reducing body fat in adults or children"},
                //         ],
                "texas2": [
                    { id: 0, text: "Immigration",},
                    { id: 1, text: "Abortion",},
                    { id: 2, text: "Cancer"},
                    { id: 3, text: "Contraceptives",},
                    { id: 4, text: "Unsure", }
                        ],  
                "hrt2": [
                    { id: 0, text: "No effect",},
                    { id: 1, text: "Increase the risk of developing cancer",},
                    { id: 2, text: "Decrease the risk of developing cancer"},
                    { id: 3, text: "Unsure", }
                        ],
                "asthma2": [
                    { id: 0, text: "The nasal flu vaccine contains activated flu viruses",},
                    { id: 1, text: "The nasal flu vaccine contains inactive flu viruses",},
                    { id: 2, text: "The nasal flu vaccine contains a mixture of active and inactive flu viruses"},
                    { id: 3, text: "Unsure",},
                        ],

                "formula2": [
                    { id: 0, text: "Toddler formula includes all nutrients needed under FDA regulation for toddler development past the age of 12 months."},
                    { id: 1, text: "Toddler formula does not provide a nutritional advantage over a well-balanced diet that includes cows milk"},
                    { id: 2, text: "Toddler formula contains certain nutrients needed for a well-balanced diet for toddlers under 12 months and should be introduced in combination with regular breastfeeding."},
                    { id: 3, text: "Unsure", }
                        ],

                "lung2": [
                    { id: 0, text: "Yes, every year."},
                    { id: 1, text: "Yes, every few years"},
                    { id: 2, text: "No, they should wait until 25 years since quitting to begin "},
                    { id: 3, text: "No, they should not get screened because they do not meet all the criteria"},
                    { id: 4, text: "Unsure", }
                        ],
                "likert2": [
                    { id: 0, text: "Very Trustworthy"},
                    { id: 1, text: "Somewhat Trustworthy"},
                    { id: 2, text: "Unsure"},
                    { id: 3, text: "Somewhat Untrustworthy", },
                    { id: 4, text: "Not Trustworthy", }
                        ],
                "rsv2": [
                    { id: 0, text: "Pneumonia"},
                    { id: 1, text: "Immunodeficiencies"},
                    { id: 2, text: "Decreased bone mass"},
                    { id: 3, text: "Nearsightedness", },
                    { id: 4, text: "Cystic Fibrosis", },
                    { id: 5, text: "Unsure", },
                        ],
                // "covid-vax2": [
                //     { id: 0, text: "6 Months"},
                //     { id: 1, text: "12 Months"},
                //     { id: 2, text: "5 years"},
                //     { id: 3, text: "10 years"},
                //     { id: 4, text: "12 years"},
                //     { id: 5, text: "Unsure"}
                //         ],
                // "sweets2": [
                //     { id: 0, text: "Sugar in food"},
                //     { id: 1, text: "Non sugar sweetners for weight loss"},
                //     { id: 2, text: "Weight loss techniques"},
                //     { id: 3, text: "Healthy diets"},
                //     { id: 4, text: "Cardiovascular diseases", },

                        // ],
                // "outside": [
                //     { id: 0, text: "Relied entirely on outside knowledge to answer"},
                //     { id: 1, text: "Relied mostly on outside knowledge to answer"},
                //     { id: 2, text: "Relied on outside knowledge and the AI conversation equally to answer"},
                //     { id: 3, text: "Relied mostly on the AI conversation to answer"},
                //     { id: 4, text: "Relied entirely on the AI conversation to answer", }
                //         ],
                "cosleep": [
                    { id: 0, text: "Generally yes"},
                    { id: 1, text: "No, it is not safe"},
                    { id: 2, text: "Yes, but only with a firm mattress"},
                    { id: 3, text: "Yes, but only with close adult supervision"},
                    { id: 4, text: "Unsure", }

                        ],
                "cosleep2": [
                    { id: 0, text: "Ideally 3 months"},
                    { id: 1, text: "Ideally 6 months"},
                    { id: 2, text: "Ideally 12 months"},
                    { id: 3, text: "Ideally 18 months"},
                    { id: 4, text: "Unsure", }
                        ],
                "papsmear": [
                    { id: 0, text: "Every year"},
                    { id: 1, text: "Every two years"},
                    { id: 2, text: "Every three years"},
                    { id: 3, text: "Every five years"},
                    { id: 4, text: "Unsure", }
                        ],
                "papsmear2": [
                    { id: 0, text: "21"},
                    { id: 1, text: "30"},
                    { id: 2, text: "40"},
                    { id: 3, text: "65"},
                    { id: 4, text: "Unsure", }
                        ],

                "mammogram": [
                    { id: 0, text: "40"},
                    { id: 1, text: "45"},
                    { id: 2, text: "55"},
                    { id: 3, text: "75"},
                    { id: 4, text: "Unsure", }
                        ],

                "mammogram2": [
                    { id: 0, text: "Cancer"},
                    { id: 1, text: "Mammograms"},
                    { id: 2, text: "Aging"},
                    { id: 3, text: "The USPSTF"},
                    { id: 4, text: "Unsure", }
                        ],

                "IUD": [
                    { id: 0, text: "IUDs have no effect on contraception"},
                    { id: 1, text: "IUDs do not cause long term infertility"},
                    { id: 2, text: "IUDs are associated with long term infertility"},
                    { id: 3, text: "Unsure", }
                        ],
                "IUD2": [
                    { id: 0, text: "Progestin"},
                    { id: 1, text: "Estrogen"},
                    { id: 2, text: "Oxytocin"},
                    { id: 3, text: "Thyroxine"},
                    { id: 4, text: "Unsure", }
                        ],

                "breastfeed": [
                    { id: 0, text: "Only breastmilk"},
                    { id: 1, text: "Water and breastmilk"},
                    { id: 2, text: "Breastmilk and milk alternatives"},
                    { id: 3, text: "Baby formula, water, and breastmilk"},
                    { id: 4, text: "Unsure", }
                        ],

                "breastfeed2": [
                    { id: 0, text: "3 months"},
                    { id: 1, text: "6 months"},
                    { id: 2, text: "12 months"},
                    { id: 3, text: "24 months"},
                    { id: 4, text: "Unsure", }
                        ],   
                }

        return answers_mapped
    }

export default getAnswersMapped;