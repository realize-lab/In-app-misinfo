const getDemoAnswersMapped = () => {

    const demo_answers_mapped = 
            {
            "ethnicity": [
              { id: 0, text: "White",},
              { id: 1, text: "Black",},
              { id: 3, text: "Asian",},
              { id: 4, text: "2 or more races", },
              { id: 2, text: "Other", }
                ],
            "gender": [
                { id: 0, text: "Male",},
                { id: 1, text: "Female",},
                { id: 2, text: "Other"},
                    ],
            "kids": [
                { id: 0, text: "0 children",},
                { id: 1, text: "1 child",},
                { id: 2, text: "2 children"},
                { id: 3, text: "3+ children",},
                    ],
            "tech": [
                { id: 0, text: "Strongly Agree",},
                { id: 1, text: "Agree",},
                { id: 2, text: "Neither Agree Nor Disagree"},
                { id: 3, text: "Strongly Disagree",},
                { id: 4, text: "Disagree", }
                    ],
            "ai_familiarity": [
                { id: 0, text: "I am involved with designing AI systems or engage in AI research",},
                { id: 1, text: "I actively use AI systems on a regular basis",},
                { id: 2, text: "I have used ChatGPT a few times or have read about AI a fair amount"},
                { id: 3, text: "I have seen a couple article headlines about AI",},
                { id: 4, text: "I have little to no knowledge about AI", }
                    ],
            "chatgpt": [
                { id: 0, text: "Yes, I have asked ChatGPT or another Large Language Model a medical question in the past",},
                { id: 1, text: "No, I have never asked ChatGPT or another Large Language Model a medical question in the past",},
                    ],

            "chatgpt2": [
                { id: 0, text: "Daily",},
                { id: 1, text: "Weekly",},
                { id: 2, text: "Monthly",},
                { id: 3, text: "Occasionally",},
                { id: 4, text: "I have asked ChatGPT or another Large Language Model a medical question once or twice",},
                { id: 5, text: "I do not ask ChatGPT or another Large Language Model medical questions",},
                    ],
            "hispanic": [
                { id: 0, text: "Yes",},
                { id: 1, text: "No",},
                    ],
            "employment": [
                { id: 0, text:  "Employed for wages"},
                { id: 1, text: "Self-employed"},
                { id: 2, text: "Homemaker"},
                { id: 6, text: "Student", },
                { id: 3, text: "Unemployed",},
                { id: 4, text: "Retired", },
                { id: 5, text: "Unable to work", }
                    ],
            "media": [
                { id: 0, text: "Radio",},
                { id: 1, text: "TV",},
                { id: 2, text: "Internet Search: Personal Blogs"},
                { id: 3, text: "Internet Search: Government Websites"},
                { id: 4, text: "Internet Search: Medical Speciality Group Websites"},
                { id: 5, text: "Social Media",},
                { id: 6, text: "TikTok", },
                { id: 7, text: "Other", }
                    ],
            "lung_demo": [
                { id: 0, text: "Yes, self",},
                { id: 1, text: "Yes, family",},
                { id: 2, text: "Yes, friends"},
                { id: 3, text: "No",},
                { id: 4, text: "Other", }
                            ],
            "breast": [
                { id: 0, text: "Yes, self",},
                { id: 1, text: "Yes, family",},
                { id: 2, text: "Yes, friends"},
                { id: 3, text: "No",},
                { id: 4, text: "Other", }
                            ],
            "diet": [
                { id: 0, text: "Strongly Agree",},
                { id: 1, text: "Agree",},
                { id: 2, text: "Neither Agree Nor Disagree"},
                { id: 3, text: "Strongly Disagree",},
                { id: 4, text: "Disagree", }
                    ],  
            }
    return demo_answers_mapped
}

export default getDemoAnswersMapped;
