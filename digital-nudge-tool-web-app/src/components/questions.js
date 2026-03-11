import getAnswers from './answers.js'
import getAnswersMapped from './answers_mapped.js';

const getQuestions = () => {
    const answers = getAnswers();
    const answers_mapped = getAnswersMapped();

    // console.log(answers_mapped)

    const questions = [
        {
            code: "texas",
            text: "What is the latest week that it is legal to get an abortion in South Carolina?",
            options: answers_mapped["texas"],
            text2: "What is this passage about?",
            options2: answers_mapped["texas2"],
        },
        {
            code: "hrt",
            text: "How does taking combined estrogen and progestin Hormone Replacement Therapy (HRT) affect a healthy individual's chance of developing breast cancer?",
            options: answers_mapped["hrt"],
            text2: "How does estrogen only Hormone Replacement Therapy (HRT) affect a healthy individual's chance of developing any form of cancer?",
            options2: answers_mapped["hrt2"],
        },
        // {
        //     code: "asthma",
        //     text: "Can a child with athsma can get a nasal flu vaccine?",
        //     options: answers_mapped["asthma"],
        //     text2: "The nasal flu vaccine differs from an injectable flu vaccine in what way?",
        //     options2: answers_mapped["asthma2"],

        // },
        {
            code: "formula",
            text: "When should parents begin feeding their child toddler formula to ensure a balanced and complete diet?",
            options: answers_mapped["formula"],
            text2: "What nutritional benefit does toddler formula contain over cows milk?",
            options2: answers_mapped["formula2"],
        },
        {
            code: "lung",
            text: "Which of the following smoker profiles should receive the yearly recommended screening for lung cancer?",
            options: answers_mapped["lung"],
            text2: "Should a 55 year old man, who has smoked 30 pack-years, but who quit 20 years ago, get screened for lung cancer?",
            options2: answers_mapped["lung2"],

        },
        // {
        //     code: "covid-vax",
        //     text: "What is the youngest age a child can get a covid-19 vaccine?",
        //     options: answers_mapped["covid-vax"],
        //     text2: "Placeholder",
        //     options2: answers_mapped["texas2"],

        // },
        // {
        //     code: "sweets",
        //     text: "What role do non-sugar sweetners play in weight-loss?",
        //     options: answers_mapped["sweets"],
        //     text2: "What is this passage about?",
        //     options2: answers_mapped["sweets2"],

        // },
        {
            code: "rsv",
            text: "What is the best way to prevent RSV in infants born during RSV season?",
            options: answers_mapped["rsv"],
            text2: "RSV in infants can be a serious condition. Infants with RSV are more likely to develop which of the following conditions?",
            options2: answers_mapped["rsv2"], 

        },

        {
            code: "cosleep",
            text: "Is it safe to sleep in the same bed as your baby at 6 months?",
            options: answers_mapped["cosleep"],
            text2: "How long should babies sleep in the same room as parents?",
            options2: answers_mapped["cosleep2"], 

        },

        {
            code: "papsmear",
            text: "How often should a 40 year old without any history of abnormal pap-smears get a Pap-smear and HPV co-test done?",
            options: answers_mapped["papsmear"],
            text2: "At what age should an individual begin to get pap-smears?",
            options2: answers_mapped["papsmear2"], 

        },

        {
            code: "mammogram",
            text: "At what age could someone begin mammography?",
            options: answers_mapped["mammogram"],
            text2: "What is this passage most closely about?",
            options2: answers_mapped["mammogram2"], 

        },

        {
            code: "IUD",
            text: "Do IUDs (Intrauterine Devices) cause infertility?",
            options: answers_mapped["IUD"],
            text2: "What hormone, found in hormonal IUDs, limits sperm movement and survival?",
            options2: answers_mapped["IUD2"], 

        },

        {
            code: "breastfeed",
            text: "What kind of liquids can a 3 month old baby safely drink?",
            options: answers_mapped["breastfeed"],
            text2: "At what age is it recommended to introduce complementary foods into an infant's diet?",
            options2: answers_mapped["breastfeed2"], 

        },

        
    ] 
    return questions;
}

export default getQuestions;