import getDemoAnswersMapped from './demographic_answers_mapped.js';

const getDemographicQuestions = () => {
    const demo_answers_mapped = getDemoAnswersMapped();

    // console.log(answers_mapped)

    const demo_questions = [
        {
            code: "gender",
            text: "What is your gender?",
            options: demo_answers_mapped["gender"]
        },
        {
            code: "ethnicity",
            text: "What do you consider to be your race?",
            options: demo_answers_mapped["ethnicity"]
        },
        {
            code: "hispanic",
            text: "Are you of Hispanic, Latino, or Spanish origin?",
            options: demo_answers_mapped["hispanic"]

        },
        {
            code: "kids",
            text: "Do you have any children? If so, how many?",
            options: demo_answers_mapped["kids"]

        },
        {
            code: "tech",
            text: "How much do you agree with the following statement: Artificial intelligence can have a beneficial impact in medicine",
            options: demo_answers_mapped["tech"]

        },
        {
            code: "ai_familiarity",
            text: "How familiar are you with artificial intelligence?",
            options: demo_answers_mapped["ai_familiarity"]

        },
        {
            code: "employment",
            text: "Which of the following best describes your current employment status?",
            options: demo_answers_mapped["employment"]

        },
        {
            code: "media",
            text: "What medium do you most normally use for receiving health information? ",
            options: demo_answers_mapped["media"]

        },
        {
            code: "chatgpt",
            text: "Have you ever asked ChatGPT or another Large Language Model about medical questions?",
            options: demo_answers_mapped["chatgpt"]

        },
        {
            code: "chatgpt2",
            text: "If you said yes to the previous question, how often do you ask ChatGPT or another Large Language Model about medical questions?",
            options: demo_answers_mapped["chatgpt2"]

        },
        {
            code: "lung_demo",
            text: "Do you know any family or friends that have or have had lung cancer?",
            options: demo_answers_mapped["lung_demo"]

        },
        {
            code: "breast",
            text: "Do you know any family or friends that have or have had breast cancer?",
            options: demo_answers_mapped["breast"]

        },
    ] 
    return demo_questions;
}

export default getDemographicQuestions;