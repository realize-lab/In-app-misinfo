const getTrialAnswersMapped = () => {

    const trial_answers_mapped = 
            {
            "trial1": [
              { id: 0, text: "The flu",},
              { id: 1, text: "RSV",},
              { id: 2, text: "Proper Carseat Usage"},
              { id: 3, text: "Postpartum care",},
              { id: 4, text: "Infertility", },
                ],
            "trial2": [
                { id: 0, text: "The vaccine was pulled due to peer-reviewed studies that proved that it caused arthritis in humans",},
                { id: 1, text: "The vaccine was pulled due to the possibility of causing arthritis in hamsters",},
                { id: 3, text: "The vaccine was pulled due to negative backlash since widespread side effects were shown in humans."},
                { id: 2, text: "The vaccine was pulled due to negative backlash from concerns related to possibly causing arthritis in humans"},
                { id: 4, text: "Unsure"},

                    ],
            "trial3": [
                { id: 2, text: "Yes",},
                { id: 1, text: "No",},
                { id: 0, text: "Unsure"},
                    ],
            }
    return trial_answers_mapped
}

export default getTrialAnswersMapped;