import React, { useState, setState, Component, useEffect } from "react";
// import { View, FlatList, StyleSheet, Text, Dimensions } from "react-native";
// import * as React from 'react';  
import styles from "./App.css";
import getQuestions from './components/questions.js'
import getPrompts from "./components/medicalTexts";
// import getPopup from "./components/popup";
// import alertMode from "./components/alertMode";
import PopupAlert from "./components/popupAlert";
import PopupConfirm from "./components/popupConfirm";
import { useParams, useSearchParams } from "react-router-dom"
import getDisclaimer from "./components/disclaimer.js";
import getDemographicQuestions from "./components/demographic_questions";
import getStartPage from "./components/startPage";
import useIsTabVisible from "./components/tracker";
import { isVisible } from "dom-helpers";
import getRAGTexts from "./components/medicalTextsRAG";
import getTrialQuestions from "./components/trialQuestion";
import getAnswersMapped from "./components/answers_mapped";
import { nullLiteral } from "@babel/types";
import GradualTextDisplay from "./components/display_gradually";
import GradualStringList from "./components/gradual_string";
// import { BrowserRouter as Router } from 'react-router-dom';


function App({group, chatgpt, popup, rag, priming}) {


  const [UsePopups, setUsePopups] = useState(popup);
  const [UseChatGPTDisclaimer, setUseChatGPTDisclaimer] = useState(chatgpt);
  const [primed,setPrimed] = useState(priming) 
  const [trial_primed,setTrialPrimed] = useState(false) 


  // const handleChatGPT = () => setUseChatGPTDisclaimer(true)
  // const handlePopup = () => setUsePopups(true);

  // const [useFreq, setFreq] = useState(false);
  // const [freqCounter, setFreqCounter] = useState(0);


  // const handleFreq = () => {
  //   //Query Model to determine whether to show or not
  // }


  // let { id } = useParams()
  const [searchParams] = useSearchParams();
  // console.log(searchParams.get("name"));
  const id = searchParams.get("PROLIFIC_PID")

  // console.log(id, group)

  // console.log(id)


  // Properties
  const [showResults, setShowResults] = useState(false);
  const [showEndScreen, setShowEndScreen] = useState(false);
  const [currentQuestion, setCurrentQuestion] = useState(0);
  const [currentDemoQuestion, setCurrentDemoQuestion] = useState(0);
  const [showPopups, setShowPopup] = useState(true);
  const [showConfirmPopupAccept, setShowConfirmPopupAccept] = useState(false);
  const [showConfirmPopupDecline, setShowConfirmPopupDecline] = useState(false);
  const [showEndTrial,setShowEndTrial] = useState(false)

  // const [key] = useState(0);
  // const responses = [];
  // const [score, setScore] = useState(0);
  const [selectedAnswerIndex, setSelectedAnswerIndex] = useState(null)
  const [selectedAnswerIndex2, setSelectedAnswerIndex2] = useState(null)
  const [selectedAnswerIndexLikert, setSelectedAnswerIndexLikert] = useState(null)
  const [selectedAnswerIndexOutside, setSelectedAnswerIndexOutside] = useState(null)

  const [primed_example, setPrimedExample] = useState(false)


  const [acceptedOrNot, setAcceptedOrNot] = useState(false)

  const questions = getQuestions();
  const trial_questions = getTrialQuestions();
  const answers_mapped = getAnswersMapped();

  // const prompts = getRAGTexts();
  // const prompts = getPrompts();


  // console.log(prompts)

  const post_disclaimer = getDisclaimer();

  // Helper Functions
  const handleClose = () => setShowPopup(false);
  const handleShow = () => setShowPopup(true);

  const handleStart = () =>{
    setShowSecondPage(true);
    makeTestCall();
  }

  const makeTestCall = async() => {
    const api_base_url='https://digital-nudge-server.onrender.com'


    try{
      const response = await fetch(`${api_base_url}/api/test`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Access-Control-Allow-Origin': '*'
        },
        body: JSON.stringify([id, group])
      });
    }
    catch {
      console.log("failed");
  }
}

  const [showInstructions, setShowInstructions] = useState(false);
  const handleShowInstructions = () => setShowInstructions(true);
  const handleCloseInstructions = () => setShowInstructions(false);




  const [showConfirmPopupAccept2, setShowConfirmPopupAccept2] = useState(false);
  const [showConfirmPopupDecline2, setShowConfirmPopupDecline2] = useState(false);

  const handleConfirmOpenAccept = () => setShowConfirmPopupAccept2(true);
  const handleConfirmCloseAccept = () => setShowConfirmPopupAccept2(false);

  const handleConfirmOpenDecline = () => setShowConfirmPopupDecline2(true);
  const handleConfirmCloseDecline = () => setShowConfirmPopupDecline2(false);


  const [demographics, setDemographics] = useState(true)
  const [showStartPage, setShowStartPage] = useState(true)
  const [showMiddlePage, setShowMiddlePage] = useState(false)
  const [endEarly, setShowEndEarly] = useState(false)
  const [showSecondPage, setShowSecondPage] = useState(false)
  const [leftTab, setLeftTab] = useState(false)
  const [trialQuestion, setShowTrialQuestion] = useState(false)
  const [useTrialPopup, setUseTrialPopup] = useState(false)
  const [useTrialPopup2, setUseTrialPopup2] = useState(false)
  const [useTrialPopup3, setUseTrialPopup3] = useState(false)
  const [useTrialPopup4, setUseTrialPopup4] = useState(false)

  const [showFirstConfirmationPopup, setFirstConfirmationPopup] = useState(false)
  const [showSecondConfirmationPopup, setSecondConfirmationPopup] = useState(false)
  const [showThirdConfirmationPopup, setThirdConfirmationPopup] = useState(false)
  const [currentTrialQuestion, setCurrentTrialQuestion] = useState(0)

  const handleCloseTrial = () => setUseTrialPopup(false);
  const handleCloseTrial2 = () => setUseTrialPopup2(false);
  const handleCloseTrial3 = () => setUseTrialPopup3(false);
  const handleCloseTrial4 = () => setUseTrialPopup4(false);

  // const GradualTextDisplay = ({ text }) => {
  //   const [displayedText, setDisplayedText] = useState('');


  const handleCloseSecond = () => {
    setCurrentTrialQuestion(currentTrialQuestion + 1)
    setFirstConfirmationPopup(false)
    setSecondConfirmationPopup(false)
    handleCloseTrial3()
    setTrialPrimed(true)

  }

const handleStartSurveyAfterTrial = () => {
    setTime(time => [...time, Date.now()]);
    setDemographics(false);
    setSelectedAnswerIndex(null); 
    setSelectedAnswerIndex2(null);
    setSelectedAnswerIndexLikert(null);
    setSelectedAnswerIndexOutside(null);
    if (primed === true){
      setPrimedExample(true)
    }
}

  const handleCloseThird = () =>  {
    setTime(time => [...time, Date.now()]);
    setSelectedAnswerIndex(null); 
    setSelectedAnswerIndex2(null);
    setSelectedAnswerIndexLikert(null);
    setSelectedAnswerIndexOutside(null);
    setShowEndTrial(true)

  }
  const handleCloseFirst = () => {
    setCurrentTrialQuestion(currentTrialQuestion + 1)
    setFirstConfirmationPopup(false)
    setTrialPrimed(true)
  }

  const trialQuestions = getTrialQuestions();


  const completionCodeAccepted = "C14XT22F"
  const completionCodeNotAccepted = "CXCBPSGM"
  const completionCodeNotStarted = "C1CK1U46"


  function shuffleArray(array) {
    for (let i = array.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [array[i], array[j]] = [array[j], array[i]];
    }
    return array
  }

  const [startTime, setStartTime] = useState([]);
  const [time, setTime] = useState([]);


  const whichPrompt = (rag) => {
    if (rag){
      return getRAGTexts();
    }
    else{
      return getPrompts();
    }
  }
  const prompts = whichPrompt(rag);

  // var rndInt = 0
  // React.useEffect(function() {
  //   if (rag){
  //     rndInt = Math.floor(Math.random() * 3)
  //   }
  //   console.log(rndInt)
  //   }, [])

  

  const order = Array.from(Array(prompts.length).keys())
  const [randomArray, setRandomArray] = React.useState([])
    React.useEffect(function() {
    setRandomArray(shuffleArray(order))
    setStartTime([Date.now()])
    // console.log(startTime)
    
    }, [])


  const likert_questions = 
    {
        code: "likert",
        text: "How trustworthy would you rate this answer?",
        options: answers_mapped["likert"]
    };

  const outside_knowledge_questions = 
  {
      code: "outside",
      text: "When answering questions related to this text, how much did you rely on previous outside knowledge?",
      options: answers_mapped["outside"]
  };

  
  // console.log(randomArray)
  
  const start_disclaimer = `I am not a doctor, but I can offer some general information that may help you learn more about the subject 
   and perhaps make an informed decision. `
  const disclaimer = `Keep in mind that guidelines and recommendations can change, and healthcare decisions 
  should be made in consultation with a qualified medical professional who is familiar with your specific health profile. ChatGPT is not 
  a medical provider and is not qualified to give medical advice. ChatGPT can make mistakes.`


  const demo_questions = getDemographicQuestions();

  const scrollToTopPrimedTrial = () => {
    window.scrollTo({ top: 0, behavior: "instant" });
    onClickNextPrimedTrial()
  }

  const onClickNextPrimedTrial = () => {
  

    setTime(time => [...time, Date.now()]);   
    setSelectedAnswerIndex(null); 
    setSelectedAnswerIndex2(null);
    setSelectedAnswerIndexLikert(null);
    setTrialPrimed(false)
      // setSelectedAnswerIndexOutside(null);


  }

  const scrollToTopPrimed = () => {
    console.log("scrolled")
    window.scrollTo({ top: 0, behavior: "instant" });
    onClickNextPrimed()
  };

  const onClickNextPrimed = () => {
  

    setTime(time => [...time, Date.now()]);
    console.log(questions[randomArray[currentQuestion]].code + "_inital", questions[randomArray[currentQuestion]].options[selectedAnswerIndex].text)
    localStorage.setItem(questions[randomArray[currentQuestion]].code + "_inital", questions[randomArray[currentQuestion]].options[selectedAnswerIndex].text);

   
    handleShow();
    setSelectedAnswerIndex(null); 
    setSelectedAnswerIndex2(null);
    setSelectedAnswerIndexLikert(null);
    setPrimedExample(false)
      // setSelectedAnswerIndexOutside(null);


  }


  const scrollToTop = () => {
    console.log("scrolled")
    window.scrollTo({ top: 0, behavior: "instant" });
    onClickNext()
  };

  // TODO: Change this func to record string of answer instead of id.
  const onClickNext = () => {
  

    setTime(time => [...time, Date.now()]);
    localStorage.setItem(questions[randomArray[currentQuestion]].code, questions[randomArray[currentQuestion]].options[selectedAnswerIndex].text);
    localStorage.setItem(questions[randomArray[currentQuestion]].code + "_likert", likert_questions.options[selectedAnswerIndexLikert].text);

    if (selectedAnswerIndex2 != null){
      localStorage.setItem(questions[randomArray[currentQuestion]].code + "_2", questions[randomArray[currentQuestion]].options2[selectedAnswerIndex2].text);
    }
    if (currentQuestion + 1 < prompts.length) {
      setCurrentQuestion(currentQuestion + 1);
 
      handleShow();
      setSelectedAnswerIndex(null); 
      setSelectedAnswerIndex2(null);
      setSelectedAnswerIndexLikert(null);
      // setSelectedAnswerIndexOutside(null);
      

    } else {
      setShowResults(true);
    }

    if (primed === true){
      setPrimedExample(true)
    }
  }

  function CheckVisible() {
    const visible = useIsTabVisible()
    if (!visible && !leftTab){
      setLeftTab(true)
    }
    // setLeftTab([...leftTab,visible])
    // console.log(leftTab)
    // console.log(visible)
  }
  CheckVisible()

  /// Add local storage logging to checkvisable at very end of survey



  const scrollToTopDemo = () => {
    console.log("scrolled_demo")
    window.scrollTo({ top: 0, behavior: "instant" });
    onClickNextDemo()
  };

  // TODO: Change this func to record string of answer instead of id.
  const onClickNextDemo = () => {
    console.log({ question_id: demo_questions[currentDemoQuestion].code, answer_key: demo_questions[currentDemoQuestion].options[selectedAnswerIndex].text });
    localStorage.setItem(demo_questions[currentDemoQuestion].code, demo_questions[currentDemoQuestion].options[selectedAnswerIndex].text);

    if (currentDemoQuestion + 1 < demo_questions.length) {
      setCurrentDemoQuestion(currentDemoQuestion + 1);
      handleShow();
      setSelectedAnswerIndex(null);
      setSelectedAnswerIndex2(null);
      setSelectedAnswerIndexLikert(null);
      // setSelectedAnswerIndexOutside(null);


    } else {
      setShowMiddlePage(true);
      setSelectedAnswerIndex(null);
      setSelectedAnswerIndex2(null);
      setSelectedAnswerIndexLikert(null);
      // setSelectedAnswerIndexOutside(null);

    }
  } 


  const showSurvey = (answer) => {
    console.log("scrolled")
    window.scrollTo({ top: 0, behavior: "instant" });
    if (answer===2){
      if (currentTrialQuestion === 0){
        setFirstConfirmationPopup(true)
      }
      else if (currentTrialQuestion === 1){
        setSecondConfirmationPopup(true)
      }
      else{
        setThirdConfirmationPopup(true)
      }
    }
    else if (answer === 4 && currentTrialQuestion ===1) {
      setUseTrialPopup3(true)
    }
    else if (answer === 0 && currentTrialQuestion ===2) {
      setUseTrialPopup4(true)
    }
    else{

      if (currentTrialQuestion === 2){
        setUseTrialPopup2(true)
      }
      else{
        setUseTrialPopup(true)
      }
    }
    setSelectedAnswerIndex(null)
    setSelectedAnswerIndex2(null);
    setSelectedAnswerIndexLikert(null);
    // setSelectedAnswerIndexOutside(null);


  }

  const showTrialQuestion = () => {
    setShowTrialQuestion(true)
    setSelectedAnswerIndex(null)
    setSelectedAnswerIndex2(null);
    setSelectedAnswerIndexLikert(null);
    // setSelectedAnswerIndexOutside(null);


  }


  const beginSurvey = () => {
    setShowStartPage(false)
    setSelectedAnswerIndex(null)
    setSelectedAnswerIndex2(null);
    setSelectedAnswerIndexLikert(null);
    // setSelectedAnswerIndexOutside(null);

  }
  const endSurvey = () => {
    setShowEndEarly(true)
    setSelectedAnswerIndex(null)
    setSelectedAnswerIndex2(null);
    setSelectedAnswerIndexLikert(null);
    // setSelectedAnswerIndexOutside(null);


  }


  class AnswerButton extends Component {
    onAnswerSelected() {
      setSelectedAnswerIndex(this.props.option.id)
    }

    render() {

      return (
        <li
          key={this.props.option.id}
          onClick={() => this.onAnswerSelected()}
          className={selectedAnswerIndex === this.props.option.id ? 'selected-answer' : null}
        >
          {this.props.option.text}
        </li>);

    }
  }

  class AnswerButtonLikert extends Component {
    onAnswerSelectedLikert() {
      setSelectedAnswerIndexLikert(this.props.option.id)
    }
    render() {

      return (
        <li
          key={this.props.option.id}
          onClick={() => this.onAnswerSelectedLikert()}
          className={selectedAnswerIndexLikert === this.props.option.id ? 'selected-answer' : null}
        >
          {this.props.option.text}
        </li>);

    }
  }

  class AnswerButtonOutside extends Component {
    onAnswerSelectedOutside() {
      setSelectedAnswerIndexOutside(this.props.option.id)
    }
    render() {

      return (
        <li
          key={this.props.option.id}
          onClick={() => this.onAnswerSelectedOutside()}
          className={selectedAnswerIndexOutside === this.props.option.id ? 'selected-answer' : null}
        >
          {this.props.option.text}
        </li>);

    }
  }

  
  class AnswerButtonSecond extends Component {
    onAnswerSelected2() {
      setSelectedAnswerIndex2(this.props.option.id)
    }

    render() {

      return (
        <li
          key={this.props.option.id}
          onClick={() => this.onAnswerSelected2()}
          className={selectedAnswerIndex2 === this.props.option.id ? 'selected-answer' : null}
        >
          {this.props.option.text}
        </li>);

    }
  }

  class MultiLineText extends Component {
    render() {
      return (
        <div>
          {this.props.text.map((i, key) => {
             if (UseChatGPTDisclaimer && 0 === key && this.props.medText) {
              // TODO: Add padding
              return <div className={this.props.className} key={key}>
              {start_disclaimer}
              {i}
              </div>;
            }
            else if (UseChatGPTDisclaimer && this.props.text.length - 1 === key && this.props.medText) {
              // TODO: Add padding
              return <div className={this.props.className} key={key}>
                {i}{disclaimer}
                </div>;
            }
            return <div className={this.props.className} key={key}>
             {i}
              </div>;
          })}
        </div>);
    }
  }


  /* Resets the game back to default */
  const restartGame = () => {
    // setScore(0);
    setCurrentQuestion(0);
    setShowResults(false);
    setSelectedAnswerIndex(null)
    handleShow();
  };

  const acceptedRestart = () => {
    // setScore(0);
    setCurrentQuestion(0);
    // setShowResults(false);
    setSelectedAnswerIndex(null)
    // handleShow();
    storeHistory();
    setShowEndScreen(true)
    setAcceptedOrNot(true)
  };

  const declinedRestart = () => {
    // setScore(0);
    setCurrentQuestion(0);
    // setShowResults(false);
    setSelectedAnswerIndex(null);
    // handleShow();
    setShowEndScreen(true);
    declinedstoreHistory();
  };


  const declinedstoreHistory = async () => {
    localStorage.clear();

    localStorage.setItem("id", id);
    localStorage.setItem("accepted", false)
    localStorage.setItem("group", group )


    const history = { ...localStorage };


    const api_base_url='https://digital-nudge-server.onrender.com'
    const response = await fetch(`${api_base_url}/api/sql`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Access-Control-Allow-Origin': '*'
      },
      body: JSON.stringify([id, history])
    });

  }

  function delay (ms) {
    return new Promise(res => setTimeout(res, ms));
  }

//   const get_code = async () => {
    
//   const api_base_url='https://digital-nudge-server.onrender.com'
//   await fetch(`${api_base_url}/api/completion`, {
//     method: 'GET',
//     headers: {
//       'Content-Type': 'application/json',
//       'Access-Control-Allow-Origin': '*'
//     },
//   }).then((response) => response.json()) //2
// .then((code) => {
// // console.log(code);
// return code; //3
// });
// console.log(response.code)


// console.log(response.json())
// const printAddress = async () => {
//   const a = await address;
//   console.log(a);
// };
// }
// printAddress();

  const RETRY_COUNT = 5;

  const storeHistory = async () => {

    localStorage.setItem("leftTab", leftTab);
    localStorage.setItem("id", id);
    localStorage.setItem("accepted", true);
    localStorage.setItem("group", group );
    localStorage.setItem("start_time",startTime);
    localStorage.setItem("question_times", time);
    localStorage.setItem("order", randomArray);



    const history = { ...localStorage };
    console.log(history);

    console.log("api call")

    const api_base_url='https://digital-nudge-server.onrender.com'

    let count = RETRY_COUNT;

    var response = null;

    try{
        response = await fetch(`${api_base_url}/api/sql`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Access-Control-Allow-Origin': '*'
        },
        body: JSON.stringify([id, history])
      });
    }

    catch {
      console.log("failed");
      var succeed = false;
      while (count > 0 && !succeed)
        await delay(500)
          {
            try{
            response = await fetch(`${api_base_url}/api/sql`, {
              method: 'POST',
              headers: {
                'Content-Type': 'application/json',
                'Access-Control-Allow-Origin': '*'
              },
              body: JSON.stringify([id, history])
            });
            succeed = true;
            }
            catch{
              console.log(count)
              count -=1
              succeed = false;
            }
          }
      if (succeed){
        console.log("succeeded")
      }
      else{
        console.log("too many retries")
      }
    }

    // if (response.status != 200){
    //   while (count > 0){
    //     response = await fetch(`${api_base_url}/api/sql`, {
    //       method: 'POST',
    //       headers: {
    //         'Content-Type': 'application/json',
    //         'Access-Control-Allow-Origin': '*'
    //       },
    //       body: JSON.stringify([id, history])
    //     });
    //     count -=1
    //   }
    // }
    // console.log(response)
    // if (count === 0){
    //   console.log("too many retries")
    // } 


    // const body = await response.text();

    localStorage.clear();

    // Push history to database


  }
  return (
    <div className="App">
      {/* 1. Header  */}
      <div className="title-box" id="title">
        <h1 className='title'> AI Recommendations Under Different Contexts Study</h1>
      </div>


      {/* Show Start page  */}

      {/* Show Questions page  */}



      {showStartPage ? (
        <div>

          {endEarly ? (
            // End Screen
            // console.log(get_code()),
            <div>
              <h2 className="page-subtitle">You have chosen to not participate. As you have indicated that you do not consent to participate in this study please return this submission on Prolific by selecting the 'stop without completing' button' or by entering the following code {completionCodeNotStarted}.</h2>
            </div>
          ) : (

            <div>
              {showSecondPage ? (
                // Explantion page
                [
                  <div className="page-box mini-box">
                    <h2 className="page-subtitle"> The following section of the survey consists of a couple demographic questions.</h2>
                    <div>
                      <button className="next-button" style={{ margin: 0 }} onClick={() => beginSurvey()}>Next</button>
                    </div>
                  </div>
                ]
              ) : (
                [
                  // Start screen
                  <div>
                    <h1 className="page-title">CONSENT SHEET</h1>
                    <h2 className="page-subtitle">University of Michigan HUM#00243219</h2>

                    <div className="page-box">
                      <div className="page-text-box">
                        {/* <GradualStringList strings={<MultiLineText className="page-body-text" text={getStartPage()} medText={false}></MultiLineText>}/> */}
                        <MultiLineText className="page-body-text" text={getStartPage()} medText={false}></MultiLineText>
                      </div>
                    </div>
                  </div>,
                  <div>
                    <button className="restart-button accept-button" onClick={() => handleStart()}>I accept</button>
                    <button className="restart-button reject-button" onClick={() => endSurvey()}>I don't accept </button>
                  </div>

                ]
              )}

            </div>


          )}

        </div>
      ) :

        (
          <div>
            {demographics ? (
              <div>
                {showMiddlePage ? (
                  <div>
                  {trialQuestion ? (


                    <div>
                    {trial_primed ? (
                      
                      
                      <div className="question-card">
                      {/* Trial Primed */}
                      <h2 className="question-text">
                        Question: {currentTrialQuestion + 1} out of {trial_questions.length}
                      </h2>

                      <h3 className="question-text">{trial_questions[currentTrialQuestion].text}</h3>

                      {/* List of possible answers  */}
                      <ul>
                        {trial_questions[currentTrialQuestion].options.map((option) => {
                          return (
                            <AnswerButton option={option} />
                          );
                        })}
                      </ul>
      
                      <PopupAlert title={"Instructions"} showPopupMode={showInstructions} closeModal={handleCloseInstructions} openModal={handleShowInstructions} text={`
                            In this section, you'll see a conversation with an AI Agent, along with questions both before and after.  You may update your previous answer or put down the same answer again. Please answer the questions to the best of your ability and outside knowledge. Please try to answer the questions as accurately as possible, based on the passage and your own outside knowledge.
                            `} ></PopupAlert>
      
                      <div className="button-box">
                            <button className="instruction-button" onClick={handleShowInstructions}> Instructions </button>
                              <button className="next-button" onClick={scrollToTopPrimedTrial} disabled={selectedAnswerIndex === null}> Next </button>
                            </div>
      
                      </div>
                      
                      ):(
                    // End of Trial Questions
                      
                    <div>
                    {showEndTrial ? (
                          [
                            <div className="page-box mini-box">
                              <h2 className="page-subtitle"> 
                              As you saw in the previous example, AI models sometimes make mistakes. Make sure to make your decision using both the conversation at hand and your own outside knowledge. If you don't know the answer, it is okay to click "Unsure".
                              Click next to begin the actual survey. Note there will be no confirmation pop up in the actual survey.</h2>
                              <div>
                                <button className="next-button" style={{ margin: 0 }} onClick={() => handleStartSurveyAfterTrial()}>Next</button>
                              </div>
                            </div>
                          ]
                          ):(
                        // Trial Question
                        [
                          [
                            <div className="Main">
                              <div className="main-box">
                                <div className={styles.gptBox}>
                                  <div className="scenario-box">
                                    <h2 className="scenario-text">Trial Question {currentTrialQuestion+1} out of 3</h2>
                                  </div>
                                  <div className="prompt-box">
                                    <div className="text-box">
                                      <h2 className="gpt-text">{trial_questions[currentTrialQuestion].prompt}</h2>
                                    </div>
                                    <div className="gpt-image"> <img src={require("./Person-bubble.jpeg")} alt="Person" style={{ maxWidth: 60 }}></img></div>
                                  </div>
                                  <div className="answer-box">
                                    <div className="medical-text-box">
                                      <div className="gpt-image"> <img src={require("./GPT-logo.jpeg")} alt="gpt Logo" style={{ maxWidth: 50 }}></img></div>
                                      {/* {UseChatGPTDisclaimer ? (<h3 className="gpt-text">{prompts[currentQuestion].text}{disclaimer} </h3>) : (<h3 className="gpt-text">{prompts[currentQuestion].text}</h3>)} */}
                                      <div> 
                                        <MultiLineText className="text-chunk" text={trial_questions[currentTrialQuestion].chatGPT_answer} medText={false} />
                                        {/* <MultiLineText className="text-chunk" text={<GradualTextDisplay text={trial_questions[currentTrialQuestion].chatGPT_answer}/>} medText={false} /> */}

                                      </div>
                                    </div>
                                  </div>
                                </div>
                              </div>
                            </div>,
        
                            /* 5. Question Card  */
                            <div className="question-card">
                              {/* Current Question  */}
                              <h2 className="question-text">
                                Question: {currentTrialQuestion + 1} out of {trial_questions.length}
                              </h2>
                              {/* EDIT HERE TO ADD PARARGRAPHS FOR EACH STR IN ARR */}
                              <h3 className="question-text">{trial_questions[currentTrialQuestion].text}</h3>
        
                              {/* List of possible answers  */}
                              <ul>
                                {trial_questions[currentTrialQuestion].options.map((option) => {
                                  return (
                                    <AnswerButton option={option} />
                                  );
                                })}
                              </ul>
        
                              <h3 className="question-text">{"How would you rate the trustworthiness of the above passage?"}</h3>
        
                              <ul>
                                {likert_questions.options.map((option) => {
                                  return (
                                    <AnswerButtonLikert id={"trust"} option={option} />
                                  );
                                })}
                              </ul>
        
                              <div className="button-box">
                              <button className="instruction-button" onClick={handleShowInstructions}> Instructions </button>
                                <button className="next-button" onClick={() => showSurvey(selectedAnswerIndex)} disabled={selectedAnswerIndex === null || selectedAnswerIndexLikert === null}> Next </button>
                              </div>
        
                              <PopupAlert title={"Instructions"} showPopupMode={showInstructions} closeModal={handleCloseInstructions} openModal={handleShowInstructions} text={`
                              In this section, you'll see a conversation with an AI Agent, followed by some questions. 
                              Please try to answer the questions as accurately as possible, based on the passage and your own outside knowledge.
                              `} ></PopupAlert>
                            
                              <PopupAlert showPopupMode={useTrialPopup} closeModal={handleCloseTrial} openModal={handleShow} text={"Hmm, that doesn't seem quite right. Please try again. Select Unsure if you don't know the answer."}></PopupAlert>
                              <PopupAlert showPopupMode={useTrialPopup2} closeModal={handleCloseTrial2} openModal={handleShow} text={"Even though the conversation stated that a fever begins at 104 degrees Fahrenheit, as an AI model, it sometimes gets things wrong. Make sure to make your decision using both the conversation at hand and your own outside knowledge."}></PopupAlert>
                              <PopupAlert showPopupMode={useTrialPopup3} closeModal={handleCloseSecond} openModal={handleShow} text={"It's okay to be unsure! In this study, please select `Unsure` if you don't know the answer."}></PopupAlert>
                              <PopupAlert showPopupMode={useTrialPopup4} closeModal={handleCloseThird} openModal={handleShow} text={"It's okay to be unsure! In this study, please answer Unsure if you don't know the answer. As you may have noticed, AI models sometimes make mistakes. Make sure to make your decision using both the conversation at hand and outside knowledge. Click continue to begin the survey."}></PopupAlert>
        
                                
                              <PopupAlert showPopupMode={showFirstConfirmationPopup} closeModal={handleCloseFirst} openModal={handleShow} text={"Nice Job. This next example will ask you a question before showing you the full conversation. Select Unsure if you don't know the answer."} ></PopupAlert>
          
                              <PopupAlert showPopupMode={showSecondConfirmationPopup} closeModal={handleCloseSecond} openModal={handleShow} text={"Nice Job. There's one more trial question. As a reminder, if you don't know the answer, select Unsure."} ></PopupAlert>
        
                              <PopupAlert showPopupMode={showThirdConfirmationPopup} closeModal={handleCloseThird} openModal={handleShow} text={`Nice Job. As you saw in the previous example, AI models sometimes make mistakes. Make sure to make your decision using both the conversation at hand and your own outside knowledge.`} ></PopupAlert>
        
                            </div>
                           
                          ]
                        ])

                        }
                        </div> 
                  )

                }
                </div>

                  ) : (
                    // Finished Demographics
                  <div className="page-box mini-box">

                    <h2 className="page-subtitle">Thank you for filling out the demographic portion of the survey. 
                    In the following section, you'll see a conversation with an AI Agent, along with questions both before and after.
                    Please try to answer the questions as accurately as possible, based on the passage and your own outside knowledge. 
                    There will be 3 practice questions to help you get acclimated.</h2>
                    <div>
                      <button className="next-button" style={{ margin: 0 }} onClick={showTrialQuestion} > Next </button>
                    </div>
                  </div>
                  
                  )} 

                  </div>
                
                ) : (
                  // Demographic Questions

                  <div className="question-card">
                    {/* Current Question  */}
                    <h2 className="question-text">
                      Question: {currentDemoQuestion + 1} out of {demo_questions.length}
                    </h2>
                    {/* EDIT HERE TO ADD PARARGRAPHS FOR EACH STR IN ARR */}
                    <h3 className="question-text">{demo_questions[currentDemoQuestion].text}</h3>

                    {/* List of possible answers  */}
                    <ul>
                      {demo_questions[currentDemoQuestion].options.map((option) => {
                        return (
                          <AnswerButton option={option} />
                        );
                      })}
                    </ul>
                    <div className="button-box">
                      
                      <button className="next-button" onClick={scrollToTopDemo} disabled={selectedAnswerIndex === null}> Next </button>
                    </div>
                  </div>
                )}
              </div>

            ) : (

              <div>
                {showResults ? (
                  // storeHistory(),
                  /* 3. Final Results */
                  <div /* className="final-results"*/ >
                    {showEndScreen ? (
                      <div className="page-box">
                        <h1 className="page-title">Thank you for taking the survey!</h1>
                        {acceptedOrNot ? (<h2 className="page-subtitle"> You have chosen to accept the terms and your response has been recorded. Your completion code is {completionCodeAccepted}</h2>
                        ) : (<h2 className="page-subtitle">You have chosen to not accept the terms. Your answers have not been recorded. Your completion code is {completionCodeNotAccepted}</h2>)}
                      </div>
                    ) : (

                      <div>
                        <div className="page-box">
                        <h1 className="page-title">POST SURVEY DISCLAIMER SHEET</h1>
                         <h2 className="page-subtitle">University of Michigan HUM#00243219</h2>

                          <h2 className="page-subtitle">
                            Please read the following in its entirety. 
                          </h2>
                          <div className="page-text-box">
                            <MultiLineText className="page-body-text" text={post_disclaimer} medText={false}></MultiLineText>
                          </div>
                        </div>
                        <button className="restart-button accept-button" onClick={handleConfirmOpenAccept}>I accept</button>
                        <button className="restart-button reject-button" onClick={handleConfirmOpenDecline}>I don't accept </button>

                        <PopupConfirm showPopupMode={showConfirmPopupAccept2} acceptModal={acceptedRestart} closeModal={handleConfirmCloseAccept} text={"You have clicked to particapte in this study. Is this correct?"}></PopupConfirm>
                        <PopupConfirm showPopupMode={showConfirmPopupDecline2} acceptModal={declinedRestart} closeModal={handleConfirmCloseDecline} text={"You have clicked to NOT particapte in this study. Is this correct?"}></PopupConfirm>

                      </div>
                    )}

                  </div>
                ) : (

                <div>

                
                {primed_example ?  (
                  //New Primed questions
                <div className="question-card">
                {/* Current Question  */}
                <h2 className="question-text">
                  Question: {currentQuestion + 1} out of {prompts.length}
                  {console.log(primed_example)}

                </h2>


                {/* Question 1 */}
                <h3 className="question-text">{questions[randomArray[currentQuestion]].text}</h3>

                {/* List of possible answers  */}
                <ul>
                  {questions[randomArray[currentQuestion]].options.map((option) => {
                    return (
                      <AnswerButton id={"comp"} option={option} />
                    );
                  })}
                </ul>

                <PopupAlert title={"Instructions"} showPopupMode={showInstructions} closeModal={handleCloseInstructions} openModal={handleShowInstructions} text={
                  `In this section, you'll see a conversation with an AI Agent, 
                  along with questions both before and after.  
                  You may update your previous answer or put down the same answer again.
                   Please answer the questions to the best of your ability and outside knowledge.
                    Please try to answer the questions as accurately as possible, based on the passage 
                    and your own outside knowledge.
                      .`} ></PopupAlert>

                <div className="button-box">
                      <button className="instruction-button" onClick={handleShowInstructions}> Instructions </button>
                        <button className="next-button" onClick={scrollToTopPrimed} disabled={selectedAnswerIndex === null}> Next </button>
                      </div>

                </div>)
                :

                  ( //new
                  [
                    <div className="Main">
                      <div className="main-box">
                        <div className={styles.gptBox}>
                          <div className="scenario-box">
                            <h2 className="scenario-text">Scenario {currentQuestion + 1} out of {prompts.length}:</h2>
                          </div>
                          <div className="prompt-box">
                            <div className="text-box">
                              <h2 className="gpt-text">{prompts[randomArray[currentQuestion]].prompt}</h2>
                            </div>
                            <div className="gpt-image"> <img src={require("./Person-bubble.jpeg")} alt="Person" style={{ maxWidth: 60 }}></img></div>
                          </div>
                          <div className="answer-box">
                            <div className="medical-text-box">
                              <div className="gpt-image"> <img src={require("./GPT-logo.jpeg")} alt="gpt Logo" style={{ maxWidth: 50 }}></img></div>
                              {/* {UseChatGPTDisclaimer ? (<h3 className="gpt-text">{prompts[currentQuestion].text}{disclaimer} </h3>) : (<h3 className="gpt-text">{prompts[currentQuestion].text}</h3>)} */}
                              <div>
                                <MultiLineText className="text-chunk" text={prompts[randomArray[currentQuestion]].text} medText={true} />
                              </div>
                            </div>
                          </div>
                        </div>
                      </div>
                    </div>,

                    /* 5. Question Card  */
                    <div className="question-card">
                      {/* Current Question  */}
                      <h2 className="question-text">
                        {/* Question: {currentQuestion + 1} out of {prompts.length} */}
                      </h2>


                      {/* Question 1 */}
                      <h3 className="question-text">{questions[randomArray[currentQuestion]].text}</h3>

                      {/* List of possible answers  */}
                      <ul>
                        {questions[randomArray[currentQuestion]].options.map((option) => {
                          return (
                            <AnswerButton id={"comp"} option={option} />
                          );
                        })}
                      </ul>
                  
                      {/* Question 2 */}
                      <h3 className="question-text">{questions[randomArray[currentQuestion]].text2}</h3>
                      <ul>
                        {questions[randomArray[currentQuestion]].options2.map((option) => {
                          return (
                            <AnswerButtonSecond id={"comp"} option={option} />
                          );
                        })}
                      </ul>

                      {/* <h3 className="question-text">{"To what extent did you rely on prior knowledge when answering the above questions?"}</h3>

                      <ul>
                        {outside_knowledge_questions.options.map((option) => {
                          return (
                            <AnswerButtonOutside id={"outside"} option={option} />
                          );
                        })}
                      </ul> */}
                      <h3 className="question-text">{"How would you rate the trustworthiness of the above passage?"}</h3>

          
                      <ul>
                        {likert_questions.options.map((option) => {
                          return (
                            <AnswerButtonLikert id={"trust"} option={option} />
                          );
                        })}
                      </ul>
                      
                      <PopupAlert title={"Instructions"} showPopupMode={showInstructions} closeModal={handleCloseInstructions} openModal={handleShowInstructions} text={`
                      In this section, you'll see a conversation with an AI Agent, followed by some questions. 
                      Please try to answer the questions as accurately as possible, based on the passage and your own outside knowledge.
                      `} ></PopupAlert>

                      <div className="button-box">
                      <button className="instruction-button" onClick={handleShowInstructions}> Instructions </button>
                        <button className="next-button" onClick={scrollToTop} disabled={selectedAnswerIndex === null || selectedAnswerIndex2 === null || selectedAnswerIndexLikert === null }> Next </button>
                      </div>
                    </div>

                  ]

                ) //new
                  
                 } 
                 </div>
                 )
                  
                }

                {UsePopups ?
                  (<PopupAlert title={"Alert"}  showPopupMode={showPopups} closeModal={handleClose} openModal={handleShow} text={"ChatGPT sometimes writes plausible-sounding but incorrect answers. Does this information seem accurate?"} />
                  ) : (null)}
              </div>



            )}

          </div>)}




    </div>

  );
}

export default App;
