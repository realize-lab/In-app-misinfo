// const getMainBox = () => {
//     const box = 
//     <>
// [
//     <div className="Main">
//       <div className="main-box">
//         <div className={styles.gptBox}>
//           <div className="scenario-box">
//             <h2 className="scenario-text">Scenario {currentQuestion + 1} out of {prompts.length}:</h2>
//           </div>
//           <div className="prompt-box">
//             <div className="text-box">
//               <h2 className="gpt-text">{prompts[randomArray[currentQuestion]].prompt}</h2>
//             </div>
//             <div className="gpt-image"> <img src={require("./Person-bubble.jpeg")} alt="Person" style={{ maxWidth: 60 }}></img></div>
//           </div>
//           <div className="answer-box">
//             <div className="medical-text-box">
//               <div className="gpt-image"> <img src={require("./GPT-logo.jpeg")} alt="gpt Logo" style={{ maxWidth: 50 }}></img></div>
//               {/* {UseChatGPTDisclaimer ? (<h3 className="gpt-text">{prompts[currentQuestion].text}{disclaimer} </h3>) : (<h3 className="gpt-text">{prompts[currentQuestion].text}</h3>)} */}
//               <div>
//                 <MultiLineText className="text-chunk" text={prompts[randomArray[currentQuestion]].text} medText={true} />
//               </div>
//             </div>
//           </div>
//         </div>
//       </div>
//     </div>,

//     /* 5. Question Card  */
//     <div className="question-card">
//       {/* Current Question  */}
//       <h2 className="question-text">
//         Question: {currentQuestion + 1} out of {prompts.length}
//       </h2>


//       {/* Question 1 */}
//       <h3 className="question-text">{questions[randomArray[currentQuestion]].text}</h3>

//       {/* List of possible answers  */}
//       <ul>
//         {questions[randomArray[currentQuestion]].options.map((option) => {
//           return (
//             <AnswerButton id={"comp"} option={option} />
//           );
//         })}
//       </ul>
  
//       {/* Question 2 */}
//       <h3 className="question-text">{questions[randomArray[currentQuestion]].text2}</h3>
//       <ul>
//         {questions[randomArray[currentQuestion]].options2.map((option) => {
//           return (
//             <AnswerButtonSecond id={"comp"} option={option} />
//           );
//         })}
//       </ul>

//       {/* <h3 className="question-text">{"To what extent did you rely on prior knowledge when answering the above questions?"}</h3>

//       <ul>
//         {outside_knowledge_questions.options.map((option) => {
//           return (
//             <AnswerButtonOutside id={"outside"} option={option} />
//           );
//         })}
//       </ul> */}
//       <h3 className="question-text">{"How would you rate the trustworthyness of the above passage?"}</h3>


//       <ul>
//         {likert_questions.options.map((option) => {
//           return (
//             <AnswerButtonLikert id={"trust"} option={option} />
//           );
//         })}
//       </ul>
      
//       <PopupAlert title={"Instructions"} showPopupMode={showInstructions} closeModal={handleCloseInstructions} openModal={handleShowInstructions} text={`
//       In this section, you'll see conversation with an AI Agent,
//       followed by some questions. The conversation shown contains AI generated information. 
//       Please answer the questions to the best of your ability and outside knowledge.`} ></PopupAlert>

//       <div className="button-box">
//       <button className="instruction-button" onClick={handleShowInstructions}> Instructions </button>
//         <button className="next-button" onClick={scrollToTop} disabled={selectedAnswerIndex === null || selectedAnswerIndex2 === null || selectedAnswerIndexLikert === null }> Next </button>
//       </div>
//     </div>
//   ]
//   </>

//   return box;
// }

// export default getMainBox;
