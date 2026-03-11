import App from './App';

import { Route, Routes, useParams} from "react-router-dom"



function Main() {
    // console.log(useParams())


    return (
        <Routes>
        <Route path="/1/" element={<App group={1} chatgpt={false} popup={false} rag={false} priming={false} />} />
        <Route path="/2/" element={<App group={2} chatgpt={true} popup={false} rag={false} priming={false} />} />  
        <Route path="/3/" element={<App group={3} chatgpt={false} popup={true} rag={false}  priming={false}/>} />
        <Route path="/4/" element={<App group={4} chatgpt={false} popup={false} rag={true}  priming={false}/>}/>
        <Route path="/5/" element={<App group={5} chatgpt={false} popup={false} rag={false}  priming={true}/>}/>
      </Routes>
      );
  }

  export default Main;