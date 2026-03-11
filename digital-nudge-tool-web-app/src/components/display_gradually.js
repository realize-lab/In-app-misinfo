import React, { useState, useEffect } from 'react';

const GradualTextDisplay = ({ text2 }) => {
  const [displayedText, setDisplayedText] = useState('');

  useEffect(() => {
    let index = 0;
    const intervalId = setInterval(() => {
      if (index <= text2.length) {
        setDisplayedText((prevText) => prevText + text2.charAt(index));
        // console.log(displayedText);
        index++;
      } 
    //   else {
    //     clearInterval(intervalId);
    //   }
    }, 10); // Adjust the delay as needed

    // Cleanup the interval on component unmount
    return () => clearInterval(intervalId);
  }, [text2]);

  return (
    <div>
      <p>{displayedText}</p>
    </div>
  );
};

export default GradualTextDisplay;
