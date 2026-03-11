import React, { useState, useEffect } from 'react';

const GradualStringMap = ({ stringMap, interval = 100 }) => {
  const [displayedText, setDisplayedText] = useState('');

  useEffect(() => {
    let currentKeyIndex = 0;
    let currentKey = '';
    let currentValue = '';

    const keys = Object.keys(stringMap);

    const intervalId = setInterval(() => {
      if (currentKeyIndex < keys.length) {
        currentKey = keys[currentKeyIndex];
        currentValue = stringMap[currentKey];
        const newText = `${currentKey}: ${currentValue}`;

        setDisplayedText((prevText) => prevText + newText.charAt(prevText.length));
      } else {
        clearInterval(intervalId);
      }

      if (displayedText === `${currentKey}: ${currentValue}`) {
        currentKeyIndex++;
      }
    }, interval);

    return () => clearInterval(intervalId); // Cleanup the interval on component unmount

  }, [stringMap, interval, displayedText]);

  return (
    <div>
      <p>{displayedText}</p>
    </div>
  );
};

export default GradualStringMap;