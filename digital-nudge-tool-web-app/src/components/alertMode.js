import React, { Component } from 'react';
import Mode1 from './Mode1'
import Mode2 from './Mode2'
// import Mode3 from './Mode3'
 
// Creating enum object to store modes
const Enumobj = {
    Mode1: <Mode1 />,
    Mode2: <Mode2 />,
    // Mode3: <Mode3 />
};
 
// Creating enum function to return
// Particular component according to state value
function alertMode({ state }) {
    return <div>{Enumobj[state]}</div>;
}

export default alertMode
