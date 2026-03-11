import React, { Component } from 'react';
import styles from './modestyles.css'
import getAlertMessage from './alertMessage';

// render text with alert message at bottom
class Mode2 extends Component {
    static alertMessage = getAlertMessage()
    constructor() {
        super();
    }
    render() {
        return (
            <div>
                <h3 className={styles.medicalText}>{super.prompts[super.currentQuestion].text}</h3> 
                <h3 className={styles.medicalText}>{this.alertMessage}</h3>
            </div>
        );
    }
}
 
export default Mode2;