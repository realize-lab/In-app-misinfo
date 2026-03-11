import React, { Component } from 'react';
import styles from './modestyles.css'

// render just text
class Mode1 extends Component {
    constructor() {
        super();
    }
    render() {
        return (
            <div>
                <h3 className={styles.medicalText}>{super.prompts[super.currentQuestion].text}</h3> 
            </div>
        );
    }
}
 
export default Mode1;