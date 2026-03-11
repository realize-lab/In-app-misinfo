// import Modal from "react-bootstrap/Modal";

import { useState } from 'react';
import Button from 'react-bootstrap/Button';
import Modal from 'react-bootstrap/Modal';
import styles from "./popupstyles.css"

function PopupAlert({title, showPopupMode, closeModal, openModal, text}) {


  return (
    <>


      <Modal
        show={showPopupMode}
        onHide={closeModal}
        backdrop="static"
        centered={true}
        keyboard={false}
        dialogClassName="popup-modal"
      >
        <Modal.Header closeButton>
          <Modal.Title style={{color: 'black'}}>{title}</Modal.Title>
        </Modal.Header>
        <Modal.Body>
          <p style={{color: 'black'}}>{text}</p>
        </Modal.Body>
        <Modal.Footer>
          <Button variant="secondary" onClick={closeModal}>
            Continue
          </Button>
          {/* <Button variant="primary">Understood</Button> */}
        </Modal.Footer>
      </Modal>
    </>
  );
}

export default PopupAlert;