// import Modal from "react-bootstrap/Modal";

import { useState } from 'react';
import Button from 'react-bootstrap/Button';
import Modal from 'react-bootstrap/Modal';
import styles from "./popupstyles.css"

function PopupConfirm({showPopupMode, acceptModal, closeModal, text}) {


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
          <Modal.Title style={{color: 'black'}}>Alert</Modal.Title>
        </Modal.Header>
        <Modal.Body>
          <p style={{color: 'black'}}>{text}</p>
        </Modal.Body>
        <Modal.Footer>
          <Button variant="secondary" onClick={closeModal}>
            Go Back
          </Button>
          <Button variant="secondary" onClick={acceptModal}>
            Confirm
          </Button>
          {/* <Button variant="primary">Understood</Button> */}
        </Modal.Footer>
      </Modal>
    </>
  );
}

export default PopupConfirm;