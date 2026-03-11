const getEndScreen = () => {

const get_code = async () => {
    
    const api_base_url='https://digital-nudge-server.onrender.com'
    const response = fetch(`${api_base_url}/api/completion`, {
      method: 'GET',
      headers: {
        'Content-Type': 'application/json',
        'Access-Control-Allow-Origin': '*'
      },
    }).then((response) => response.json()) //2
  .then((code) => {
  console.log(code); //3
  });

  console.log(code)
  endScreen = <h2>You have chosen to accept the terms and your response has been recorded. Your completion code is {completionCodeAccepted}</h2>

  return endScreen;
    }
}

export default endScreen;
