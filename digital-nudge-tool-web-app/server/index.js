// server/index.js
const sqlite3 = require('sqlite3').verbose();
const password = encodeURIComponent("nOrGl73icv0Aea3O");
const username = encodeURIComponent("realizelab");
var bodyParser = require('body-parser')
const cors = require("cors");


collec = "pilot";
collec2 = "back_fill"
collec3 = "test"
collec4 = "Increased_Sample"

const { MongoClient, ServerApiVersion } = require('mongodb');
let uri =
  `mongodb+srv://${username}:${password}@cluster0.8lo7e.mongodb.net/?retryWrites=true&w=majority`;

const client = new MongoClient(uri, {
  serverApi: {
    version: ServerApiVersion.v1,
    strict: true,
    deprecationErrors: true,
  }
});


async function update(client, id, answers){

    try {
        // Connect to the MongoDB cluster
        await client.connect();

        // Make the appropriate DB calls

        // Create a single new listing

        const new_user = {
            user_id: id ,
            answers: answers,
        }
        const result = await client.db("aitool").collection(collec4).replaceOne(
            {"user_id": id},
            {"user_id": id, "answers": answers},
            {upsert: true});
        console.log(`New listing created with the following id: ${id}`);
        
    } 
    finally {
        // Close the connection to the MongoDB cluster
        await client.close();
        console.log("connection closed")
    }
    }

  async function update_initial(client, id, answers){

      try {
          // Connect to the MongoDB cluster
          await client.connect();
  
          // Make the appropriate DB calls
  
          // Create a single new listing
  
          const new_user = {
              user_id: id ,
              answers: answers,
          }
          const result = await client.db("aitool").collection("test").replaceOne(
              {"user_id": id},
              {"user_id": id, "answers": answers},
              {upsert: true});
          console.log(`New initial test for: ${id}`);
          
      } 
      finally {
          // Close the connection to the MongoDB cluster
          await client.close();
          console.log("connection closed")
      }
      }


const express = require("express");


const PORT = process.env.PORT || 3001;

const app = express();
app.use(cors());
app.use(bodyParser.json())

app.use(express.static('build'));


app.post('/api/sql', (req, res) => {
    console.log(req.url);
    update(client, req.body[0],req.body[1])
    res.send("logged")
});

app.get('/api/completion', (req, res) => {
  // console.log(req.url);
  res.status(200).send({code:"C14XT22F"})
});


app.post('/api/test', (req, res) => {
  // console.log(req.url);
  console.log(req.url);
  update_initial(client, req.body[0],req.body[1])
  res.send("logged")
});

app.listen(PORT, () => {
  console.log(`Server listening on ${PORT}`);
});
