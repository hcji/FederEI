# FederEI
## create with docker


## create with conda
We recommend using conda to create an environment, then installing the packages using pip.
```
cd FederEI  
conda create -n FederEI -y  
conda activate FederEI  
conda install python=3.8 -y
```
Then,install the dependencies using pip:
```
pip install -r env/requirements.txt
```
### Center server run:
```
python app/master.py 
```
### Client server run:
```
python app/server.py center_server_addr center_server_message_port center_server_file_port database_path
```
for example:
```
python app/server.py "127.0.0.1" 5000 5001 "database.mgf"
```
