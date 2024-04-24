# FederEI
## create env with docker
```
docker pull FederEI
```

```
docker run -p 5000:5000 -p 5001:5001 my-container
```


## create env with conda
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
export port
```
sudo iptables -A INPUT -p tcp --dport 5000 -j ACCEPT
sudo iptables -A INPUT -p tcp --dport 5001 -j ACCEPT
```
```
python app/master.py 
```
### Client server run:
export port
```
sudo iptables -A INPUT -p tcp --dport 6000 -j ACCEPT
sudo iptables -A INPUT -p tcp --dport 6001 -j ACCEPT
```
get ip
```
ifconfig
```

```
python app/server.py center_addr center_message_port center_file_port Client_addr Client_message_port Client_file_port database_path
```
for example:
```
python app/server.py '127.0.0.1' 5000 5001 '127.0.0.1' 6000 6001 "database.mgf"
```
