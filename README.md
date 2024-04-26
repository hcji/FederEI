# FederEI
## create env with docker

### center 
open two port
```
sudo iptables -A INPUT -p tcp --dport port1 -j ACCEPT
sudo iptables -A INPUT -p tcp --dport port2 -j ACCEPT
```
check the ip address
```
ifconfig
```
pull the docker image
```
docker pull kpbl1/federei:v1
```
run image and center server
```
docker run -i -p port1:port1 -p port2:port2 federei:v1
python3.8 /app/center.py yourip port1 port2
```

### client
client should be run after center

open two port
```
sudo iptables -A INPUT -p tcp --dport port1 -j ACCEPT
sudo iptables -A INPUT -p tcp --dport port2 -j ACCEPT
```
check the ip address
```
ifconfig
```
pull the docker image
```
docker pull kpbl1/federei:v1
```
run image 
```
docker run -d -p port1:port1 -p port2:port2 federei:v1
```
check container id
```
docker ps
```
copy database into container
```
docker cp your_database_path container_id:/app
```

```
docker attach container_id
```
run client server
```
python3.8 /app/client.py center_address center_port1 center_port2 client_address client_port1 client_port2
```

## create env with pip
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
python app/center.py 
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
python app/client.py '127.0.0.1' 5000 5001 '127.0.0.1' 6000 6001 "database.pickle"
```
