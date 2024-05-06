# Installation and Deployment

## Deployment of central-server
Pull the Docker image, which contains the environment and code required for execution:

```
docker pull kpbl1/federei:v1
```

Retrieve the IP address of the current machine for later use:

For Linux:
```
ifconfig
```

For Windows:
```
ipconfig
```

Open two TCP ports on the machine for listening, for example, ports 5000 and 5001. 
Port 1 is designated for client registration, while port 2 is for users to send files.

If your machine's firewall is managed by iptables, execute:

```
sudo iptables -A INPUT -p tcp --dport port1 -j ACCEPT
sudo iptables -A INPUT -p tcp --dport port2 -j ACCEPT
```

If the firewall is managed by firewalld, execute:

```
sudo firewall-cmd --zone=public --add-port=port1/tcp --permanent
sudo firewall-cmd --zone=public --add-port=port2/tcp --permanent
sudo firewall-cmd --reload
```

Run the central server:

```
docker run -i -p port1:port1 -p port2:port2 federei:v1
python3.8 /app/center.py your_ip_address port1 port2
```

## Deployment of client-server
The client server should be launched after the central server. You'll need the central IP and 
port information to start it. Once the client server is running, it will be automatically added 
to the registry table.

Pull the Docker image containing the environment and code:

```
docker pull kpbl1/federei:v1
```

Retrieve the IP address of the current machine for later use:

For Linux:

```
ifconfig
```

For Windows:

```
ipconfig
```

Open two TCP ports on the machine for listening, for example, ports 5000 and 5001.

If your machine's firewall is managed by iptables, execute:

```
sudo iptables -A INPUT -p tcp --dport port1 -j ACCEPT
sudo iptables -A INPUT -p tcp --dport port2 -j ACCEPT
```

If the firewall is managed by firewalld, execute:

```
sudo firewall-cmd --zone=public --add-port=port1/tcp --permanent
sudo firewall-cmd --zone=public --add-port=port2/tcp --permanent
sudo firewall-cmd --reload
```

Run the client:
Launch the container in the background and copy the database into it:

```
docker run -d -p port1:port1 -p port2:port2 federei:v1
```

Identify the ID of the container that was started:

```
docker ps
```

Then copy the database into it:

```
docker cp your_database_path container_id:/app
```

Attach to this container and run the server. The database_name is the last item of your_database_path:

```
docker attach container_id
python3.8 /app/client.py central_ip center_port1 center_port2 client_ip client_port1 client_port2 database_name
```

## Installation of user-interface

Download the software from the following link:

https://github.com/wustjie/wustjie.github.io/releases/tag/FederEI