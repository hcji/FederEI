# Advanced Usage

## Adding a Client-Server

To add the Client-Server,please follow the steps outlined in the Installation document after starting the Center-Server. Once the Client is launched, it will be automatically registered in the Center-Server's registry table and begin functioning.

## Changing the Library

To change the library, follow these steps:

If it's already running, you'll need to stop it first and detach the container.
```
docker attach container_id
crtl+c
crtl+p,crtl+q
```

Then, copy the library into the container and restart it.

```
docker cp your_library_path container_id:/app
docker attach container_id
python3.8 /app/client.py central_ip center_port1 center_port2 client_ip client_port1 client_port2 database_name
```

### Library format
Library data supports *.pickle* and *.mgf* types.The program will automatically generate the vectors required for searching based on the provided Library. The following information for each compound in the Library is required when returning results,ensure that each compound in the Library contains at least this information: 

    SMILES
    Compound Name

If you want to modify the information returned in the search results, you can follow the format below to modify the information returned in Client.py for mgf results:

    mgf_data['compound_name']=sorted_cleaned_pkl_data[I[index1][index2]].metadata['compound_name']
    mgf_data['smiles']=sorted_cleaned_pkl_data[I[index1][index2]].metadata['smiles']

