import socket
import json
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import rsa, padding
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
from cryptography.hazmat.backends import default_backend
import os
from multiprocessing import Process, Manager
import sys
import struct
import threading
import time
import fcntl
from data_process import spec
from matchms.exporting import save_as_mgf
from rdkit.DataStructs import TanimotoSimilarity
from rdkit import Chem
from rdkit.Chem import AllChem

server_ip=sys.argv[1]
server_port=int(sys.argv[2])
file_port=int(sys.argv[3])



private_key = rsa.generate_private_key(
    public_exponent=65537,
    key_size=2048,
    backend=default_backend()
)

public_key = private_key.public_key()


public_pem = public_key.public_bytes(
    encoding=serialization.Encoding.PEM,
    format=serialization.PublicFormat.SubjectPublicKeyInfo
)

lock=threading.Lock()


def send_public_recv_symkey(client_socket,public_pem):
    client_socket.send(public_pem)
    encrypted_sym_key = client_socket.recv(256)
    sym_key = private_key.decrypt(
        encrypted_sym_key,
        padding.OAEP(
            mgf=padding.MGF1(algorithm=hashes.SHA256()),
            algorithm=hashes.SHA256(),
            label=None
        )
    )
    return sym_key

def recv_public_send_symkey(client_socket):
    public_pem = client_socket.recv(1024)
    public_key = serialization.load_pem_public_key(public_pem, backend=default_backend())

    sym_key = os.urandom(32)
    encrypted_sym_key = public_key.encrypt(
        sym_key,
        padding.OAEP(
            mgf=padding.MGF1(algorithm=hashes.SHA256()),
            algorithm=hashes.SHA256(),
            label=None
        )
    )
    client_socket.send(encrypted_sym_key)
    return sym_key

def message(server_file_addr,server_message_addr,file_email): 
    server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server.bind((server_ip, server_port))
    server.listen(10)

    while True:
        client_socket, client_address = server.accept()
        print(f"Connection from {client_address}")

        try:
            sym_key=send_public_recv_symkey(client_socket,public_pem)

            received_data = receive_json(client_socket, sym_key)

            if received_data['identity']=='client' :
                
                print("accept client message")
            elif received_data['identity']=='server':
                fileaddr={}
                messaddr={}
                fileaddr[received_data['file_addr']]=received_data['file_port']
                messaddr[received_data['message_addr']]=received_data['message_port']
                if fileaddr in server_file_addr :
                    print('file addr already exist')
                elif messaddr in server_message_addr :
                    print('message addr already exist')
                else:
                    server_file_addr.append(fileaddr)
                    server_message_addr.append(messaddr)
                    print(f"regist success,file addr {received_data['file_addr']},port {received_data['file_port']},message addr{received_data['message_addr']},port {received_data['message_port']}")
            else:
                print('illegal messege type')

        finally:
            client_socket.close()




def socket_service(server_file_addr,server_message_addr,file_email):
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        s.bind((server_ip, file_port))

        s.listen(10)
    except socket.error as msg:
        print (msg)
        sys.exit(1)
    print ('Waiting connection...')
 
    while 1:
        client_socket, client_address = s.accept()

        t = threading.Thread(target=deal_data, args=(client_socket, client_address,server_file_addr,server_message_addr,file_email))
        t.start()

def jiegouxiangsi(smiles1,smiles2): 
    mol1 = Chem.MolFromSmiles(smiles1) 
    mol2 = Chem.MolFromSmiles(smiles2) 
 
    # 生成分子指纹 
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048) 
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048) 
 
    # 计算 Tanimoto 系数 
    tanimoto_score = TanimotoSimilarity(fp1, fp2) 
    return tanimoto_score         
 
def deal_data(client_socket, client_address,server_file_addr,server_message_addr,file_email):

    print(f"Connection from {client_address}")
    sym_key=send_public_recv_symkey(client_socket,public_pem)
    file_name, file_data = receive_file(client_socket, sym_key)
    if not os.path.exists(file_name) or file_email[os.path.splitext(file_name)[0]]==0:
        with open(file_name, "wb") as f:
            f.write(file_data)

    else:
        with open(file_name, "ab") as f:
            f.write(file_data)

    lock.acquire()
    name=os.path.splitext(file_name)[0] 
    if name in file_email :
        file_email[name]+=1
    else:
        file_email[name]=0
    print(f"get file {file_name},now num of {name} is {file_email[name]}")
    lock.release()

    if file_email[name]==0 :
        for i in server_file_addr:
            for address,port in i.items() :
                #socket_client(str(fn),address,port)
                try:
                    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                    s.connect((address, port))
                except socket.error as msg:
                    print (msg)
                    sys.exit(1)
                sym_key_temp=recv_public_send_symkey(s)
                send_file(s, sym_key_temp, file_name)
                s.close()

        while 1:
            if file_email[name] == len(server_file_addr) :
                break
            else:
                time.sleep(1)
        
        result = list(spec.load_from_mgf(file_name))
        
        top_scores = {}
        for data in result:
            id = data.metadata["origin_name"]
            if id not in top_scores:
                top_scores[id] = []
            top_scores[id].append(data)
        sorted_result=[]
        for id,data_list in top_scores.items():
            data_list=sorted(data_list, key=lambda x: x.metadata['distance'])
            data_list=data_list[0:min(10,len(data_list))]
            sorted_result+=data_list

        save_as_mgf(sorted_result,file_name)
        
        
        send_file(client_socket, sym_key, file_name)

        try :
            os.remove(file_name)
        except:
            print("remove cache fail")


def socket_client(filepath,address,port):
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((address, port))
    except socket.error as msg:
        print (msg)
        sys.exit(1)



    if os.path.isfile(filepath):

        fileinfo_size = struct.calcsize('128sl')

        fhead = struct.pack('128sl', os.path.basename(filepath).encode('utf-8'), os.stat(filepath).st_size)

        s.send(fhead)

        fp = open(filepath, 'rb')
        while 1:
            data = fp.read(1024)
            if not data:
                print ('{0} file send over...'.format(os.path.basename(filepath)))
                break
            s.send(data)

        s.close()


def receive_file(client_socket, sym_key):
    json_data = receive_json(client_socket, sym_key)
    file_content = json_data["file_content"].encode('utf-8')
    file_name = json_data["file_name"]

    return file_name, file_content

def send_file(client_socket, sym_key, file_path):
    with open(file_path, "rb") as f:
        file_content = f.read()
    file_name = os.path.basename(file_path)
    json_data = {
        "file_name": file_name,
        "file_content": file_content.decode('utf-8')
    }
    send_json(client_socket, sym_key, json_data)


def receive_json(client_socket, sym_key):
    iv = client_socket.recv(16)
    msg_size = int.from_bytes(client_socket.recv(8), 'big')
    encrypted_data = b""
    while len(encrypted_data) < msg_size:
        packet = client_socket.recv(4096)
        if not packet:
            break
        encrypted_data += packet

    decryptor = Cipher(
        algorithms.AES(sym_key),
        modes.CFB(iv),
        backend=default_backend()
    ).decryptor()
    data = decryptor.update(encrypted_data) + decryptor.finalize()
    
    return json.loads(data.decode('utf-8'))

def send_json(client_socket, sym_key, data):
    iv = os.urandom(16)
    json_data = json.dumps(data).encode('utf-8')
    encryptor = Cipher(
        algorithms.AES(sym_key),
        modes.CFB(iv),
        backend=default_backend()
    ).encryptor()
    encrypted_data = encryptor.update(json_data) + encryptor.finalize()

    client_socket.send(iv)
    client_socket.send(len(encrypted_data).to_bytes(8, 'big'))
    client_socket.sendall(encrypted_data)


def main():



    with Manager() as manager:
        server_file_addr = manager.list()
        server_message_addr = manager.list()
        file_email = manager.dict()

        p1 = Process(target=message, args=(server_file_addr,server_message_addr,file_email))
        p2 = Process(target=socket_service, args=(server_file_addr,server_message_addr,file_email))

        p1.start()
        p2.start()

        p1.join()
        p2.join()

if __name__ == '__main__':
    main()
