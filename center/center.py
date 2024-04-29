import socket
import json
from multiprocessing import Process, Manager
import sys
import struct
import os
import threading
import time
import fcntl
from data_process import spec
from matchms.exporting import save_as_mgf


server_ip=sys.argv[1]
server_port=int(sys.argv[2])
file_port=int(sys.argv[3])

lock=threading.Lock()

def message(server_file_addr,server_message_addr,file_email): 
    server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server.bind((server_ip, server_port))
    server.listen(10)

    while True:
        print("server waiting")
        
        client, addr = server.accept()
        data = client.recv(1024).decode()  
        received_data = json.loads(data)
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
        conn, addr = s.accept()

        t = threading.Thread(target=deal_data, args=(conn, addr,server_file_addr,server_message_addr,file_email))
        t.start()
 
def deal_data(conn, addr,server_file_addr,server_message_addr,file_email):
    print ('Accept new connection from {0}'.format(addr))
    
    conn.send('Hi, Welcome to the server!'.encode('utf-8'))
    while 1:
        
        fileinfo_size = struct.calcsize('128sl')
        
        buf = conn.recv(fileinfo_size)
        
        if buf:
            
            filename, filesize = struct.unpack('128sl', buf)
            fn = filename.strip(b'\00')
            fn = fn.decode()
            print ('file new name is {0}, filesize if {1}'.format(str(fn),filesize))

            recvd_size = 0  

            if not os.path.exists(str(fn)) or file_email[os.path.splitext(str(fn))[0]]==0:

                fp = open('./' + str(fn), 'wb')
                fcntl.flock(fp, fcntl.LOCK_SH)
            else:

                fp = open('./' + str(fn), 'ab')
                fcntl.flock(fp, fcntl.LOCK_SH)

            print ('start receiving...')
            
            while not recvd_size == filesize:
                if filesize - recvd_size > 1024:
                    data = conn.recv(1024)
                    recvd_size += len(data)
                else:
                    data = conn.recv(filesize - recvd_size)
                    recvd_size = filesize
                fp.write(data)

            fcntl.flock(fp, fcntl.LOCK_UN)
            fp.close()

            print ('end receive...')

        break

    lock.acquire()
    name=os.path.splitext(str(fn))[0] 
    if name in file_email :
        file_email[name]+=1
    else:
        file_email[name]=0
    print(f"get file {str(fn)},now num of {name} is {file_email[name]}")
    lock.release()

    if file_email[name]==0 :
        for i in server_file_addr:
            for address,port in i.items() :
                socket_client(str(fn),address,port)

        while 1:
            if file_email[name] == len(server_file_addr) :
                break
            else:
                time.sleep(1)
        
        result = list(spec.load_from_mgf(str(fn)))

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

        save_as_mgf(sorted_result,str(fn))

        filepath=name+'.mgf' 

        fileinfo_size = struct.calcsize('128sl')

        fhead = struct.pack('128sl', os.path.basename(filepath).encode('utf-8'), os.stat(filepath).st_size)

        conn.send(fhead)

        fp = open(filepath, 'rb')
        while 1:
            data = fp.read(1024)
            if not data:
                print ('{0} client file send over...'.format(os.path.basename(filepath)))
                break
            conn.send(data)


        try :
            os.remove(str(fn))
        except:
            print("remove cache fail")

    conn.close()
    



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
