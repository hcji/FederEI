import socket
import json
from multiprocessing import Process, Manager
import sys
import struct
import os
import threading

import pickle
import numpy as np
from tqdm import tqdm
import sys
sys.path.append("..")
from data_process import spec
from data_process.spec_to_wordvector import spec_to_wordvector
import os
import gensim
from scipy.sparse import load_npz
import hnswlib
from scipy.sparse import csr_matrix, save_npz


master_ip=sys.argv[1]
master_message_port=int(sys.argv[2])
master_file_port=int(sys.argv[3])
server_ip=sys.argv[4]
message_port=int(sys.argv[5])
file_port=int(sys.argv[6])
database_path=sys.argv[7]





def main():

    spectovec,p,sorted_cleaned_pkl_data=fastei_init(database_path)
    print('fastei init')
    
    client = socket.socket()
    
    addr = (master_ip, master_message_port)


    client.connect(addr)

    data = {}
    data['identity']='server'
    data['file_addr']=server_ip
    data['file_port']=file_port
    data['message_addr']=server_ip
    data['message_port']=message_port

    json_data = json.dumps(data)
    
    client.send(json_data.encode())
    client.close()

    

    with Manager() as manager:
        file_email = manager.dict()

        p1 = Process(target=message, args=(file_email,))
        p2 = Process(target=socket_service, args=(file_email,spectovec,p,sorted_cleaned_pkl_data,))

        
        p1.start()
        p2.start()

        
        p1.join()
        p2.join()


def message(file_email):
    server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server.bind((server_ip, message_port))
    server.listen(10)

    while True:
        print("server waiting")
        
        client, addr = server.accept()
        data = client.recv(1024).decode()  
        received_data = json.loads(data)
        file_email.update(received_data)
        print(f"accept {received_data}")
            



def socket_service(file_email,spectovec,p,sorted_cleaned_pkl_data):
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
        
        t = threading.Thread(target=deal_data, args=(conn, addr,file_email,spectovec,p,sorted_cleaned_pkl_data))
        t.start()

def deal_data(conn, addr,file_email,spectovec,p,sorted_cleaned_pkl_data):
    print ('Accept new connection from {0}'.format(addr))
    
    while 1:
        
        fileinfo_size = struct.calcsize('128sl')
        
        buf = conn.recv(fileinfo_size)
        
        if buf:
            
            filename, filesize = struct.unpack('128sl', buf)
            fn = filename.strip(b'\00')
            fn = fn.decode()
            print ('file new name is {0}, filesize if {1}'.format(str(fn),filesize))
 
            recvd_size = 0  
            fp = open('./' + str(fn), 'wb')
            print ('start receiving...')
            
            
            while not recvd_size == filesize:
                if filesize - recvd_size > 1024:
                    data = conn.recv(1024)
                    recvd_size += len(data)
                else:
                    data = conn.recv(filesize - recvd_size)
                    recvd_size = filesize
                fp.write(data)
            fp.close()
            print ('end receive...')
        
        conn.close()
        break

    fastei(str(fn),spectovec,p,sorted_cleaned_pkl_data)

    socket_client(str(fn),master_ip, master_file_port)
    try :
        os.remove(str(fn))
    except:
        print("remove cache fail")
    

def fastei_init(database_path):
    model_file =os.path.join("Word2vec.model")
            
    model = gensim.models.Word2Vec.load(model_file)
    spectovec = spec_to_wordvector(model=model, intensity_weighting_power=0.5)

    if os.path.exists(database_path):
        if os.path.exists('sorted_cleaned_pkl_data.pickle') and os.path.exists('sorted_precursor_mz_list.pickle') and os.path.exists('all_predicted_word_embeddings.npz'):
            with open('sorted_cleaned_pkl_data.pickle', 'rb') as file:
                sorted_cleaned_pkl_data = pickle.load(file)
            with open('sorted_precursor_mz_list.pickle', 'rb') as file:
                sorted_precursor_mz_list = pickle.load(file)
        else:
            if database_path.lower().endswith('.mgf') :
                pkl_data = list(spec.load_from_mgf(database_path))
            elif database_path.lower().endswith('.pickle') :
                with open(database_path, 'rb') as file:
                    pkl_data = pickle.load(file)

            sorted_pkl_data = sorted(pkl_data, key=lambda x: x.metadata['precursor_mz'])
            word2vectors=[]
            word_smiles=[]
            sorted_cleaned_pkl_data=[]
            sorted_pkl_data = [s for s in sorted_pkl_data if s is not None and s.metadata['smiles']!='' and s.metadata['precursor_mz']!=0 ]
            for i in tqdm(range(len(sorted_pkl_data))):
                spectrum_in = spec.SpectrumDocument(sorted_pkl_data[i], n_decimals=0)
                try:
                    vetors=spectovec._calculate_embedding(spectrum_in)
                    word_smiles.append(spectrum_in.metadata['smiles'])
                    sorted_cleaned_pkl_data.append(sorted_pkl_data[i])
                    word2vectors.append(vetors)
                except:
                    pass
            sorted_precursor_mz_list=[]
            for i in sorted_cleaned_pkl_data:
                sorted_precursor_mz_list.append(i.metadata['precursor_mz'])

            with open('sorted_cleaned_pkl_data.pickle', 'wb') as file:
                pickle.dump(sorted_cleaned_pkl_data,file)
            with open('sorted_precursor_mz_list.pickle', 'wb') as file:
                pickle.dump(sorted_precursor_mz_list,file)
            np.save('all_smiles.npy', word_smiles)
            word_vec=csr_matrix(np.array(word2vectors))
            save_npz('all_predicted_word_embeddings.npz', word_vec)
    else:
        print("database path error")
        
    if os.path.exists('references_index.bin'):
        dim = 500
        p = hnswlib.Index(space='l2', dim=dim) 
        p.load_index('references_index.bin', max_elements =2166721)
        
    else:
        xb= load_npz('all_predicted_word_embeddings.npz').todense().astype('float32')
        xb_len =  np.linalg.norm(xb, axis=1, keepdims=True)
        xb = xb/xb_len
        dim = 500
        num_elements =len(xb)
        ids = np.arange(num_elements)
        p = hnswlib.Index(space = 'l2', dim = dim) 
        p.init_index(max_elements = num_elements, ef_construction = 800, M = 64)
        p.add_items(xb, ids)
        p.set_ef(300) 
        p.save_index('references_index.bin')
        
    
    return spectovec,p,sorted_cleaned_pkl_data

def fastei(path,spectovec,p,sorted_cleaned_pkl_data):
    
    spectrums=list(spec.load_from_mgf(path))
    
    
    word2vectors=[]
    spectrums = [s for s in spectrums if s is not None]
    for i in tqdm(range(len(spectrums))):
        spectrum_in = spec.SpectrumDocument(spectrums[i], n_decimals=0)
        vetors=spectovec._calculate_embedding(spectrum_in)

        word2vectors.append(vetors)
    
    
    word_vec=csr_matrix(np.array(word2vectors))
    save_npz('query.npz', word_vec)
    xq=load_npz('query.npz').todense().astype('float32')
    xq_len = np.linalg.norm(xq, axis=1, keepdims=True)
    xq = xq/xq_len
    p.set_ef(300) 
    k=10
    I, D = p.knn_query(xq, k )

    
    with open(path, 'w') as mgf_file:
        for index1,lis in enumerate(I):
            for index2,i in enumerate(lis):
                mgf_data = {}
                mgf_data['origin_name']=spectrums[index1].metadata['compound_name']
                mgf_data['distance']=D[index1][index2]
                mgf_data['compound_name']=sorted_cleaned_pkl_data[I[index1][index2]].metadata['compound_name']
                mgf_data['smiles']=sorted_cleaned_pkl_data[I[index1][index2]].metadata['smiles']

                
                mgf_file.write("BEGIN IONS\n")
                for key, value in mgf_data.items():
                    mgf_file.write(f"{key}={value}\n")
                for mz, intensity in zip(sorted_cleaned_pkl_data[I[index1][index2]].mz, sorted_cleaned_pkl_data[I[index1][index2]].intensities):
                    mgf_file.write(f"{mz} {intensity}\n")
                mgf_file.write("END IONS\n")
                mgf_file.write("\n")

def socket_client(filepath,address,port):
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((address, port))
    except socket.error as msg:
        print (msg)
        sys.exit(1)
    print (s.recv(1024))


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


if __name__ == '__main__':
    main()