# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 07:45:36 2022

@author: jihon
"""
import csv
import random
import socket
import struct
import datetime
import hnswlib
import numpy as np
import pandas as pd
from platform import system

import rdkit
import rdkit.Chem
import rdkit.Chem.Draw
import PyQt5

import PyQt5.Qt
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.backends.backend_qt5agg
import matplotlib.backends.backend_qt5
import matplotlib.figure

import hnswlib
import os
import sys
sys.path.append("../..")
from sklearn.preprocessing import minmax_scale, maxabs_scale
from pyteomics.mgf import MGF
from typing import Optional


from FastEI_ import Ui_Form

def load_from_mgf(filename: str) :
    """Load spectrum(s) from mgf file.
    """

    for pyteomics_spectrum in MGF(filename, convert_arrays=1):

        metadata = pyteomics_spectrum.get("params", None)
        mz = pyteomics_spectrum["m/z array"]
        intensities = maxabs_scale(pyteomics_spectrum["intensity array"])
        
        # Sort by mz (if not sorted already)
        if not np.all(mz[:-1] <= mz[1:]):
            idx_sorted = np.argsort(mz)
            mz = mz[idx_sorted]
            intensities = intensities[idx_sorted]

        yield Spectrum(mz=mz, intensities=intensities, metadata=metadata)


class Spikes:
    """
    Adapted from: https://github.com/matchms/matchms
	Copyright (c) 2021, Netherlands eScience Center
	See licenses/MATCHMS_LICENSE
    Stores arrays of intensities and M/z values, with some checks on their internal consistency.
    """
    def __init__(self, mz=None, intensities=None):
        assert isinstance(mz, np.ndarray), "Input argument 'mz' should be a np.array."
        assert isinstance(intensities, np.ndarray), "Input argument 'intensities' should be a np.array."
        assert mz.shape == intensities.shape, "Input arguments 'mz' and 'intensities' should be the same shape."
        assert mz.dtype == "float", "Input argument 'mz' should be an array of type float."
        assert intensities.dtype == "float", "Input argument 'intensities' should be an array of type float."

        self._mz = mz
        self._intensities = intensities

        assert self._is_sorted(), "mz values are out of order."

    def __eq__(self, other):
        return \
            self.mz.shape == other.mz.shape and \
            np.allclose(self.mz, other.mz) and \
            self.intensities.shape == other.intensities.shape and \
            np.allclose(self.intensities, other.intensities)

    def __len__(self):
        return self._mz.size

    def __getitem__(self, item):
        return [self.mz, self.intensities][item]

    def _is_sorted(self):
        return np.all(self.mz[:-1] <= self.mz[1:])

    def clone(self):
        return Spikes(self.mz, self.intensities)

    @property
    def mz(self):
        """getter method for mz private variable"""
        return self._mz.copy()

    @property
    def intensities(self):
        """getter method for intensities private variable"""
        return self._intensities.copy()

    @property
    def to_np(self):
        """getter method to return stacked np array of both peak mz and
        intensities"""
        return np.vstack((self.mz, self.intensities)).T


class Spectrum:

    def __init__(self, mz: np.array, intensities: np.array, metadata: Optional[dict] = None):
        """
        
    	Adapted from: https://github.com/matchms/matchms
	    Copyright (c) 2021, Netherlands eScience Center
	    See licenses/MATCHMS_LICENSE
        
        Parameters
        ----------
        mz
            Array of m/z for the peaks
        intensities
            Array of intensities for the peaks
        metadata
            Dictionary with for example the scan number of precursor m/z.
        """
        self.peaks = Spikes(mz=mz, intensities=maxabs_scale(intensities))
        if metadata is None:
            self.metadata = {}
        else:
            self.metadata = metadata

    def __eq__(self, other):
        return \
            self.peaks == other.peaks and \
            self.__metadata_eq(other.metadata)

    def __metadata_eq(self, other_metadata):
        if self.metadata.keys() != other_metadata.keys():
            return False
        for i, value in enumerate(list(self.metadata.values())):
            if isinstance(value, np.ndarray):
                if not np.all(value == list(other_metadata.values())[i]):
                    return False
            elif value != list(other_metadata.values())[i]:
                return False
        return True

    def clone(self):
        """Return a deepcopy of the spectrum instance."""
        clone = Spectrum(mz=self.peaks.mz,
                         intensities=self.peaks.intensities,
                         metadata=self.metadata)
        clone.losses = self.losses
        return clone

    def plot(self, intensity_from=0.0, intensity_to=None, with_histogram=False):
        """To visually inspect a spectrum run ``spectrum.plot()``

        """

        def plot_histogram():
            """Plot the histogram of intensity values as horizontal bars, aligned with the spectrum axes"""

            def calc_bin_edges_intensity():
                """Calculate various properties of the histogram bins, given a range in intensity defined by
                'intensity_from' and 'intensity_to', assuming a number of bins equal to 100."""
                edges = np.linspace(intensity_from, intensity_to, n_bins + 1)
                lefts = edges[:-1]
                rights = edges[1:]
                middles = (lefts + rights) / 2
                widths = rights - lefts
                return edges, middles, widths

            bin_edges, bin_middles, bin_widths = calc_bin_edges_intensity()
            counts, _ = np.histogram(self.peaks.intensities, bins=bin_edges)
            histogram_ax.set_ylim(bottom=intensity_from, top=intensity_to)
            matplotlib.pyplot.barh(bin_middles, counts, height=bin_widths, color="#047495")
            matplotlib.pyplot.title(f"histogram (n_bins={n_bins})")
            matplotlib.pyplot.xlabel("count")

        def plot_spectrum():
            """plot mz v. intensity"""

            def make_stems():
                """calculate where the stems of the spectrum peaks are going to be"""
                x = np.zeros([2, self.peaks.mz.size], dtype="float")
                y = np.zeros(x.shape)
                x[:, :] = np.tile(self.peaks.mz, (2, 1))
                y[1, :] = self.peaks.intensities
                return x, y

            spectrum_ax.set_ylim(bottom=intensity_from, top=intensity_to)
            x, y = make_stems()
            matplotlib.pyplot.plot(x, y, color="#0f0f0f", linewidth=1.0, marker="")
            matplotlib.pyplot.title("Spectrum")
            matplotlib.pyplot.xlabel("M/z")
            matplotlib.pyplot.ylabel("intensity")

        if intensity_to is None:
            intensity_to = self.peaks.intensities.max() * 1.05

        n_bins = 100
        fig = matplotlib.pyplot.figure()

        if with_histogram:
            spectrum_ax = fig.add_axes([0.2, 0.1, 0.5, 0.8])
            plot_spectrum()
            histogram_ax = fig.add_axes([0.72, 0.1, 0.2, 0.8])
            plot_histogram()
            histogram_ax.set_yticklabels([])
        else:
            spectrum_ax = fig.add_axes([0.2, 0.1, 0.7, 0.8])
            plot_spectrum()
            histogram_ax = None

        return fig

    def get(self, key: str, default=None):

        return self._metadata.copy().get(key, default)

    def set(self, key: str, value):

        self._metadata[key] = value
        return self

    @property
    def metadata(self):
        return self._metadata.copy()

    @metadata.setter
    def metadata(self, value):
        self._metadata = value

    @property
    def peaks(self) -> Spikes:
        return self._peaks.clone()

    @peaks.setter
    def peaks(self, value: Spikes):
        self._peaks = value


class MakeFigure(matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg):
    def __init__(self,width=5, height=5, dpi=300):
        self.fig = matplotlib.figure.Figure(figsize=(width, height), dpi=dpi)
        self.fig.subplots_adjust(top=0.95,bottom=0.2,left=0.15,right=0.95)
        super(MakeFigure,self).__init__(self.fig) 
        self.axes = self.fig.add_subplot(111)
        self.axes.spines['bottom'].set_linewidth(0.5)
        self.axes.spines['left'].set_linewidth(0.5)
        self.axes.spines['right'].set_linewidth(0.5)
        self.axes.spines['top'].set_linewidth(0.5)
        self.axes.tick_params(labelsize=5)
        matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg.setSizePolicy(self, PyQt5.QtWidgets.QSizePolicy.Expanding, PyQt5.QtWidgets.QSizePolicy.Expanding)
        matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg.updateGeometry(self)
        
        
    def PlotSpectrum(self, spectrum, reference = None):
        self.axes.cla()
        mz, abunds = spectrum.peaks.mz, spectrum.peaks.intensities
        abunds /= np.max(abunds)
        self.axes.vlines(mz, ymin=0, ymax=abunds, color='r', lw = 1)
        if reference is not None:
            mz1, abunds1 = reference.peaks.mz, reference.peaks.intensities
            abunds1 /= np.max(abunds1)
            self.axes.vlines(mz1, ymin = 0, ymax = -abunds1, color='b', lw = 1)
        self.axes.axhline(y=0,color='black', lw = 1)
        self.axes.set_xlabel('m/z', fontsize = 5)
        self.axes.set_ylabel('abundance', fontsize = 5)
        self.draw()

        
    def PlotComparsion(self, vector, reference):
        self.axes.cla()
        vector, reference = np.array(vector), np.array(reference)
        baseline = min(np.min(vector), np.min(reference))
        ind = np.argsort(vector)
        self.axes.plot(vector[ind] - baseline, color = 'r', lw = 1, label = 'query')
        self.axes.plot(-(reference[ind] - baseline), color = 'b', lw = 1, label = 'reference')
        self.axes.axhline(y=0,color='black', lw = 1)
        self.axes.set_xlabel('index', fontsize = 5)
        self.axes.set_ylabel('value', fontsize = 5)
        self.draw()
    


class TableModel(PyQt5.QtCore.QAbstractTableModel):
    def __init__(self, data, showAllColumn=False):
        PyQt5.QtCore.QAbstractTableModel.__init__(self)
        self.showAllColumn = showAllColumn
        self._data = data


    def rowCount(self, parent=None):
        return self._data.shape[0]


    def columnCount(self, parent=None):
        return self._data.shape[1]


    def data(self, index, role=PyQt5.QtCore.Qt.DisplayRole):
        if index.isValid():
            if role == PyQt5.QtCore.Qt.DisplayRole:
                return str(self._data.iloc[index.row(), index.column()])
        return None


    def headerData(self,col,orientation,role):
        if orientation == PyQt5.QtCore.Qt.Horizontal and role == PyQt5.QtCore.Qt.DisplayRole:
            if type(self._data.columns[col]) == tuple:
                return self._data.columns[col][-1]
            else:
                return self._data.columns[col]
        elif orientation == PyQt5.QtCore.Qt.Vertical and role == PyQt5.QtCore.Qt.DisplayRole:
            return (self._data.axes[0][col])
        return None


class Thread_LoadDatabase(PyQt5.Qt.QThread): 
    _compounds = PyQt5.QtCore.pyqtSignal(pd.DataFrame)
    _database = PyQt5.QtCore.pyqtSignal(list)
    
    def __init__(self, db_path):
        super().__init__()
        self.db_path = db_path

    def run(self):       
        pass
        


class Thread_LoadIndex(PyQt5.Qt.QThread): 
    _index = PyQt5.QtCore.pyqtSignal(hnswlib.Index)
    
    def __init__(self, spec_path):
        super().__init__()
        self.spec_path = spec_path

    def run(self):       
        pass


class FastEI(PyQt5.QtWidgets.QWidget, Ui_Form):
    
    def __init__(self, parent=None): 
        super(FastEI, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("FederEI")
        self.setWindowIcon(PyQt5.QtGui.QIcon("FederEI.ico"))
        self.ProcessBar(0, 'ready')
        
        self.figSpe = MakeFigure(1.8, 1.2, dpi = 200)
        self.figSpe_ntb = matplotlib.backends.backend_qt5.NavigationToolbar2QT(self.figSpe, self)
        self.gridlayoutfigSpec = PyQt5.QtWidgets.QGridLayout(self.groupBoxSpe)
        self.gridlayoutfigSpec.addWidget(self.figSpe)
        self.gridlayoutfigSpec.addWidget(self.figSpe_ntb)
        
        self.Labelmol = PyQt5.QtWidgets.QLabel()
        self.gridlayoutMol = PyQt5.QtWidgets.QGridLayout(self.groupBoxMol)
        self.gridlayoutMol.addWidget(self.Labelmol)


        self.pushButtonIndex.clicked.connect(self.retrieve_targets)
        self.pushButtonQue.clicked.connect(self.InputQuery)
        

        self.listWidgetQue.itemClicked.connect(self.ViewResult)
        self.tableWidgetRes.itemClicked.connect(self.PlotResult)

        
        self.allButtons = [self.pushButtonIndex,
                           self.pushButtonQue,
                           ]
        
        self.QueryList = []
        self.Database = []
        self.SpectrumList = []
        self.VectorList = []
        self.SpectrumDB = None
        self.spectovec = None
        self.compounds = None
        self.spec_bin = None
        self.ResultIndex = None
        self.ResultScore = None
        self.Finished = False
        self.resultpath = None
        self.address='127.0.0.1'
        self.port=5001

        self.Thread_LoadDatabase = None
        self.Thread_LoadIndex = None
        
        
        try:
            with open('ip.txt', 'r') as file:
                ip = file.read()
            self.plainTextEdit_prot_inp.setText(ip)
            split=ip.split(':')
            self.address=split[0]
            self.file_port=int(split[1])
        except:
            self.WarnMsg('ip.txt not found')
            
        
        self.ProcessBar(100, 'Ready!')
        
        
        
    def Load_default(self):
        self.ProcessBar(50, 'loading database...')
        
        
        
    def WarnMsg(self, Text):
        msg = PyQt5.QtWidgets.QMessageBox()
        msg.resize(550, 200)
        msg.setIcon(PyQt5.QtWidgets.QMessageBox.Warning)
        msg.setText(Text)
        msg.setWindowTitle("Warning")
        msg.exec_()    
    
    
    def ErrorMsg(self, Text):
        msg = PyQt5.QtWidgets.QMessageBox()
        msg.resize(550, 200)
        msg.setIcon(PyQt5.QtWidgets.QMessageBox.Critical)
        msg.setText(Text)
        msg.setWindowTitle("Error")
        msg.exec_()
        
        
    def InforMsg(self, Text):
        msg = PyQt5.QtWidgets.QMessageBox()
        msg.resize(550, 200)
        msg.setIcon(PyQt5.QtWidgets.QMessageBox.Information)
        msg.setText(Text)
        msg.setWindowTitle("Information")
        msg.exec_()
        
        
    def ProcessBar(self, msg, info):
        self.progressBar.setValue(int(msg))
        self.labelStatus.setText(info)


    def SetCompound(self):
        self.Finished = False
        options = PyQt5.QtWidgets.QFileDialog.Options()
        options |= PyQt5.QtWidgets.QFileDialog.DontUseNativeDialog
        
        
    def _SetCompound(self, msg):
        self.compounds = msg
    
    
    def _SetDatabase(self, msg):
        self.Database = msg
    
    
    def _SetCompoundFinished(self):
        for s in self.allButtons:
            s.setEnabled(True)
        self.ProcessBar(100, 'Ready!')
        
            

    def retrieve_targets(self):
        gene_input = self.plainTextEdit_prot_inp.text()
        if gene_input == '':
            self.ErrorMsg('No valid targets were retrieved')
            return      
        try:
            split=gene_input.split(':')
            self.address=split[0]
            self.file_port=int(split[1])   
        except:
            self.WarnMsg("ip type error")   
            return 
        self.InforMsg(f"ip set to {gene_input}")
        with open('ip.txt', 'w') as file:
            file.write(gene_input)
        

    def SetIndex(self):
        #self.Finished = False
        options = PyQt5.QtWidgets.QFileDialog.Options()
        options |= PyQt5.QtWidgets.QFileDialog.DontUseNativeDialog
        
            
            
    def _SetIndex(self, msg):
        self.spec_bin = msg    
    
    
    def _SetIndexFinished(self):
        for s in self.allButtons:
            s.setEnabled(True)
        self.ProcessBar(100, 'Finished!')            


    def SetModel(self):
        #self.Finished = False
        options = PyQt5.QtWidgets.QFileDialog.Options()
        options |= PyQt5.QtWidgets.QFileDialog.DontUseNativeDialog
        #fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self,"Load", "","Model Files (*.model)", options=options)
        #if fileName:
        self.textBrowserMod.setText("0")

    def InputQuery(self):
        self.Finished = False
        options = PyQt5.QtWidgets.QFileDialog.Options()
        options |= PyQt5.QtWidgets.QFileDialog.DontUseNativeDialog
        filepath_list, _ = PyQt5.QtWidgets.QFileDialog.getOpenFileNames(self,"Load", "","CSV Files (*.csv);;TXT Files (*.txt);;MGF Files (*.mgf)", options=options)
        
        self.QueryList = []
        self.SpectrumList = []
        self.listWidgetQue.clear()
        
        for s in self.allButtons:
            s.setEnabled(False)

        if len(filepath_list) == 0:
            pass
        else:
            
            self.QueryList, self.SpectrumList = self.ReadMSFiles(filepath_list)
            if len(self.QueryList) != len(filepath_list):
               #self.WarnMsg('Ignore invalid files')
               pass
            for fileName in self.QueryList:
                self.listWidgetQue.addItem(fileName)

            current_time = datetime.datetime.now()
            time_str = current_time.strftime("%Y-%m-%d-%H-%M-%S")
            original_filename=str(random.randint(1, 100000))
            mgfname=f"{original_filename}_{time_str}.mgf"
            with open(mgfname, 'w') as mgf_file:
                for filepath in filepath_list:
                    if filepath.lower().endswith('.csv') :
                        with open(filepath,newline='') as csvfile:
                                reader = csv.reader(csvfile)
                                mgf_data = {
                                    'compound_name': os.path.basename(filepath)
                                }
                                mgf_file.write("BEGIN IONS\n")
                                for key, value in mgf_data.items():
                                    mgf_file.write(f"{key}={value}\n")
                                for row in reader:
                                    mgf_file.write(f"{float(row[0])} {float(row[1])}\n")
                                mgf_file.write("END IONS\n")
                                mgf_file.write("\n")
                    if filepath.lower().endswith('.txt') :
                        with open(filepath,'r') as txtfile:
                            mgf_data = {
                                'compound_name': os.path.basename(filepath)
                            }
                            mgf_file.write("BEGIN IONS\n")
                            for key, value in mgf_data.items():
                                mgf_file.write(f"{key}={value}\n")
                            for line in txtfile:
                                line=line.strip()
                                values = line.split()
                                mgf_file.write(f"{float(values[0])} {float(values[1])}\n")
                            mgf_file.write("END IONS\n")
                            mgf_file.write("\n")
                    if filepath.lower().endswith('.mgf') :
                        spectrums=list(load_from_mgf(filepath));
                        for spectrum in spectrums:
                            mgf_data = spectrum.metadata
                            mgf_file.write("BEGIN IONS\n")
                            for key, value in mgf_data.items():
                                mgf_file.write(f"{key}={value}\n")
                            for mz, intensity in zip(spectrum.peaks.mz, spectrum.peaks.intensities):#查询数据的质谱
                                mgf_file.write(f"{mz} {intensity}\n")
                            mgf_file.write("END IONS\n")
                            mgf_file.write("\n")
                    
            self.ProcessBar(30, 'sending file...')
            
            
            
            try:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.connect((self.address, self.file_port))
            except socket.error as msg:
                print (msg)
                return 
                
            
            print (s.recv(1024))

            
            filepath=mgfname

            if os.path.isfile(filepath):
                
                fhead = struct.pack('128sq', mgfname.encode('utf-8'), os.stat(filepath).st_size)
                
                s.send(fhead)

                
                fp = open(filepath, 'rb')
                while 1:
                    data = fp.read(1024)
                    if not data:
                        print ('{0} file send over...'.format(mgfname))
                        break
                    s.send(data)
                
                self.ProcessBar(50, 'sending sucess,waiting response...')

                
                while 1:
                    
                    fileinfo_size = struct.calcsize('128sq')
                    
                    buf = s.recv(fileinfo_size)
                    
                    if buf:
                        
                        filename, filesize = struct.unpack('128sq', buf)
                        fn = filename.strip(b'\00')
                        fn = fn.decode()
                        print ('file new name is {0}, filesize if {1}'.format(str(fn),filesize))
            
                        recvd_size = 0  
                        fp = open('./' + str(fn), 'wb')
                        print ('start receiving...')
                        
                        
                        while not recvd_size == filesize:
                            if filesize - recvd_size > 1024:
                                data = s.recv(1024)
                                recvd_size += len(data)
                            else:
                                data = s.recv(filesize - recvd_size)
                                recvd_size = filesize
                            fp.write(data)
                        fp.close()
                        print ('end receive...')
                    
                    s.close()
                    break
            
            self.resultpath=mgfname
            self.Finished = True

            self.ProcessBar(100, 'receive file success')

        for s in self.allButtons:
            s.setEnabled(True)

    def RunProgram(self):
        options = PyQt5.QtWidgets.QFileDialog.Options()
        options |= PyQt5.QtWidgets.QFileDialog.DontUseNativeDialog
        filepath_list, _ = PyQt5.QtWidgets.QFileDialog.getOpenFileNames(self,"Load", "","CSV Files (*.csv);;TXT Files (*.txt)", options=options)
        
        if len(filepath_list) == 0 :
            return 
        
        self.ProcessBar(30, 'sending file...')

        current_time = datetime.datetime.now()
        time_str = current_time.strftime("%Y-%m-%d-%H-%M-%S")
        original_filename=str(random.randint(1, 100000))
        mgfname=f"{original_filename}_{time_str}.mgf"

        with open(mgfname, 'w') as mgf_file:
            for csvpath in filepath_list :
                with open(csvpath,newline='') as csvfile:
                    reader = csv.reader(csvfile)
                    
                    mgf_data = {
                        'compound_name': os.path.basename(csvpath)
                    }
                    mgf_file.write("BEGIN IONS\n")
                    for key, value in mgf_data.items():
                        mgf_file.write(f"{key}={value}\n")

                    for row in reader:
                        mgf_file.write(f"{float(row[0])} {float(row[1])}\n")
                    
                    mgf_file.write("END IONS\n")
        
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.connect((self.address, self.file_port))
        except socket.error as msg:
            print (msg)
            return 
            
        
        print (s.recv(1024))

        
        filepath=mgfname

        if os.path.isfile(filepath):
            
            fhead = struct.pack('128sq', mgfname.encode('utf-8'), os.stat(filepath).st_size)
            
            s.send(fhead)

            
            fp = open(filepath, 'rb')
            while 1:
                data = fp.read(1024)
                if not data:
                    print ('{0} file send over...'.format(mgfname))
                    break
                s.send(data)
            
            self.ProcessBar(50, 'sending sucess,waiting response...')

            
            

            
            while 1:
                
                fileinfo_size = struct.calcsize('128sq')
                
                buf = s.recv(fileinfo_size)
                
                if buf:
                    
                    filename, filesize = struct.unpack('128sq', buf)
                    fn = filename.strip(b'\00')
                    fn = fn.decode()
                    print ('file new name is {0}, filesize if {1}'.format(str(fn),filesize))
        
                    recvd_size = 0  
                    fp = open('./' + str(fn), 'wb')
                    print ('start receiving...')
                    
                    
                    while not recvd_size == filesize:
                        if filesize - recvd_size > 1024:
                            data = s.recv(1024)
                            recvd_size += len(data)
                        else:
                            data = s.recv(filesize - recvd_size)
                            recvd_size = filesize
                        fp.write(data)
                    fp.close()
                    print ('end receive...')
                
                text=s.recv(1024)
                self.textBrowserMod.setText(text)

                s.close()
                break
        
        self.ProcessBar(100, 'receive file success')


    def ViewResult(self):
        if not self.Finished:
            self.ErrorMsg('Please run program first!')
            return
        
        selectItem = self.listWidgetQue.currentItem()
        if not selectItem:
            self.ErrorMsg('No item is selected!')
            return
        else:
            selectItem = selectItem.text()


        
        spectrums = list(load_from_mgf(self.resultpath))
        
        data={}
        data['Rank']=1 + np.arange(len(spectrums))
        data['Distance']=[s.metadata['distance'][0:5] for s in spectrums]
        data['CompID']=[s.metadata['compound_name'] for s in spectrums]
        data['SMILES']=[s.metadata['smiles'] for s in spectrums]
        
        
        df = pd.DataFrame(data)

        


        self.FillResultWidget(df)
        

    def PlotComparsion(self, vector, reference):
        self.axes.cla()
        vector, reference = np.array(vector), np.array(reference)
        baseline = min(np.min(vector), np.min(reference))
        ind = np.argsort(vector)
        self.axes.plot(vector[ind] - baseline, color = 'r', lw = 1, label = 'query')
        self.axes.plot(-(reference[ind] - baseline), color = 'b', lw = 1, label = 'reference')
        self.axes.axhline(y=0,color='black', lw = 1)
        self.axes.set_xlabel('index', fontsize = 5)
        self.axes.set_ylabel('value', fontsize = 5)
        self.draw()


    

    def PlotResult(self):
        selectItem = self.listWidgetQue.currentItem()
        if not selectItem:
            self.ErrorMsg('No item is selected!')
            return
        else:
            selectItem = selectItem.text()        
        wh = self.QueryList.index(selectItem)
        spectrum = self.SpectrumList[wh]
        self.figSpe.PlotSpectrum(spectrum)
        
        header = [self.tableWidgetRes.horizontalHeaderItem(i).text() for i in range(self.tableWidgetRes.columnCount())]
        try:
            i = self.tableWidgetRes.selectedIndexes()[0].row()#选中的行号
        except:
            return
        
        

        spectrums = list(load_from_mgf(self.resultpath))
        spectrums=[s for s in spectrums if s.metadata['origin_name']==os.path.basename(selectItem)]
        spectrums = sorted(spectrums, key=lambda x: x.metadata['distance'])
        spectrums=spectrums[0:min(10,len(spectrums))]
        smiles=spectrums[i].metadata['smiles']

        mol = rdkit.Chem.MolFromSmiles(smiles)
        if mol is None:
            return
        
        #PILmol = Draw.MolToQPixmap(mol)
        #self.Labelmol.setPixmap(PILmol)

        rdkit.Chem.Draw.MolToFile(mol, 'temp/{}.png'.format(spectrums[i].metadata['compound_name']), wedgeBonds=False)

        pic2show = PyQt5.QtGui.QPixmap('temp/{}.png'.format(spectrums[i].metadata['compound_name']))#.scaled(new_width, new_height, Qt.KeepAspectRatioByExpanding, Qt.SmoothTransformation)
        self.Labelmol.setPixmap(pic2show)

        
        
        
        refspec=Spectrum(mz=spectrums[i].peaks.mz,intensities=spectrums[i].peaks.intensities)

        self.figSpe.PlotSpectrum(spectrum, refspec)


    def FillResultWidget(self, data):
        self.tableWidgetRes.setRowCount(data.shape[0])
        self.tableWidgetRes.setColumnCount(data.shape[1])
        self.tableWidgetRes.setHorizontalHeaderLabels(data.columns)
        self.tableWidgetRes.setVerticalHeaderLabels(data.index.astype(str))
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if type(data.iloc[i,j]) == np.float64:
                    item = PyQt5.QtWidgets.QTableWidgetItem()
                    item.setData(PyQt5.QtCore.Qt.EditRole, PyQt5.QtCore.QVariant(float(data.iloc[i,j])))
                else:
                    item = PyQt5.QtWidgets.QTableWidgetItem(str(data.iloc[i,j]))
                self.tableWidgetRes.setItem(i, j, item)


    def ReadMSFiles(self, all_files):
        spectrums_m, valid_files = [], []
        for i in range(len(all_files)):
            try:
                f = all_files[i]                
                data = pd.read_csv(f, header = None)
                M = np.array(data.iloc[:,0])
                I = np.array(data.iloc[:,1])
                I = I / np.max(I)
                keep = np.where(I > 0.001)[0]
                M = M[keep].astype(float)
                I = I[keep].astype(float)
                spectrum = Spectrum(mz=M,intensities=I, metadata={'compound_name': str(all_files[i])})
            except:
                continue
            
            valid_files.append(f)
            spectrums_m.append(spectrum)
        return valid_files, spectrums_m

    

if __name__ == '__main__':
    import sys
    
    app = PyQt5.QtWidgets.QApplication(sys.argv)
    if system() == 'Darwin':
        import _sysconfigdata__darwin_darwin
        app.setWindowIcon(PyQt5.QtGui.QIcon("FastEI.ico"))
    ui = FastEI()
    ui.show()
    sys.exit(app.exec_())
