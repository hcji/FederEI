import gensim
import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem

from matchms import calculate_scores
from matchms.similarity import CosineGreedy
from matchms.importing import load_from_mgf
from spec2vec import SpectrumDocument
from spec2vec import Spec2Vec
import bisect


def main(query_path,output_path='resluts.pickle',method='Cosine',peak_tolerance=0.8,mz_tolerance=1):
    reference_pos,reference_neg,sorted_mz_list_pos,sorted_mz_list_neg=Database_clean()
    spectrums_pos,spectrums_neg=Query_clean(query_path)
    if method == 'Cosine':
        pos_anno=Cosine(spectrums_pos,reference_pos,sorted_mz_list_pos,peak_tolerance,mz_tolerance)
        neg_anno=Cosine(spectrums_neg,reference_neg,sorted_mz_list_neg,peak_tolerance,mz_tolerance)
        anno=pos_anno+neg_anno
        Save_results(output_path,anno)
        
    elif method == 'Spec2vec':
        reference_formula_pos = np.array([AllChem.CalcMolFormula(Chem.MolFromSmiles(s.get('smiles'))) for s in tqdm(reference_pos)])
        reference_formula_neg = np.array([AllChem.CalcMolFormula(Chem.MolFromSmiles(s.get('smiles'))) for s in tqdm(reference_neg)])
        
        pos_anno=Spec2Vec(spectrums_pos,reference_formula_pos,'pos',sorted_mz_list_pos,peak_tolerance,mz_tolerance)
        neg_anno=Spec2Vec(spectrums_neg,reference_formula_neg,'neg',sorted_mz_list_neg,peak_tolerance,mz_tolerance)
        anno=pos_anno+neg_anno
        Save_results(output_path,anno)
    else :
        raise Exception('the method does not exist')
    
    return anno

def Database_clean():
    with open('/data/chenjie/github/MetDNA/data/references_spectrums_positive.pickle', 'rb') as file:
        reference_pos = pickle.load(file)
    with open('/data/chenjie/github/MetDNA/data/references_spectrums_negative.pickle', 'rb') as file:
        reference_neg = pickle.load(file)

    reference_pos = np.array([s for s in reference_pos if Chem.MolFromSmiles(s.get('smiles')) is not None and s.metadata['inchikey']!= ''])
    reference_neg = np.array([s for s in reference_neg if Chem.MolFromSmiles(s.get('smiles')) is not None and s.metadata['inchikey']!= ''])

    sorted_mz_list_pos=[]
    for i in reference_pos:
        sorted_mz_list_pos.append(i.metadata['precursor_mz'])
    sorted_mz_list_neg=[]
    for i in reference_neg:
        sorted_mz_list_neg.append(i.metadata['precursor_mz'])

    # formula calculate 
    #reference_formula_pos = np.array([AllChem.CalcMolFormula(Chem.MolFromSmiles(s.get('smiles'))) for s in tqdm(reference_pos)])
    #reference_formula_neg = np.array([AllChem.CalcMolFormula(Chem.MolFromSmiles(s.get('smiles'))) for s in tqdm(reference_neg)])

    return reference_pos,reference_neg,sorted_mz_list_pos,sorted_mz_list_neg

def Query_clean(query_path):
    # need peaks mz name ionmode
    spectrums = [s for s in load_from_mgf(query_path)]
    spectrums_pos=[s for s in spectrums if s.metadata['ionmode']=='positive']
    spectrums_neg=[s for s in spectrums if s.metadata['ionmode']=='negative']
    
    return spectrums_pos,spectrums_neg

def Cosine(spectrums,lib,sorted_mz_list,peak_tolerance=0.8,mz_tolerance=1):
    anno=[]
    for s in spectrums :
        spectrum=[]
        spectrum.append(s)
        try:
            left_index = bisect.bisect(sorted_mz_list,s.metadata['precursor_mz']-mz_tolerance)
            right_index = bisect.bisect(sorted_mz_list,s.metadata['precursor_mz']+mz_tolerance)
            scores = calculate_scores(spectrum, lib[left_index:right_index], CosineGreedy())
            max_score=0.0
            one_anno={}
            for (reference, query, score) in scores:
                if score[0]>=peak_tolerance :#and abs(reference.get('precursor_mz')-query.get('precursor_mz')) < mz_tolerance:
                    if score[0]>max_score:
                        one_anno['query_name']=reference.get('compound_name')
                        one_anno['lib_name']=query.get('compound_name')
                        one_anno['score']=score[0]
                        one_anno['query_mz']=reference.get('precursor_mz')
                        one_anno['lib_mz']=query.get('precursor_mz')
                        max_score=score[0]
            if one_anno != {} :
                anno.append(one_anno)
        except:
            pass
    return anno


def Spec2vec(spectrums,lib,pol,sorted_mz_list,peak_tolerance=0.8,mz_tolerance=1):
    if pol == 'neg':
        model = gensim.models.Word2Vec.load("D:/DeepMASS2_GUI_20231025/DeepMASS2_GUI/model/Ms2Vec_allGNPSnegative.hdf5")
        spec2vec_similarity = Spec2Vec(model=model, intensity_weighting_power=1, allowed_missing_percentage=100)
    elif pol == 'pos':
        model = gensim.models.Word2Vec.load("D:/DeepMASS2_GUI_20231025/DeepMASS2_GUI/model/Ms2Vec_allGNPSpositive.hdf5")
        spec2vec_similarity = Spec2Vec(model=model, intensity_weighting_power=1, allowed_missing_percentage=100)
    else:
        raise Exception('error pol type')

    anno=[]
    for s in spectrums :
        spectrum=[]
        spectrum.append(s)
        try:
            left_index = bisect.bisect(sorted_mz_list,s.metadata['precursor_mz']-mz_tolerance)
            right_index = bisect.bisect(sorted_mz_list,s.metadata['precursor_mz']+mz_tolerance)
            scores = calculate_scores(spectrum, lib[left_index:right_index], similarity_function=spec2vec_similarity)
            max_score=0.0
            one_anno={}
            for (reference, query, score) in scores:
                if score[0]>=peak_tolerance :#and abs(reference.get('precursor_mz')-query.get('precursor_mz')) < mz_tolerance:
                    if score[0]>max_score:
                        one_anno['query_name']=reference.get('compound_name')
                        one_anno['lib_name']=query.get('compound_name')
                        one_anno['score']=score[0]
                        one_anno['query_mz']=reference.get('precursor_mz')
                        one_anno['lib_mz']=query.get('precursor_mz')
                        max_score=score[0]
            if one_anno != {} :
                anno.append(one_anno)
        except:
            pass
    return anno

def Save_results(output_path,res):
    with open(output_path, 'wb') as file:
        pickle.dump(res,file)