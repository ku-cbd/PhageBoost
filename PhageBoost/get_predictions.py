#!/usr/bin/env python
# Created: 2020-06-09 17:51:48
# Last changed: Time-stamp: <Last changed 2020-07-02 12:26:41 by Thomas Sicheritz-Pontén, thomas>

import os
import pandas as pd
import pickle
import gzip
from scipy.stats import kruskal
import warnings
warnings.filterwarnings("ignore")


def parse_p(series,neighbouring=2,gaps=2):
    """
    input series has to be boolean series against a threshold
    """
    if neighbouring > 0:
        # define how many neighbouring genes forward have to be true
        sers = pd.concat([series.diff(periods=(1+i)).fillna(-1) for i in range(neighbouring)],axis=1)
        sers = ((sers.sum(1) ==0) & (series ==1)).astype('int')
        # fill the reverse neighbors
        sers = ((pd.concat([sers.shift((-1-i),fill_value=0) for i in range(neighbouring)],axis=1).sum(1) + sers)>0).astype('int')
    else:
        sers = series
    if gaps > 0:
        # close gaps backwards
        sers = ((pd.concat([sers.shift((-1-i),fill_value=0) for i in range(gaps)],axis=1).sum(1) + sers)>0).astype('int')
    else:
        sers = sers
    return sers

def get_p(series,threshold):
    regions = []
    lengths = []
    start = []
    length = 0
    for j,i in enumerate(series.values):
        if i == 1:
            length = length + 1
            start.append(j)
        elif i ==0:
            if (len(start) > 0) and (length>=threshold):
                regions.append(start)
                lengths.append(length)
                length = 0
                start = []
            else:
                length = 0
                start = []
        else:
            print('something fishy')
    if (i==1) and (length>=threshold):
        regions.append(start)
        lengths.append(length)
    startstop = [(min(region),max(region)+1) for region in regions]
    return startstop

def get_phage_contig(genecalls,threshold,length,gaps,neighbouring):
    genecalls = genecalls.copy()
    ser = (genecalls['preds'] > threshold).astype('int')
    series = parse_p(ser,neighbouring,gaps)
    startstop = get_p(series,length)
    genecalls['regions'] = '0'
    contig = genecalls['contig'].unique()[0]
    for j,i in enumerate(startstop):
        genecalls.loc[genecalls.iloc[i[0]:i[1]].index,  'regions'] = 'phage{}_{}'.format(j,contig)
    return genecalls

def get_phage(genecalls,threshold,length,gaps,neighbouring):
    """
    threshold: prediction threshold ('knob for less/more conservative')
    length: how many genes in the whole region with gaps ('less length more smaller regions')
    gaps: how many gaps allowed when parsing subregions ('allow mismatch')
    neighbouring:  how many neighbors for parsing the region (makes the preds more conservative / less probability for false positives)    
    """
    genecalls = genecalls.copy()    
    genecalls = pd.concat([get_phage_contig(genecalls[genecalls['contig']==contig],threshold=threshold,length=length,gaps=gaps,neighbouring=neighbouring) for contig in genecalls['contig'].unique()],axis=0,sort=False)
    # rename phages
    phages = [i for i in genecalls['regions'].unique() if i !='0']
    if len(phages)>0:
        mapper = pd.Series(['phage{}'.format(i) for i in range(len(phages))],index=phages)
        mapper['0'] = '0'
        genecalls['regions'] = genecalls['regions'].map(mapper)
    return genecalls

def regionise(genecalls,alpha = 0.001):
    regions = [i for i in genecalls.regions.unique() if i != '0']
    kruskall_wallis  = []
    rejection_threshold = genecalls['preds20'].mean() + genecalls['preds20'].std()
    if rejection_threshold < 0.5:
        rejection_threshold = 0.5 
    grouper = genecalls.groupby('regions')['preds20'].mean() 
    idx = grouper[grouper > rejection_threshold].index.tolist()
    if len(idx) > 0:
        idx = [i for i in idx if i !='0']
        data2 = genecalls[genecalls['regions']=='0']['preds20'].copy()
        for phage in idx:
                data1 = genecalls[genecalls['regions']==phage]['preds20'].copy()    
                stat, p = kruskal(data1, data2)
                if p <= alpha:
                    kruskall_wallis.append(phage)
    series = [i if i in kruskall_wallis else '0' for i in genecalls['regions'].values]
    return series

def get_deltas(data):
    data = (data - data.mean()) / data.std()
    data.columns = ['{}-delta'.format(i) for i in data.columns]
    return data

def load_model(model_file,n_jobs):
    model = pickle.load(gzip.open(model_file, "rb")) if model_file.endswith('.gz') else pickle.load(open(model_file, "rb")) 
    model.set_param({'predictor':'cpu_predictor','nthread': n_jobs})  # set predictions for 1CPU
    return model

def prepare_predicted_phage_regions(res, fasta,look_for_repeat_flag, att_size):
    genomes = {}
    fasta_dict = {x.name:x.seq for x in fasta}
    
    for _, phage in res.iterrows():
        contig = phage['contig']
        seq = fasta_dict.get(contig, '')
        assert len(seq) >0, "This should not happen ... why does this sequence not exists???" 
        if look_for_repeat_flag == 1:
            start, stop = (phage['start'], phage['stop']) if (phage['matchsize'] < att_size) else (phage['updated_start'], phage['updated_stop'])
        else:
            start, stop = (phage['start'], phage['stop'])
        phagename = '_'.join([contig.split('.')[0], phage['phage']])
        if look_for_repeat_flag == 1:
            phage_metadata = '_'.join([phagename, str(start), str(stop), str(phage['genes']), str(round(phage['probability'],2)), phage['string']])
        else:
            phage_metadata = '_'.join([phagename, str(start), str(stop), str(phage['genes']), str(round(phage['probability'],2))])
        #genomes.setdefault(phage['phage'], []).append((phagename, contig, phage_metadata,  start, stop, str(seq)[start:stop+1]))
        genomes.setdefault(contig, []).append((phagename, phage, phage_metadata,  start, stop, str(seq)[start:stop+1]))
    return genomes

def save_phage_as_fasta(res, fasta, look_for_repeat_flag, att_size, output, zipfile=None):
    genomes = prepare_predicted_phage_regions(res, fasta, look_for_repeat_flag, att_size)
    files = []
    for contig in genomes:
        for phagename, phage, phage_metadata, start, stop, seq in genomes[contig]:
            filename = '{}/{}.fasta'.format(output, phagename)
            with open(filename,'wt+') as fid:
                print('>{}'.format(phage_metadata), file=fid)
                print(seq, file=fid)
                #print('{} is saved as {}.fasta.gz'.format(phage, phagename))
            files.append(filename)
    if zipfile:
        from zipfile import ZipFile
        zipObj = ZipFile(zipfile, 'w')
        [zipObj.write(file, arcname=os.path.basename(file)) for file in files]
        zipObj.close()
            
if __name__ == '__main__':
    print("not working currently, checkout the notebooks or use the main")
