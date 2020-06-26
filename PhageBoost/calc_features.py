#!/usr/bin/env python
# Created: Thu Jun  4 22:34:47 2020
# Last changed: Time-stamp: <Last changed 2020-06-26 16:39:24 by Kimmo Siren>

import string, re
import os, sys
import glob, gzip
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data import IUPACData
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex, CodonsDict
from Bio.SeqUtils.MeltingTemp import Tm_NN
from Bio.SeqUtils import GC, GC_skew, GC123
from Bio.Seq import Seq
from zlib import compress
import pandas as pd
import numpy as np
from math import log
import itertools
from more_itertools import chunked
from collections import namedtuple

#from .bbcache import bbcachier


property_residues  = {'Polar': ['D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T', 'Z'], 'Aliphatic': ['A', 'I', 'L', 'V'], 'Aromatic': ['F', 'H', 'W', 'Y'], 'Basic': ['H', 'K', 'R'], 'Small': ['A', 'B', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V'], 'Acidic': ['B', 'D', 'E', 'Z'], 'Charged': ['B', 'D', 'E', 'H', 'K', 'R', 'Z'], 'Tiny': ['A', 'C', 'G', 'S', 'T'], 'Non-polar': ['A', 'C', 'F', 'G', 'I', 'L', 'M', 'P', 'V', 'W', 'Y']}
dayhoff_freq = {'A': 8.6, 'C': 2.9, 'E': 6.0, 'D': 5.5, 'G': 8.4, 'F': 3.6, 'I': 4.5, 'H': 2.0, 'K': 6.6, 'M': 1.7, 'L': 7.4, 'N': 4.3, 'Q': 3.9, 'P': 5.2, 'S': 7.0, 'R': 4.9, 'U': 0.1, 'T': 6.1, 'W': 1.3, 'V': 6.6, 'Y': 3.4}
murphy_10_tab = {'A': 'A', 'C': 'C', 'E': 'E', 'D': 'E', 'G': 'G', 'F': 'F', 'I': 'L', 'H': 'H', 'K': 'K', 'M': 'L', 'L': 'L', 'N': 'E', 'Q': 'E', 'P': 'P', 'S': 'S', 'R': 'K', 'T': 'S', 'W': 'F', 'V': 'L', 'Y': 'F'}

def biopython_proteinanalysis_seq(seq, scaling=False):
    res = ProteinAnalysis(seq)
    d = {}
    flex = np.array(res.flexibility())
    d['flex:min'], d['flex:max'], d['flex:std'] = flex.min(), flex.max(), flex.std()
    d['gravy'] = res.gravy()
    d['instability_index'] = res.instability_index()
    d['isoelectric_point'] = res.isoelectric_point()
    r, c = res.molar_extinction_coefficient()
    d['molar_extinction_coefficient_reduced'], d['molar_extinction_coefficient_cysteines'] = r, c
    d['molecular_weight'] = res.molecular_weight()
    d['percent_helix_naive'],d['percent_turn_naive'],d['percent_strand_naive']   = res.secondary_structure_fraction()

    aap = res.get_amino_acids_percent()
    aas = sorted(aap.keys())
    d.update({'percent:%s' % aa:aap[aa] for aa in aas})
    d.update({'prop_res_%s' % key:sum([aap.get(x, 0) for x in value]) for key, value in list(property_residues.items())})
    return d

class myCodonAdaptationIndex(CodonAdaptationIndex):
    def __init__(self, verbose = False):
        CodonAdaptationIndex.__init__(self)
        self.verbose = verbose
        
    def _count_codons(self, sequences): 
        self.codon_count = CodonsDict.copy() 
        for seq in sequences:
            if str(seq).islower(): 
                dna_sequence = str(seq).upper() 
            else: 
                dna_sequence = str(seq) 
                for i in range(0, len(dna_sequence), 3): 
                  codon = dna_sequence[i:i + 3]
                  if len(codon) != 3: continue
                  if codon in self.codon_count: 
                    self.codon_count[codon] += 1 
                  else:
                    if self.verbose:
                      print("illegal codon %s in gene" % (codon), dna_sequence, Seq(dna_sequence).translate(), file=sys.stderr)

def calculate_CAI(entries):
    name, sequences = list(zip(*entries))
    CAI = myCodonAdaptationIndex()
    CAI.generate_index(sequences)    
    cais = []
    for seq in sequences:
        try:
            cai = CAI.cai_for_gene(str(seq))
        except TypeError:
            cai = 1
        cais.append(cai)
    _cais = [x for x in cais if x < 1]
    if not len(_cais): return
    medium_score = sum(_cais)/len(_cais)
    return pd.Series([medium_score if x == 1 else x for x in cais], name='CAI')

def entropy(seq, is_protein = False):
    alphabet = list(dayhoff_freq.keys()) if is_protein else "ACGT" 
    cnt = [seq.count(i) for i in alphabet]
    d = sum(cnt)
    ent = []
    for i in [float(i)/d for i in cnt]:
        if i == 0:
            i = 1
        ent.append(i * log(i, 2))
    return -1 * sum(ent)

def calcultate_compressebility(entries):
    name, sequences = list(zip(*entries))        
    return pd.Series([len(compress(x[:50].encode('ASCII')))/50.0 for x in sequences], name='compressebility_minsize_50')

def calculate_entropy(entries, is_protein=False):
    name, sequences = list(zip(*entries))        
    return pd.Series([entropy(x, is_protein=is_protein)/2.0 for x in sequences], name='entropyAA' if is_protein else 'entropy')

def calculcate_melting_temperature(entries):
    name, sequences = list(zip(*entries))    
    return pd.Series([Tm_NN(x, check=True)/100.0 for x in sequences], name='melting_temp')

def calculate_GC(entries):
    gc, gc1, gc2, gc3 = [], [], [], []
    name, sequences = list(zip(*entries))
    for seq in sequences:
        _gc, _gc1, _gc2, _gc3 = GC123(seq)
        gc.append(_gc)
        gc1.append(_gc1)
        gc2.append(_gc2)
        gc3.append(_gc3)
    df = pd.DataFrame([[x/100.0 for x in gc], [x/100.0 for x in gc1], [x/100.0 for x in gc2], [x/100.0 for x in gc3]]).T
    df.columns = "gc_content gc1_content gc2_content gc3_content".split()
    return df

def calculate_GCskew(entries):
    name, sequences = list(zip(*entries))
    total_seq = ''.join(sequences)
    L = len(total_seq)
    window = L//len(sequences)
    gc_skew_values = GC_skew(total_seq, window=window)
    return pd.Series(gc_skew_values[:len(sequences)], name='gc_skew')

def calculcate_gene_lengths(entries, locations):
    name, sequences = list(zip(*entries))
    s = [locations[x+1][0]-locations[x][1] for x in range(len(sequences)-1)]
    s.append(np.mean(s))
    s = pd.Series(s)
    diffs = np.sign(s)*np.log(abs(s)+0.001)/10
    df = pd.DataFrame([[log(len(x))/14.0 for x in sequences], [x[2] for x in locations], list(diffs)]).T
    df.columns = "gene_lengths gene_strands gene_diffs".split()
    return  df

def oligopetide_frequencies(seq, n=2, alphabet="ACEDGFIHKMLNQPSRTWVY"):
    L = float(len(seq))
    nmers = sorted([''.join(x) for x in list(itertools.product(list(alphabet), repeat=n))])
    return nmers, [seq.count(x)/L for x in nmers]

def dipeptide_frequencies(entries, prefix="DIPEP", alphabet="ACEDGFIHKMLNQPSRTWVY"):
    m = []
    for n, (name, seq) in enumerate(entries):
        nmers, freqs = oligopetide_frequencies(seq, n=2, alphabet=alphabet)
        values = freqs
        m.append(values)
    df = pd.DataFrame(m)
    df.columns = ['%s:%s' % (prefix, x) for x in nmers]
    return df

def tripeptide_frequencies(entries, prefix = "TRIPEP", alphabet="ACEDGFIHKMLNQPSRTWVY"):
    m = []
    for n, (name, seq) in enumerate(entries):
        nmers, freqs = oligopetide_frequencies(seq, n=3, alphabet=alphabet)
        values = freqs
        m.append(values)
    df = pd.DataFrame(m)
    df.columns = ['%s:%s' % (prefix, x) for x in nmers]        
    return df

def biopython_proteinanalysis(entries, scaling=False):
    m = []
    for name, seq in entries:
        d = biopython_proteinanalysis_seq(seq, scaling=scaling)
        m.append([d[x] for x in list(d.keys())])
    df = pd.DataFrame(m)
    df.columns = ["BP:" + x for x in list(d.keys())]
    return df
    
def RunAA(entries, scaling=False, verbose = False, splitted=False):
    reduced_entries = [(name,''.join([murphy_10_tab.get(x,x) for x in seq])) for name,seq in entries]
    names, sequences = list(zip(*entries))
    df_entropy = calculate_entropy(entries, is_protein=True)
    df_biopython = biopython_proteinanalysis(entries, scaling=scaling)
    df_dipeptide_frequencies = dipeptide_frequencies(entries, prefix="DIPEP")
    reduced_alphabet = ''.join(set(murphy_10_tab.values()))
    df_reduced_dipeptide_frequencies = dipeptide_frequencies(reduced_entries, alphabet = reduced_alphabet, prefix="RED_DIPEP")
    df_reduced_tripeptide_frequencies = tripeptide_frequencies(reduced_entries, alphabet = reduced_alphabet, prefix="RED_TRIPEP")
    DF = pd.concat([df_entropy,df_biopython,df_dipeptide_frequencies,df_reduced_dipeptide_frequencies,df_reduced_tripeptide_frequencies,], axis = 1)
    DF.index = names
    return DF

def RunDNA(entries, locations=None, verbose = False):
    names, sequences = list(zip(*entries))
    df_CAI = calculate_CAI(entries)
    df_GC = calculate_GC(entries)
    df_GC_skew = calculate_GCskew(entries)    
    df_melting_temp = calculcate_melting_temperature(entries)
    df_gene_length = calculcate_gene_lengths(entries, locations)
    df_compressebility = calcultate_compressebility(entries)
    df_entropy = calculate_entropy(entries)
    DF = pd.concat([df_CAI,df_melting_temp,df_gene_length,df_GC,df_GC_skew,df_compressebility,df_entropy], axis = 1)
    DF.index = names
    return DF
    
def RunAAandDNA(DNA_entries, AA_entries, locations, verbose = False):
    DF_DNA = RunDNA(DNA_entries, locations, verbose = verbose)
    DF_AA = RunAA(AA_entries, verbose = verbose)
    DF = pd.concat([DF_AA, DF_DNA], axis=1)
    return DF, DF_AA, DF_DNA

def df2entries(df, name = 'contig'):
    aa_entries = list(df[[name, 'AAseq']].values)
    dna_entries = list(df[[name, 'DNAseq']].values)
    locations = list(df[['start', 'stop', 'direction']].values)
    return (dna_entries, aa_entries, locations)

#@bbcachier
def df2AAandDNAfeatures(df, name='contig'):
    dna_entries, aa_entries, locations = df2entries(df, name=name)
    DF, DF_AA, DF_DNA = RunAAandDNA(dna_entries, aa_entries, locations)
    return DF, DF_AA, DF_DNA
