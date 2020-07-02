#!/usr/bin/env python
# Created: 2020-06-09 17:51:11
# Last changed: Time-stamp: <Last changed 2020-06-26 16:38:59 by Kimmo Siren

import pyrodigal
import pandas as pd
#from .bbcache import bbcachier

def parse_genecall(contig,gene,i):
    name = contig.name
    start = gene.begin
    stop = gene.end
    direction = gene.strand
    partial = 1 if gene.partial_end or gene.partial_begin else 0
    nt = str(contig.seq[(start-1):stop].reverse_complement()) if direction == -1 else str(contig.seq)[(start-1):stop]
    aa = str(gene.translate()).strip('*')
    return [name,i,start,stop,direction,partial,nt,aa]

#@bbcachier
def get_genecalls_for_contig(contig, meta=True):
    seq = str(contig.seq)
    p = pyrodigal.Pyrodigal(meta=meta)
    if not meta: p.train(seq)
    a = p.find_genes(seq)
    genecalls = [parse_genecall(contig, gene, i) for i, gene in enumerate(p.find_genes(seq))]
    idx =['contig','id','start','stop','direction','partial','DNAseq','AAseq']
    genecalls =  pd.DataFrame(genecalls,columns=idx) 
    return genecalls
