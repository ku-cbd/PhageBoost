#!/usr/bin/env python
# Created: 2020-06-09 17:50:17
# Last changed: Time-stamp: <Last changed 2020-06-26 16:38:53 by Kimmo Siren

import pandas as pd
import numpy as np
from difflib import SequenceMatcher
import gzip

#from .bbcache import bbcachier

#@bbcachier
def find_better_regions(contig,ser,search_region,inwards):
    ser = ser.copy()
    start = int(ser['start'])
    stop = int(ser['stop'])
    backtrack_r = start-search_region
    backtrack_f = stop + search_region
    if backtrack_r >= 0:
        string1 = contig[backtrack_r:(start+inwards)]
    else:
        string1 = contig[0:(start+inwards)]
    if  backtrack_f <= len(contig):
        string2 = contig[(stop-inwards):backtrack_f]
    else:
        string2 = contig[(stop-inwards):len(contig)]

    match = SequenceMatcher(None, string1, string2,autojunk=False).find_longest_match(0, len(string1), 0, len(string2))
    truestart = (start-search_region+match.a)
    truestop  = ((stop-inwards)+match.b+match.size)
    ser['updated_start'] = truestart if truestart>0 else 0
    ser['updated_stop'] = truestop if truestop < len(contig) else len(contig)
    ser['matchsize'] = match.size
    ser['string'] = string1[match.a:(match.a+match.size)]
    return ser

#@bbcachier
def get_startstop(genecalls,phage):
    genecalls1 = genecalls[genecalls['regions']==phage]
    return phage,genecalls1['contig'].unique()[0],genecalls1['start'].min(),genecalls1['stop'].max(),genecalls1.shape[0],genecalls1['preds'].mean()

def string_search(i,fasta,search_region,inwards):
    contig = [k for k in fasta if k.name == i['contig']][0]
    res_ = find_better_regions(str(contig.seq),i,search_region,inwards)
    return res_
