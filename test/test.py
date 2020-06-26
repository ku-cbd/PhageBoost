#!/usr/bin/env python
# Created: Sat Jun 13 14:03:36 2020
# Last changed: Time-stamp: <Last changed 2020-06-15 11:42:00 by Thomas Sicheritz-Pontén, thomas>

import string, re
import os, sys
sys.path.insert(0, './PhageBoost')
import glob, gzip
import time
from PhageBoost import main as pb
from tabulate import tabulate



period = 20
win_type = 'parzen'
min_periods = 1
threshold = 0.5
min_size_of_contig = 10000
length=10
gaps = 1
neighbouring = 1
search_region = 5500
inwards = 500
att_size = 14
alpha = 0.001

look_for_repeat_flag = False

model_file = 'PhageBoost/models/model_delta_std_hacked.pickled.gz'
fasta_file = 'example/data/NC_000907.fasta.gz'
output = 'ttt'
timer = time.time()
name = os.path.basename(fasta_file).split('.')[0]

print('processing: {}'.format(name))
fasta, meta = pb.read_sequence_file(fasta_file)

print("Gene call")
genecalls = pb.call_genes(fasta, meta, min_size_of_contig)
print('time after genecalls: {}'.format(time.time()-timer))

df = pb.calculate_features(genecalls)
print('time after feature calculations: {}'.format(time.time()-timer))

model, feats, feats_, limit = pb.read_model_from_file(model_file)
df = pb.get_predictions.get_deltas(df[feats_])
genecalls, nphages, res = pb.predict(model, genecalls, df, feats, period, win_type, min_periods, limit, threshold, length, gaps, neighbouring, alpha)
print('time after predictions: {}'.format(time.time()-timer))

if nphages == 0:
    print("no phages found ")

else:
    res = pb.look_for_repeats(genecalls, fasta, search_region, inwards)
    print('time after looking for repeats: {}'.format(time.time()-timer))
    #pb.save(output, name, genecalls, res, fasta, att_size)

    fasta_dict = {x.id:x.seq for x in fasta}
    genomes = pb.get_predictions.prepare_predicted_phage_regions(res, fasta, look_for_repeat_flag, att_size)
    contigs = genomes.keys()

    genome = []
    for contig in contigs:
        coordinates = [x[3:5] for x in genomes[contig]]
        seq = fasta_dict[contig]

        df = pd.concat([x[:-1][1] for x in genomes[contig]], axis=1).drop('contig').T
        df.set_index('phage', inplace=True)
        df.index.name = name
        print("Predicted %d prophage%s for %s" % (len(df), 's' if len(df)>1 else '', contig))
        print(tabulate(df, headers='keys', tablefmt='fancy_grid'))

        print("="*80)
