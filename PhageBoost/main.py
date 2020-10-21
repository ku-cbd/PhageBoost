#!/usr/bin/env python
# Created:2020-06-09 17:50:51
# Last changed: Time-stamp: <Last changed 2020-06-26 16:42:23 by Kimmo Siren

from PhageBoost import get_genecalls
from PhageBoost import get_predictions
from PhageBoost import get_matches
from PhageBoost import calc_features
import Bio.SeqIO
import gzip
import pandas as pd

import xgboost as xgb
import PhageBoost

from argparse import ArgumentParser, RawDescriptionHelpFormatter
import time
import sys
import os

#from .bbcache import bbcachier


def call_genes(fasta_seqs, meta, min_size_of_contig):
    genecalls = pd.DataFrame()
    for contig in fasta_seqs:
        if len(contig.seq) >= min_size_of_contig:
            contig_genecalls = get_genecalls.get_genecalls_for_contig(
                contig, meta)
            genecalls = pd.concat(
                [genecalls, contig_genecalls], sort=False, axis=0)
    if genecalls.shape[0] < 1:
        print('after min_size_of_contig {}, no contigs were left. This caused the run to fail'.format(
            min_size_of_contig))
        sys.exit(1)
    else:
        genecalls = genecalls.reset_index(drop=True)
        genecalls['header'] = genecalls.contig + \
            '_' + genecalls.index.astype('str')
    return genecalls


def read_sequence_file(fasta_file):
    handle = gzip.open(fasta_file, 'rt') if fasta_file.endswith(
        '.gz') else open(fasta_file)
    fasta = [x.upper() for x in Bio.SeqIO.parse(handle, 'fasta')]
    names = {}
    for entry in fasta:
        names.setdefault(entry.name, []).append(entry.name)
        count = len(names[entry.name])
        if count > 1:
            entry.name = entry.name + '_%d' % count

    meta = len(fasta) > 2
    return fasta, meta


def calculate_features(genecalls):
    df, _, _ = calc_features.df2AAandDNAfeatures(genecalls, name='header')
    return df


def read_model_from_file(model_file):
    model = get_predictions.load_model(model_file, n_jobs=1)
    feats = model.feature_names
    feats_ = [i.replace('-delta', '') for i in feats]
    limit = model.best_ntree_limit
    return model, feats, feats_, limit


def predict(model, genecalls, df, feats, period, win_type,  min_periods, limit, threshold, length, gaps, neighbouring, alpha):
    dtest = xgb.DMatrix(df[feats], feature_names=feats)
    genecalls['preds'] = model.predict(dtest, ntree_limit=limit)
    genecalls['preds20'] = 0

    for contig in genecalls['contig'].unique():
        idx = genecalls[genecalls['contig'] == contig].index
        genecalls.loc[idx, 'preds20'] = genecalls.loc[idx, 'preds'].rolling(
            period, win_type=win_type, center=True, min_periods=min_periods).mean()
    genecalls = get_predictions.get_phage(
        genecalls, threshold=threshold, length=length, gaps=gaps, neighbouring=neighbouring)

    genecalls['regions'] = get_predictions.regionise(genecalls, alpha=alpha)
    nphages = genecalls['regions'].value_counts().shape[0] - 1

    res = [get_matches.get_startstop(
        genecalls, phage) for phage in genecalls['regions'].fillna('0').unique() if phage != '0']
    res = pd.DataFrame(
        res, columns=['phage', 'contig', 'start', 'stop', 'genes', 'probability'])
    return genecalls, nphages, res


def look_for_repeats(genecalls, fasta, search_region, inwards):
    phages = [get_matches.get_startstop(
        genecalls, phage) for phage in genecalls['regions'].fillna('0').unique() if phage != '0']
    phages = pd.DataFrame(
        phages, columns=['phage', 'contig', 'start', 'stop', 'genes', 'probability'])
    res = pd.concat([get_matches.string_search(
        i, fasta, search_region, inwards) for j, i in phages.iterrows()], axis=1).T
    return res


def save(output_dir, name, genecalls, res, fasta, look_for_repeat_flag, att_size):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    genecalls.to_csv('{}/genecalls_{}.gff3'.format(output_dir, name), sep='\t')
    get_predictions.save_phage_as_fasta(
        res, fasta, look_for_repeat_flag=look_for_repeat_flag, att_size=att_size, output=output_dir)

    res_to_gff(res, output_dir, name)
    
def res_to_gff(res, output_dir, name):
    res.insert(0, 'source', 'PhageBoost')
    res.insert(0, 'type', 'prophage_prediction')
    res.insert(0, 'strand', '.')
    res.insert(0, 'phase', '.')
    res.rename({'probability':'score', 'contig':'seqid', 'stop':'end'}, axis=1, inplace=True)
    res.insert(0, 'attributes', 'n_genes=' + res.genes.astype(str) + ';phage_id=' + res.phage)
    res = res['seqid source type start end score strand phase attributes'.split()]
    csv_file = '{}/phages_{}.gff'.format(output_dir,name)
    fid = open(csv_file, 'w+')
    fid.write('##gff-version 3\n#')
    res.to_csv(fid, sep='\t', index=False)
    fid.close()

def get_version():
    from PhageBoost import __version__
    return 'PhageBoost version: {version}'.format(version=__version__)


def main():
    usage = "%prog [options] file (or - for stdin)\n"
    parser = ArgumentParser(
        usage,
        formatter_class=RawDescriptionHelpFormatter,
        epilog="""
    Example of usage:
    PhageBoost -f example/data/NC_000907.fasta.gz -o results
    """)

    parser.add_argument("-f", "--file", action="store",
                        type=str, dest="files", default=[], nargs="+")
    parser.add_argument("-o", "--output", action="store",
                        type=str, dest="output", default=None)
    parser.add_argument("-m", "--model", action="store", type=str, dest="model",
                        default=PhageBoost.modelpath())
    parser.add_argument("-j", "--threads", action="store",
                        type=int, dest='n_jobs', default=1)
    parser.add_argument("-cs", "--mincontigsize", action="store",
                        type=int, dest='min_size_of_contig', default=10000)

    parser.add_argument("-t", "--threshold", action="store",
                        type=float, dest='threshold', default=0.9)
    parser.add_argument("-l", "--length", action="store",
                        type=int, dest='length', default=10)
    parser.add_argument("-g", "--gaps", action="store",
                        type=int, dest='gaps', default=5)
    parser.add_argument("-n ", "--neighbouring", action="store",
                        type=int, dest='neighbouring', default=0)

    parser.add_argument("-r ", "--look_for_repeats", action="store",
                        type=int, dest='look_for_repeat_flag', default=0)
    parser.add_argument("-sr", "--search-region", action="store",
                        type=int, dest='search_region', default=5500)
    parser.add_argument("-i", "--inwards", action="store",
                        type=int, dest='inwards', default=500)
    parser.add_argument("-att", "--att-size", action="store",
                        type=int, dest='attsize', default=14)

    parser.add_argument("-meta", "--meta", action="store",
                        type=int, dest='meta', default=0)

    parser.add_argument("-a ", "--alpha", action="store",
                        type=float, dest='alpha', default=0.001)
    parser.add_argument('-v', '--version', action='version',
                        version=get_version())

    args = parser.parse_args()

    fasta_file = args.files[0]
    min_size_of_contig = args.min_size_of_contig

    threshold = args.threshold
    length = args.length
    gaps = args.gaps
    neighbouring = args.neighbouring

    # repeat trigger
    look_for_repeat_flag = args.look_for_repeat_flag
    search_region = args.search_region
    inwards = args.inwards
    att_size = args.attsize

    metamode = args.meta  # forces metamode no matter what contigs

    # kruskal
    alpha = args.alpha

    n_jobs = args.n_jobs  # this can be removed in future
    model_file = args.model
    output = args.output

    # rolling params # only for kruskal/interpretation/visualization purposes therefore hardcoded
    period = 20
    win_type = 'parzen'
    min_periods = 1

    timer = time.time()
    name = os.path.basename(fasta_file).split('.')[0]

    print('processing: {}'.format(name))
    fasta, meta = read_sequence_file(fasta_file)
    if metamode == 1:
        meta = True

    genecalls = call_genes(fasta, meta, min_size_of_contig)
    print('time after genecalls: {}'.format(time.time()-timer))

    df = calculate_features(genecalls)
    print('time after feature calculations: {}'.format(time.time()-timer))

    # read model here to faster debug if it breaks
    model, feats, feats_, limit = read_model_from_file(model_file)
    df = get_predictions.get_deltas(df[feats_])

    genecalls, nphages, res = predict(model, genecalls, df, feats, period,
                                      win_type, min_periods, limit, threshold, length, gaps, neighbouring, alpha)
    print('time after predictions: {}'.format(time.time()-timer))
    print(genecalls['regions'].value_counts().drop('0').to_dict())

    if nphages == 0:
        print("no phages found ")
        return 0
    if look_for_repeat_flag == 1:
        res = look_for_repeats(genecalls, fasta, search_region, inwards)
    if output:
        save(output, name, genecalls, res, fasta,
             look_for_repeat_flag, att_size)
    else:
        print('no output specified, nothing saved')


if __name__ == '__main__':
    main()
