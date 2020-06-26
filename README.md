# PhageBoost

## todo - critical
1. make feats faster with both numba and cache
3. add better -h flag than current
4. save model path during pip to make the -m flag unnecessary
4. check that scales to full cores
6. for web-version: hard-code everything except length, and sort the res phages with descending probabilities
7. add file downloads
7. add summary file and/or genecalls to download? 
7. if x % of contig is predicted, assign the whole contig as phage


## todo - in future
2. add nicer interface (similar to the original phageboost)
2. use dash/plotly for visualistion
5. overall beautification
6. add feeling adventurous button, to disable the .50 mean probability and only keep the std and mean and kruskal

## Introduction 
![Predictions](fig1a.png)
Prophage predictor based on gene features against a background

## Installation
```
conda create -n PhageBoost python=3.7 -y &&
conda activate PhageBoost &&
git clone git@github.com:kkpsiren/PhageBoost.git &&
cd PhageBoost/ &&
python setup.py bdist_wheel &&
pip install --user . 
```

## CLI 
```
PhageBoost -h
PBMODEL="`ls $PWD/PhageBoost/models/model_delta_std.pickled`"
PhageBoost -f example/data/NC_000907.fasta.gz -m $PBMODEL -o . --threshold 0.9 --length 10 --neighbouring 0 --gap 05
```
## Notebooks
There are some notebooks examples
