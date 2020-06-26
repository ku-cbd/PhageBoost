# PhageBoost

## Introduction 
![Predictions](fig1a.png)
Prophage predictor based on gene features against a background. 

Prophages are phages integrated into prokaryotic genomes that drive many aspects of bacterial biology.  Their extreme diversity means they are challenging to detect using sequence similarity. We present a novel fast and generalizing machine learning method to facilitate novel phage discovery.

## Publications

to be added  


## Getting Started
### Installation
```
conda create -n PhageBoost python=3.7 -y &&
conda activate PhageBoost &&
git clone git@github.com:kkpsiren/PhageBoost.git &&
cd PhageBoost/ &&
python setup.py bdist_wheel &&
pip install --user . 
```

### CLI 
```
PhageBoost -h
PhageBoost -f example/data/NC_000907.fasta.gz -o results
```
### Notebooks
There are basic notebook examples in the ```notebooks/```
These notebooks provide a way how to bring your own genecalls to PhageBoost.
You can connect your PhageBoost kernel to Jupyter via ipykernel:

```
conda activate PhageBoost
pip install ipykernel
python -m ipykernel install --user --name PhageBoost --display-name "PhageBoost" 
```



## todo - critical
1. make feats faster with both numba and cache
3. add better -h flag than current
6. for web-version: hard-code everything except length, and sort the res phages with descending probabilities
7. if x % of contig is predicted, assign the whole contig as phage
8. parse regions together if close.
9. proper documentation

## todo - in future
2. add a nicer CLI (similar to the alpha phageboost)
2. use dash/plotly for visualistion
5. overall beautification
6. add feeling adventurous button, to disable the .50 mean probability and only keep the std and mean and kruskal
