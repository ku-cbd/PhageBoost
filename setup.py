#!/usr/bin/env python
# Created:2020-06-09 17:50:51
# Last changed: Time-stamp: <Last changed 2020-06-26 16:58:52 by Kimmo Siren

import setuptools
import sys
import os

try:
    if sys.version_info.major != 3:
        sys.stderr.write(
            "Still using Python 2 ('%d'). We built this for Python 3\n" % sys.version_info.major)
        sys.exit(-1)
except Exception:
    sys.stderr.write("(hopefully you are running Python 3)\n\n")
with open("README.md", "r") as fh:
    long_description = fh.read()

# with open("requirements.txt", "r") as fh:
#    requirements = fh.read()
#requirements = requirements.replace('==', '>=').split('\n')[:-1]

# tweak to autoload the model for linux and OSX.
path = os.getcwd()
absolute_model_location = '{}/PhageBoost/models/model_delta_std_hacked.pickled.silent.gz'.format(
    path)
example_data_location = '{}/example/data/NC_000907.fasta.gz'.format(path)
with open('PhageBoost/main.py', 'r') as file:
    filedata = file.read()
    filedata = filedata.replace(
        'default_model_location', absolute_model_location)
with open('PhageBoost/main.py', 'w') as file:
    file.write(filedata)

setuptools.setup(
    name='PhageBoost',
    version='v0.1.7',
    author="Kimmo Sirén and Thomas Sicheritz-Pontén",
    author_email='kkpsiren@gmail.com',
    description="a Fast Prophage and Phage Predictor",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ku-cbd/PhageBoost",
    keywords=['machine learning', 'bioinformatics', 'phage', 'prophage',
              'bacteria', 'ngs', 'metagenomics', 'wgs', 'microbiology'],
    packages=setuptools.find_packages(),
    package_data={'PhageBoost': [
        absolute_model_location, example_data_location]},
    install_requires="biopython,joblib,more-itertools,numpy,pandas,pickleshare,pyrodigal,scipy,tables,xgboost,cachier,tabulate".split(
        ','),
    entry_points={'console_scripts': ['PhageBoost = PhageBoost.main:main']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        'Development Status :: 3 - Alpha ',
        'Topic :: Scientific/Engineering',
    ],
    python_requires='>=3.6',
)
