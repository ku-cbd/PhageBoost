import os
import PhageBoost

_model_path = os.path.join(os.path.dirname(PhageBoost.__file__), 'models/model_delta_std_hacked.pickled.silent.gz')
_example_data_path = os.path.join(os.path.dirname(PhageBoost.__file__), 'example/data/NC_000907.fasta.gz')
def modelpath():
    return _model_path
def testdata():
    return _example_data_path

from PhageBoost import get_genecalls, get_predictions, get_matches, calc_features, main, bbcache
