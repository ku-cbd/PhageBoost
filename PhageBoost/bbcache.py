#!/usr/bin/env python
# Created: Sat Jun 13 12:54:44 2020
# Last changed: Time-stamp: <Last changed 2020-06-17 16:02:19 by Kimmo Siren>

import string, re
from cachier import cachier
import datetime
import hashlib
from Bio.SeqRecord import SeqRecord
import pandas as pd

def _hash_params(args, kwargs):
    def _hash(obj):
        if isinstance(obj, pd.core.frame.DataFrame):
            return hashlib.sha256(pd.util.hash_pandas_object(obj).values.tobytes()).hexdigest()
        elif isinstance(obj, pd.core.frame.Series):
            return hashlib.sha256(pd.util.hash_pandas_object(obj).values.tobytes()).hexdigest()
        elif isinstance(obj, SeqRecord):
            return hashlib.sha256((str(obj.seq) + obj.id + obj.name).encode()).hexdigest()
        return obj

    k_args = tuple(map(_hash, args))
    k_kwargs = tuple(sorted({k: _hash(v) for k, v in kwargs.items()}.items()))
    return k_args + k_kwargs

#cache_dir = '__bbcachier_dir__'
#bbcachier = cachier(hash_params=_hash_params, stale_after=datetime.timedelta(days=3), cache_dir=cache_dir)
