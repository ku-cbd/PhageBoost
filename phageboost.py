#!/usr/bin/env python3
# Created: Fri Dec 20 16:32:00 2024
# Last changed: Time-stamp: <Last changed 2024-12-20 16:41:57 by Thomas Sicheritz-Pontén, thomas>

import os, sys, re
from icecream import ic

base = os.path.realpath(os.path.dirname(__file__))
sys.path.insert(0, base)
sys.path.insert(0, os.path.join(base, 'PhageBoost'))
from PhageBoost.main import main

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw|\.exe)?$', '', sys.argv[0])
    sys.exit(main())

