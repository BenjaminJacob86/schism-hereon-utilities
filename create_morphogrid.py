import sys

sys.path.insert(0,'/p/home/jusers/jacob6/juwels/git/schism-hereon-utilities/')
from schism import *


s=schism_setup()
s.dump_gr3('imorphogrid.gr3',const=1.0)
