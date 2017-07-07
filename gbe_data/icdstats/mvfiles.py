#!/usr/bin/python -w

from __future__ import division
from __future__ import print_function

import os, sys 

f1r = open('icdinfo.txt','r').readlines()
for line in f1r[0:]:
    line = line.rstrip()
    line = line.split()
    code = line[0]
    os.system('mv /opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/icdassoc/hybrid.backup/hybrid/c*.' + str(code) + '.*rewrite.gz* /opt/biobankengine/GlobalBioBankEngineRepo/gbe_data/icdassoc/hybrid')
