import itertools
import re
import sys


import config


RE_INFO = re.compile(';(?=\w)')
FIELD = sys.argv[1]


for line in sys.stdin:
    if line.startswith('#'):
        continue

    fields = line.split('\t')
    info = dict(i.split('=', 1) if '=' in i else (i, i)
                for i in RE_INFO.split(fields[7]))

    if 'CSQ' in info:
        anns = [dict(zip(config.VARIANT_CSQ, csq.split('|')))
                for csq in info['CSQ'].split(',')]
        vep_annotations = [ann for ann in anns
                           if ('Feature' in ann and
                               ann['Feature'].startswith('ENST'))]
        for value in set(ann[FIELD] for ann in vep_annotations):
            print('\t'.join((fields[0], fields[1], value)))
