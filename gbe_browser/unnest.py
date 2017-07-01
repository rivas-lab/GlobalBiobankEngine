import fileinput
import itertools
import re


RE_INFO = re.compile(';(?=\w)')  # TODO is \w needed?
EMPTY_CSQ = ('0',) + ('',) * 67


for line in fileinput.input():
    if line.startswith('#'):
        continue

    fields = line.split('\t')
    info = dict(i.split('=', 1) if '=' in i else (i, i)
                for i in RE_INFO.split(fields[7]))

    if 'CSQ' in info:           # TODO is CSQ missing sometimes?
        for idx, csq in enumerate(info['CSQ'].split(',')):
            print('\t'.join(itertools.chain(
                fields[:7],
                (str(idx),),
                csq.split('|'))))

    else:
        print('\t'.join(itertools.chain(
            fields[:7],
            EMPTY_CSQ)))
