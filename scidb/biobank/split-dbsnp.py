#!/usr/bin/env python

import sys

skip_first = True

for line in sys.stdin:
  if skip_first:
    skip_first = False
    continue
  if ',' in line:
    parts = line.split()
    sub_parts = parts[4].split(',')
    for sub_part in sub_parts:
      print('\t'.join(parts[:4] + [sub_part]))
  else:
    print(line.strip())
