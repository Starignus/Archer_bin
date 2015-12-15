#!/usr/bin/env python
# encoding: utf-8
import sys

def parse_line(line):
  cols = line.split()
  cols.pop(0)
  return [float(x) for x in cols]
# 1 index is orignal file, and the index 2 is the relaxed file
file1 = open(sys.argv[1])
file2 = open(sys.argv[2])
print 'delta_z'

cols1_old = None
cols2_old = None
for line1, line2 in zip(file1, file2):
  if 'ATOMIC_POSITIONS' in line1:
    continue
  cols1_new = parse_line(line1)
  cols2_new = parse_line(line2)
  if cols1_old is not None:
    dist1 = cols1_old[2] - cols1_new[2]
    dist2 = cols2_old[2] - cols2_new[2]
    delta = dist2/dist1 - 1
    print '{:7.2f}%'.format(100*delta)
  cols1_old = cols1_new
  cols2_old = cols2_new
  # old is the previous layer in the file (Zn)
  # new is the next layer (Z_n-1)  (d_n=Z_n-Z_n-1)