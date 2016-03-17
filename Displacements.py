#!/usr/bin/env python
# encoding: utf-8


######################################################
#
# Calcuating interlayer spacig in the slab
# after relaxation.
#
#######################################################
#Author: Dr. Ariadna Blanca Romero
#        Postdoctoral Research Associate
#        Imperial College London
#        Thomas Young Centre-Chemestry
#        ariadna@starignus.com or starignus@gmail.com
#        https://github.com/Starignus
#######################################################


import sys

def parse_line(line):
  cols = line.split()
  cols.pop(0)
  return [float(x) for x in cols]

file1 = open(sys.argv[1])
file2 = open(sys.argv[2])
print '       Dx       Dy         Dz'

for line1, line2 in zip(file1, file2):
  if 'ATOMIC_POSITIONS' in line1:
    continue
  cols1 = parse_line(line1)
  cols2 = parse_line(line2)
  cols = [x1 - x2 for x1, x2 in zip(cols1, cols2)]
  print ' '.join(str(x) for x in cols)
