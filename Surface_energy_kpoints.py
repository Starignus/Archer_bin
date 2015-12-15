#!/usr/bin/python
# encoding: utf-8

# el comando se ejecuta desde la terminal como:
# ./SumColumns.py file1 file2 > salida,

# para leer la entrada
import sys

def parse_line(line):
  cols = line.split()
  return [float(x) for x in cols]

file_1 = open(sys.argv[1])
file_2 = open(sys.argv[2])
print (" Give me the surface 111,110 or 100:")
surface=str(input("Surface="))
print (surface)
print (" Give me the number of layers in your slab:")
N=int(input("#N layer="))

if surface=='111':
   a=float(13.8964)  #area in Ang^2
elif surface=='100':
   a=float(8.02309024)
elif surface=='110':
   a=float(5.673182)
   print('hi')
# check that the input file is given as an argument

for line1, line2 in zip(file_1, file_2):
  energy1 = parse_line(line1)[0]
  #energy1=float(line1)
  #energy2=float(line2) In case there is just one element in the line 
  energy2 = parse_line(line2)[0]
  energy =217.98741*((energy1 - N*((energy1 - energy2)/2))/(2*a)) 
  print '{:7.5f}J/m2'.format(energy)
