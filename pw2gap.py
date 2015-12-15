#! /usr/bin/python
import os, sys

# By Julen Larrucea
# www.larrucea.eu

output="out"
if len(sys.argv)>1:
 for i in sys.argv:
  if i.startswith('-'):
   option=i.split('-')[1]
   if option=="o":
     output= sys.argv[sys.argv.index('-o')+1]
   if option=="h":
    print '''
    -h display this help
    -o QE output file name 
    '''
    sys.exit()

au=0.5291772083

readline=0
nkpoints=0
bulk=[]

for line in open (output,"r"):
 if readline == 0 and " k =" in line and " bands " in line:
  readline=1
  nkpoints+=1
 elif "the Fermi energy is" in line:
  EF=float(line.split()[4])
  readline=0
 elif "occupation numbers" in line:
  readline=0
 elif "SPIN" in line:
  1
 elif readline == 1 and len(line) > 5:
  if " k =" in line :
   nkpoints+=1
  else:
   for i in range(len(line.split())):
    bulk.append(float(line.split()[i]))

#print sorted(bulk)
bulk=sorted(bulk)
for x in bulk:
 if x > EF:
  lumo=float(x)
  break

print "HOMO= ", bulk[bulk.index(lumo)-1], "  LUMO=", lumo, "eV"

print "H/L GAP = ", lumo -  bulk[bulk.index(lumo)-1]
