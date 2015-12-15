#! /usr/bin/python

# By Julen Larrucea
# www.larrucea.eu

import math,sys

output="out"
if len(sys.argv)>1:
 for i in sys.argv:
  if i.startswith('-'):
   option=i.split('-')[1]
   if option=="o":
     output= sys.argv[sys.argv.index('-o')+1]
   if option=="h":
    print '''
    ### pw2cellvec.py ###

    -h display this help
    -o QE output file name 
    '''
    sys.exit()

au=0.5291772083
pi=3.1415926535897932384

capture=4
fh=open(output,"r")
celli=[]

for line in fh:
 if capture<3:
  celli.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
  capture+=1

 if capture==3:
   a=math.sqrt(celli[0][0]**2+celli[0][1]**2+celli[0][2]**2)
   b=math.sqrt(celli[1][0]**2+celli[1][1]**2+celli[1][2]**2)
   c=math.sqrt(celli[2][0]**2+celli[2][1]**2+celli[2][2]**2)
   # alpha= b^c  beta= a^c gamma=a^b
   alpha=math.acos( (celli[1][0]*celli[2][0] + celli[1][1]*celli[2][1] + celli[1][2]*celli[2][2]) / (b*c))*180/pi
   beta=math.acos( (celli[0][0]*celli[2][0] + celli[0][1]*celli[2][1] + celli[0][2]*celli[2][2]) / (a*c))*180/pi
   gamma=math.acos( (celli[0][0]*celli[1][0] + celli[0][1]*celli[1][1] + celli[0][2]*celli[1][2]) / (a*b))*180/pi
   Vol=abs((celli[0][0]*celli[1][1]*celli[2][2])+(celli[0][1]*celli[1][2]*celli[2][0])+(celli[0][2]*celli[1][0]*celli[2][1])-(celli[0][2]*celli[1][1]*celli[2][0])-(celli[0][1]*celli[1][0]*celli[2][2])-(celli[0][0]*celli[1][2]*celli[2][1]))*(alat*au)**3

   print "AA a=",a*alat*au, " b=",b*alat*au, " c=", c*alat*au, "   Vol.=", Vol
   print "au a=",a*alat, " b=",b*alat, " c=", c*alat, "   c/a=", c/a
   print " alpha=", alpha, "  beta=",beta, "  gamma=", gamma
   capture+=1

 if "CELL_PARAMETERS" in line:
   celli=[]
   alat=float(line.split()[2].split(")")[0])
   capture=0
   print " - - "

if len(celli) > 0 :   
 print "" 
 print "final lattice vectors (Angstrom)" 
 for i in celli:
  print i[0]*alat*au,i[1]*alat*au,i[2]*alat*au

