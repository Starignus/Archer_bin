#!/usr/bin/python
# encoding: utf-8
#
####################################################################
#
# Script to multiply the primitive lattice vectors from QE for alat
# and tranform to angstroms.
# Usage: ./primvect_multip.py prim_vect.dat 
#
####################################################################
#Author: Dr. Ariadna Blanca Romero
#        Postdoctoral Research Associate
#        Imperial College London
#        Thomas Young Centre-Chemestry
#        ariadna@starignus.com or starignus@gmail.com
#        https://github.com/Starignus
####################################################################

def read_alat_line(file):
    # assumes that the next line in the file something like:
    #   'CELL_PARAMETERS (alat= 12.51221877)'
    # and splits the line in 3 words:
    #   ['CELL_PARAMETERS', '(alat=', '12.51221877)']
    words = file.readline().split()
    # from the third word (index 2), remove any trailing ')',
    # convert into a float and return that value.
    return float(words[2].rstrip(')'))

#función para que lea los vectores y los transforme a float
def read_vector(file):
    v=file.readline().split()
    for i in range(0,len(v)):
        v[i] = float(v[i])
    return v


# para leer la entrada
import fileinput
import sys
import math

file = open(sys.argv[1])
alat = read_alat_line(file)        # read lattice parameter a
A = read_vector(file)
B = read_vector(file)
C = read_vector(file)

print "lattice alat (a.u./Ang):"
print alat,",", alat*0.529177
print "In Bohor:"
print" "
print A[0]*alat, A[1]*alat, A[2]*alat
print B[0]*alat, B[1]*alat, B[2]*alat
print C[0]*alat, C[1]*alat, C[2]*alat
print" "
print "In Angstroms:"
print" "
print A[0]*alat*0.529177, A[1]*alat*0.529177, A[2]*alat*0.529177
print B[0]*alat*0.529177, B[1]*alat*0.529177, B[2]*alat*0.529177
print C[0]*alat*0.529177, C[1]*alat*0.529177, C[2]*alat*0.529177
