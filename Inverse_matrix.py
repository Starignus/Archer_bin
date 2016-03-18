#!/usr/bin/python
# encoding: utf-8
#
####################################################################
#
# Inverse of a  3x3 matrices
# Usage: ./matrix_multipl_2.py matrixA.dat matrixB.dat
#
####################################################################
#Author: Dr. Ariadna Blanca Romero
#        Postdoctoral Research Associate
#        Imperial College London
#        Thomas Young Centre-Chemestry
#        ariadna@starignus.com or starignus@gmail.com
#        https://github.com/Starignus
####################################################################


import numpy

#function that changes a line into a vector of floats
def parse_vector(line):
    return [float(v) for v in line.split()]

def read_matrix(filename):
    file = open(filename)
    line = file.readline()
    if 'alat' in line:
        alat = float(line.split()[2].rstrip(')'))
        line = file.readline()
    else:
        alat = 1.0
    A = parse_vector(line)
    B = parse_vector(file.readline())
    C = parse_vector(file.readline())
    M = numpy.array([A, B, C])
    return alat * M

def print_matrix(M):
    for row in M:
        print '\t'.join([str(v) for v in row])

# For reading inputs
import fileinput
import sys
import math
R = read_matrix(sys.argv[1]) # Reading matrix 

# Inverse 3x3
print 'Inverse'
print_matrix(numpy.linalg.inv(R))

# Transpose
print ''
print 'Transpose'
print_matrix(R.transpose())

# Determinant 
print ''
print 'Determinant'
print numpy.linalg.det(R)
