#!/usr/bin/python
# encoding: utf-8

#funci?n para que lea los vectores y los transforme a float
def read_vector(file):
    v=file.readline().split()
    for i in range(0,len(v)):
        v[i] = float(v[i])
    return v

# producto punto
def dot_prod(a,b):
    sum = 0
    for i in range(0,len(a)):
        sum= sum+ a[i]*b[i]
    return sum
    #return a[0]*b[0]+ a[1]*b[1]+a[2]*b[2] vector tama?o 3

# funci?n norma del vector
def norm(vec):
    return math.sqrt(dot_prod(vec,vec))

#angulo entre vectores       
def angle(a,b):
    a=math.degrees(math.acos(dot_prod(a,b)/(norm(a)*norm(b))))
#    if a > 90:a=180-a
    return a
# para leer la entrada
import fileinput
import sys
import math

file = open(sys.argv[1])
alat = float(sys.argv[2]) #constante de la red
file.readline()    # ignora la primera linea
A = read_vector(file)
B = read_vector(file)
C = read_vector(file)

# A, B, C, alpha, beta, gamma
#print dot_prod(A,B)/(norm(A)*norm(B))
#print "  "
#print A,"A:",alat*norm(A),dot_prod(B,C),"alpha:",angle(B,C)
#print B,"B:",alat*norm(B),dot_prod(A,C),"beta:",angle(A,C)
#print C,"C:",alat*norm(C),dot_prod(A,B),"gamma:",angle(A,B)

# A, B, C, alpha, beta, gamma: Bohr, Angtromgs
print alat, alat*0.529177
print alat*norm(A),alat*norm(B),alat*norm(C)
print alat*norm(A)*0.529177,alat*norm(B)*0.529177,alat*norm(C)*0.529177
print angle(B,C),angle(A,C),angle(A,B)

#print A,"A:",alat*norm(A),dot_prod(B,C),"alpha:",angle(B,C)
#print B,"B:",alat*norm(B),dot_prod(A,C),"beta:",angle(A,C)  
#print C,"C:",alat*norm(C),dot_prod(A,B),"gamma:",angle(A,B)

#nota:/alpha/ angulo entre los vectores *B* and *C*
#     /beta/ angulo entre los vectores *A* y *C*
#    /gamma/, angulo entre los vectores *A* and *B*
#   /alpha/ = acos( *b.c* / b*c )
