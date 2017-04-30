# This file will contain all quantum operators needed for this project.

from numpy import identity, zeros, ones
from math import sqrt

def H():
    return [[(1.0/sqrt(2)),(1.0/sqrt(2))],
            [(1.0/sqrt(2)),(-1.0/sqrt(2))]]

def I():
    return [[1,0],
            [0,1]]

def X():
    return [[0,1],
            [1,0]]

def Z():
    return [[1,0],[0,-1]]

def CNOT():
    return [[1,0,0,0],
            [0,1,0,0],
            [0,0,0,1],
            [0,0,1,0]]

def CkNOT(bit_count):
    matrix = identity(pow(2,bit_count+1), dtype=int)
    matrix[pow(2, bit_count+1)-2][pow(2, bit_count+1)-2] = 0
    matrix[pow(2, bit_count+1) - 2][pow(2, bit_count+1) - 1] = 1
    matrix[pow(2, bit_count+1) - 1][pow(2, bit_count+1) - 1] = 0
    matrix[pow(2, bit_count+1) - 1][pow(2, bit_count+1) - 2] = 1
    print matrix
    return matrix

def W(bit_count):
    matrix = ones((pow(2,bit_count),pow(2,bit_count)), dtype=float)
    for i in xrange(pow(2,bit_count)):
        for j in xrange(pow(2,bit_count)):
            if i == j:
                matrix[i][j] = (2.0/pow(2,bit_count)) - 1.0
            else:
                matrix[i][j] = (2.0/pow(2,bit_count))
    print '>> W:'
    print matrix
    return matrix