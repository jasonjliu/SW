#!/usr/bin/python
__author__ = "Jason Liu"
__email__ = "jason.j.liu@yale.edu"
__copyright__ = "Copyright 2021"
__license__ = "GPL"
__version__ = "1.0.0"


### Usage: python hw1.py -i <input file> -s <score file> > output.txt
### Example: python hw1.py -i input.txt -s blosum62.txt > output.txt
### Note: Smith-Waterman Algorithm

# import necessary libraries

import argparse
import pandas as pd
import numpy as np

# parse the arguments, and set default gap values
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
args = parser.parse_args()

openGap   = args.opengap
extGap   = args.extgap

# variable for -infty
infty = -float("inf")

# read in the input file and separate into two strings, s1 and s2
f=open(args.input)
lines=f.readlines()
s1 = lines[0].strip('\n')
s2 = lines[1].strip('\n')

# Read in the score (blosum62) file
blosum = pd.read_csv(args.score,delim_whitespace=True)
blosum.reset_index(inplace=True)
blosum.drop("index",inplace=True,axis=1)

# Helper function to return match or mismatch score
# Find the two letters to check for the match or mismatch score in the blosum62 file
def matchOrMismatch(s1, s2, i, j):
    return blosum.iloc[blosum.columns.tolist().index(s1[j-1]),blosum.columns.tolist().index(s2[i-1])]

# Helper functions for initial values for the three matrices
def xmatrix(i, j):
    if i == 0 or j == 0:
        return infty
    else:
        return 0

def ymatrix(i, j):
    if j == 0 or i == 0:
        return infty
    else:
        return 0

def dmatrix(i, j):
    return 0

# The main calculation for Smith Waterman for two strings, s1 and s2
def sw(s1, s2, openGap, extGap):
    #define N, M as row and column lengths
    N = len(s2) + 1
    M = len(s1) + 1
    # initialize the matrices
    X = [[xmatrix(i, j) for j in range(0, M)] for i in range(0, N)]
    Y = [[ymatrix(i, j) for j in range(0, M)] for i in range(0, N)]
    D = [[dmatrix(i, j) for j in range(0, M)] for i in range(0, N)]

    # recursion for populating X, Y, D matrix

    for i in range(1, N):
        for j in range(1, M):
            X[i][j] = max((openGap + D[i-1][j]), (extGap + X[i-1][j]))
            Y[i][j] = max((openGap + D[i][j-1]), (extGap + Y[i][j-1]))

            D[i][j] = max(0, (matchOrMismatch(s1, s2, i, j) + D[i-1][j-1]), X[i][j], Y[i][j])

    return X, Y, D

## Backtrace algorithm
def backtrace(s1, s2, X, Y, D):

    # find the alignment
    seq1 = ''
    seq2 = ''
    starting_i = int(np.where(D==np.amax(D))[0])
    starting_j = int(np.where(D==np.amax(D))[1])
    i = starting_i
    j = starting_j
    while (i>0 and j>0 and D[i][j] != 0):
        if (D[i][j] == D[i-1][j-1] + matchOrMismatch(s1, s2, i, j)):
            seq1 += s1[j-1]
            seq2 += s2[i-1]
            i = i-1
            j = j-1
        elif (D[i][j] == X[i][j]):
            seq1 += '-'
            seq2 += s2[i-1]
            i -= 1
        elif (D[i][j] == Y[i][j]):
            seq1 += s1[j-1]
            seq2 += '-'
            j -= 1
      
    seq1r = seq1[len(seq1)::-1]
    seq2r = seq2[len(seq2)::-1]


    # consider the edge cases for the sequences (start)
    if (j==0 or i==0):
      if (j>=i):
        start1string=s1[0:(j-i)]
        start2string=' '*(j-i)
      elif (j<i):
        start2string=s2[0:(i-j)]
        start1string=' '*(i-j)
    else:
      if (j>=i):
        start1string=s1[0:(j)]
        start2string=' '*(j-i) + s2[0:i]
      elif (j<i):
        start2string=s2[0:(i)]
        start1string=' '*(i-j) + s1[0:j]

    # consider the edge cases for the sequences (end)

    ending1 = len(s1) - starting_j
    ending2 = len(s2) - starting_i
    if ending1==0 and ending2==0:
      end1string = ""
      end2string = ""
    elif ending2==0:
      end1string = s1[-ending1:]
      end2string = ' '*ending1
    elif ending1==0:
      end1string = ' '*ending2
      end2string = s2[-ending2:]
    else:
      end1string = s1[-ending1:]
      end2string = s2[-ending2:]


    align1string = start1string + "(" + seq1r + ")" + end1string
    align2string = start2string + "(" + seq2r + ")" + end2string

    return align1string, align2string


## Helper function to deal with the pipes ("|") found in the desired output
def match_compare(align1string, align2string):
    pipe_string = ""
    for i in range(0,len(align1string)):
      if align1string[i] == align2string[i] and (align1string[i] != "(" and align1string[i] != ")"):
        pipe_string += "|"
      else:
        pipe_string += " "
    
    return pipe_string


## Running S-W algorithm on s1, s2
X, Y, D = sw(s1,s2,openGap,extGap)


## The maximum score
maxScore = np.amax(D)

## Running the backtrace to get the alignment
align1, align2 = backtrace(s1,s2,X,Y,D)


## print the desired output format

print("""-----------
|Sequences|
-----------""")
print("sequence1")
print(s1)
print("sequence2")
print(s2)
print("""--------------
|Score Matrix|
--------------""")
x="\t".join([char for char in s1])
print("\t\t" + x)
print("\t" + "\t".join([str(char) for char in D[0]]))
for i in range(0,len(s2)):
  print(s2[i] + "\t" + "\t".join([str(char) for char in D[i+1]]))
print("""----------------------
|Best Local Alignment|
----------------------""")
print("Alignment Score:%s" %(maxScore))
print("Alignment Results:")
print(align1)
print(match_compare(align1, align2))
print(align2)


