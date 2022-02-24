'''
File name: align.py
Author: Etha Hua
Date: Feb 8th, 2022
Purpose: Implementation of the global alignment algorithm
with a naive scoring function and a linear gap genalty.
'''
import sys

# score
# Purpose: compare two letters and return their similarity score
# Parameters: p and q are two letters to be comapred
# Returns: the similarity score in integer of the two letters
def score(p,q):
    if p == q:
        return match
    else:
        return mismatch

# fillMatrix
# Purpose: fill in the row th row, col th column element of the matrix, and record
#          the corresponding traceback for that element in a traceBack matrix
# Parameters: the indices indicating the element to be filled in row and col
# Returns: N/A
# Note: requires the alignMatrix and traceBackMatrix to be defined
def fillMatrix(row, col):
    diag = alignMatrix[row - 1][col - 1] + score(s1[row - 1], s2[col - 1])
    vertGap = alignMatrix[row - 1][col] + gap
    horiGap = alignMatrix[row][col - 1] + gap 
    rtnval = max(diag, vertGap, horiGap)

    # preferrence:  diagonal > vertical > horizontal
    if rtnval == diag:
        traceBackMatrix[row].append(0)
    elif rtnval == vertGap:
        traceBackMatrix[row].append(1)
    else: # horiGap
        traceBackMatrix[row].append(2)
    return rtnval

# fillAllMatrices
# Purpose: fill in the alignment matrix with the initialization of the top row and
#          the leftmost column to be zeros and fill in the traceback matrix according
#          to the alignment matrix
# Parameters: N/A
# Returns: N/A
# Note: requires the alignMatrix and traceBackMatrix to be defined
def fillAllMatrices():
    for i in range(m+1):
        alignMatrix.append([])
        traceBackMatrix.append([])
        if i == 0:
            for j in range(n+1):
                alignMatrix[i].append(gap * j)
                traceBackMatrix[i].append(2)
        else:
            for j in range(n+1):
                if j == 0:
                    alignMatrix[i].append(gap * i)
                    traceBackMatrix[i].append(1)
                else:
                    alignMatrix[i].append(fillMatrix(i,j))

# traceBack
# Purpose: print out the aligned sequences by tracing back from the 
#          rightmost bottommost element in an alignment matrix
# Parameters: the length of s1 and the length of s2.
# Returns: N/A
def traceBack(s1len, s2len):
    alignStk = []
    i = s1len
    j = s2len
    while i > 0 or j > 0: 
        if traceBackMatrix[i][j] == 0:
            alignStk.append((s1[i - 1], s2[j - 1]))
            i -= 1
            j -= 1
        elif traceBackMatrix[i][j] == 2:
            alignStk.append(("-", s2[j - 1]))
            j -= 1
        else: # traceBackMatrix[i][j] == 1
            alignStk.append((s1[i - 1], "-"))
            i -= 1 

    s1alned = ""
    s2alned = ""
    for _ in range(len(alignStk)):
        currPair = alignStk.pop()
        s1alned += currPair[0]
        s2alned += currPair[1]
    print(s1alned)
    print(s2alned)


# start of the program

# Default values for M, m and g are set as follows
match = 4
mismatch = -2
gap = -2

# Only if the user inputted the correct number of arguments 
# for M, m and g does the program accepts them as new 
# scores for M, m and g
if len(sys.argv) == 4:
    match = int(sys.argv[1])
    mismatch = int(sys.argv[2])
    gap = int(sys.argv[3])

# Sequence strings are initialized to empty strings
seq = ["",""]

# Reading in the two sequences from stdin, ignoring the 
# commenting legends above the sequence content
index = -1
for line in sys.stdin:
    if (line[0] == '>'):
        index += 1
    else:
        while line[-1] == "\n" or line[-1] == "\r":
            line = line[0:-1]
        seq[index] += line

s1 = seq[0]
s2 = seq[1]

m = len(s1)
n = len(s2)

# Initializing the matrices
alignMatrix = []
traceBackMatrix = []

fillAllMatrices()
traceBack(m,n)
