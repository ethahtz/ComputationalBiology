"""
gibbsSampler-i.py

Description: An implementation of the Gibbs Sampler algorithm with 
    the improved convergence criterion that ends whenever there are 
    no updates on motif positions of s* for k consecutive iterations
    
Created by Etha Hua, Feb 24 2022
"""

from lib2to3.pytree import convert
import sys
import random

# readSeqs
# Purpose: reads sequences from an input file with FASTA format
#          and put them into the seqs[] list
def readSeqs(sequences):
    currSeq = ""
    for line in sys.stdin:
        if (line[0] == '>'):
            if currSeq != "":
                sequences.append(currSeq)
                currSeq = ""
        else:
            while line[-1] == "\n" or line[-1] == "\r":
                line = line[0:-1]
            currSeq += line
    sequences.append(currSeq)


# initPositions
# Purpose: Initializes random positions for the motif position
#          list
def initPositions(numSeqs, motifLen, posArray, sequences):
    for i in range(numSeqs):
        randindex = random.randint(0, len(sequences[i]) - motifLen)
        posArray.append(randindex)


# encodeDNA
# Purpose: encodes ATCG to 0123
def encodeDNA(letter):
    if letter == 'A':
        return 0
    elif letter == 'T':
        return 1
    elif letter == 'C':
        return 2
    else: # letter = 'G'
        return 3


# buildPSSM
# Purpose: to build a position specific substitution matrix
#          for a specific iteration, with a skipIdx to indicate the
#          s* to be excluded temporarily
# Parameters: sequences is the list of all complete sequences
#             motifPos is the list of motif positions (integers) of 
#             each sequence
#             skipIdx is the index of the s* randomly chosen for this 
#             iteration
# Returns: a PSSM matrix based on the current motif sequences 
def buildPSSM(sequences, motifPos, skipIdx):
    # A, T, C, G
    currMtx = [[], [], [], []]

    # add pseudo count
    for i in range(4):
        for _ in range(motifLength):
            currMtx[i].append(pseudoCount)

    # add real counts
    for i in range(len(sequences)):
        if i == skipIdx: 
            continue
        else:
            currPos = motifPos[i]
            for j in range(motifLength):
                currLetter = sequences[i][currPos + j]
                currMtx[encodeDNA(currLetter)][j] += 1
                
    # divide by the total count of each column to calculate the freq
    # and then divide by 0.25 to get the odds ratio
    for i in range(4):
        for j in range(motifLength):
            currMtx[i][j] /= ((len(sequences) - 1) + 4 * pseudoCount)
            currMtx[i][j] /= 0.25
            currMtx[i][j] = round(currMtx[i][j], 3)
    
    return currMtx


# bestIdxStar 
# Purpose: given the PSSM matrix and a sequence, it returns 
#          a position index that produces the best motif 
#          sequence within the complete sequence based
#          on the PSSM.
# Parameters: starSeq is the sequence to be updated with
#             its motif index
#             pssmMtx is the PSSM matrix calculated for
#             this iteration
# Returns: an integer representing the best-motif index in
#          starSeq
def bestIdxStar(starSeq, pssmMtx):
    scores = []
    for i in range(len(starSeq) - motifLength):
        currScore = 1.0
        for j in range(motifLength):
            currLetter = starSeq[i + j]
            currScore *= pssmMtx[encodeDNA(currLetter)][j]
        scores.append(currScore)

    bestScore = (0, scores[0])
    for i in range(len(scores)):
        if scores[i] > bestScore[1]:
            bestScore = (i, scores[i])

    return bestScore[0]

# Start of the main
seqs = []
motifLength = 6
pseudoCount = 1
# Convergence criterion
conv = 20



# Account for user defined motif length
if len(sys.argv) == 2:
    motifLength = int(sys.argv[1])

# Account for user defined motif length and convergence
# criterion
if len(sys.argv) == 3:
    motifLength = int(sys.argv[1])
    conv = int(sys.argv[2])


readSeqs(seqs)
numSeqs = len(seqs)
motifPos = []
initPositions(numSeqs, motifLength, motifPos, seqs)

# Pick the initial s*, initialize the iteration count and 
# consecutive no-update count
sstarIdx = random.randint(0, numSeqs - 1)
numItn = 0
endNum = 0

# Start the iteration loop, and end when there are no updates 
# on the motif position of s* for [conv] consecutive iterations
while endNum < conv:
    numItn += 1
    mx = buildPSSM(seqs, motifPos, sstarIdx)
    iStar = bestIdxStar(seqs[sstarIdx], mx)
    if (iStar == motifPos[sstarIdx]):
        endNum += 1
    else:
        endNum = 0
        motifPos[sstarIdx] = iStar
    nextIdx = random.randint(0, numSeqs - 1)
    while nextIdx == sstarIdx:
        nextIdx = random.randint(0, numSeqs - 1)
    sstarIdx = nextIdx

# Outputting the result
print("Number of iterations:", numItn)
for i in range(numSeqs):
    print(seqs[i][motifPos[i]:motifPos[i]+motifLength])
