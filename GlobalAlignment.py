# -*- coding: utf-8 -*-
"""
Filename: GlobalAlignment.py
Date created: 2021/05/29, Sat, 12:00:00 (UTC+8)
@author: Rishika Gupta
Purpose: Step-by-step explanation of dynamic programming solution to sequence alignment.
Steps: 
https://python.plainenglish.io/global-sequence-alignment-implementation-in-python-from-scratch-5a8a611dbb1e
https://github.com/risg99/Local-And-Global-Sequence-Alignment-Python-Scratch/blob/main/GlobalAlignment.py
"""
import numpy as np
import pandas as pd
# To adjust the dataframe appearance
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option("display.width", 200)
pd.set_option('display.expand_frame_repr', False)

class ScoreParams:
    # Define scores for each parameter
    def __init__(self,gap,match,mismatch):
        self.gap = gap
        self.match = match
        self.mismatch = mismatch


    def misMatchChar(self,x,y):
        if x != y:
            return self.mismatch
        else:
            return self.match


def getMatrix(x,y,gap):
    # Create an initial matrix of zeros, with dim of len(x) by len(y).
    matrix = []
    sizeX = len(x)
    sizeY = len(y)
    for i in range(sizeX+1):
        subMatrix = []
        for j in range(sizeY+1):
            subMatrix.append(0)
        matrix.append(subMatrix)

    # Initializing the first row and first column with the gap values
    for j in range(1,sizeY+1):
        matrix[0][j] = j*gap
    for i in range(1,sizeX+1):
        matrix[i][0] = i*gap
    return matrix


def getTraceBackMatrix(x,y):
    # Create an initial matrix of '0's, with dim of len(x) by len(y).
    matrix = []
    sizeX = len(x)
    sizeY = len(y)
    for i in range(sizeX+1):
        subMatrix = []
        for j in range(sizeY+1):
            subMatrix.append('0')
        matrix.append(subMatrix)

    # Initializing the first row and first column with the up or left values
    for j in range(1,sizeY+1):
        matrix[0][j] = '←'
    for i in range(1,sizeX+1):
        matrix[i][0] = '↑'
    matrix[0][0] = '✓'
    return matrix


def globalAlign(x,y,score):
    '''
    Fill in the matrix with alignment scores
    '''
    matrix = getMatrix(x,y,score.gap)
    traceBack = getTraceBackMatrix(x,y)

    for i in range(1,len(x)+1):
        for j in range(1,len(y)+1):
            left = matrix[i][j-1] + score.gap
            up = matrix[i-1][j] + score.gap
            diag = matrix[i-1][j-1] + score.misMatchChar(x[i-1],y[j-1])
            matrix[i][j] = max(left,up,diag)
            if matrix[i][j] == left:
                traceBack[i][j] = '←'
            elif matrix[i][j] == up:
                traceBack[i][j] = '↑'
            else:
                traceBack[i][j] = '↖'
    return matrix, traceBack

# Use '-' and ' ' to handle numbers.
def getAlignedSequences(x,y,matrix,traceBack):
    # Obtain x and y globally aligned sequence arrays using the bottom-up approach.
    xSeq = []
    ySeq = []
    i = len(x)
    j = len(y)
    while(i > 0 or j > 0):
        if traceBack[i][j] == '↖':
            # Diag is scored when x[i-1] == y[j-1]
            xSeq.append(x[i-1])
            ySeq.append(y[j-1])
            i = i-1
            j = j-1
        elif traceBack[i][j] == '←':
            # Left holds true when '-' is added from x string and y[j-1] from y string
            xSeq.append('-' * len(y[j-1]))
            ySeq.append(y[j-1])
            j = j-1
        elif traceBack[i][j] == '↑':
            # Up holds true when '-' is added from y string and x[j-1] from x string
            xSeq.append(x[i-1])
            ySeq.append('-' * len(x[i-1]))
            i = i-1
        elif traceBack[i][j] == '✓':
            # Break condition when we reach the [0,0] cell of traceback matrix
            break
    return xSeq,ySeq


def printMatrix(matrix):
    # Create a custom function to print the matrix
    for i in range(len(matrix)):
        print(matrix[i])
    print()


def driver(x, y, gap=-1, match=1, mismatch=-1, debug=False):
    # Driver Code.
    score = ScoreParams(gap,match,mismatch)
    matrix, traceBack = globalAlign(x,y,score)
    xSeq,ySeq = getAlignedSequences(x,y,matrix,traceBack)
    # Change alignment from horizontal to vertical (pandas).
    if debug:
        print('Input sequences are: ')
        print(x)
        print(y)

        print('Printing the score matrix:')
        printMatrix(matrix)
        print('Printing the trace back matrix:')
        printMatrix(traceBack)
    if False:
        print('The globally aligned sequences are:')
        print(*xSeq[::-1])
        # Add pipes and crosses somehow.
        print(*ySeq[::-1])
    else:
        # aln_df = pd.DataFrame(data=[xSeq,ySeq],columns=['base','new'])
        aln_df = pd.DataFrame(data=[xSeq[::-1],ySeq[::-1]])
        aln_df = aln_df.transpose()
        # aln_df[2] = aln_df[0].equals(aln_df[1])
        aln_df[2] = np.where(aln_df[0] == aln_df[1], '=', 'x')
        aln_df = aln_df[[0,2,1]]
        aln_df = aln_df.rename(columns={0: 'base', 2:'match', 1:'new'})

        # Add pipes and crosses somehow.

        print('The globally aligned sequences are:')
        print(aln_df)

# driver('aaac', 'agc', -2, 1, -1, debug=True)
driver(['a','a','a','c'], ['a','g','c'], -2, 1, -1, debug=True)