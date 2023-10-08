# Updated print to print() but the prints are 1 char/line...
"""
align.py
 
Contains an implementation of the Needleman Wunsch algorithm.
Written by Michael Hamilton, and edited by Asa Ben-Hur
"""
 
import numpy as np
from random import randint
 
verbose = True  # in verbose mode the DP and traceback will be printed
 
# stuff for displaying the traceback matrix:
 
# to print colored arrows you will need the termcolor module
# available at https://pypi.python.org/pypi/termcolor
# if you don"t have it, arrows will be printed without color
color = True
try :
    from termcolor import colored
except :
    color = False
 
# the three directions you can go in the traceback:
DIAG = 0 
UP = 1 
LEFT = 2
# UTF-8 representations of arrow symbols
# arrows[DIAG] is a diagonal arrow
# arrows[UP] is an up arrow
# arrows[LEFT] is a left
arrows = [u"\u2196", u"\u2191", u"\u2190"]
 
 
def show_ptr_matrix(ptr, seq1, seq2) :
 
    print("\n"+"~`"*25)
    print("Traceback")
    global color
    print("       " + "  ".join(seq2))
    for i in range(len(ptr)) :
        if (i > 0) :
            print(seq1[i-1] + " ",)
        if (i == 0) :
            print("  ",)
        for j in range(len(ptr[i])) :
            if color and ptr[i,j] >= 3 :
                print(" " + colored(arrows[ptr[i,j]-3], "green" ),)
            else :
                if ptr[i,j] >=3 :
                    ptr[i,j] -=3
                print(" " + arrows[ptr[i,j]],)
        print()
 
def show_dp_matrix(s, seq1, seq2) :
 
    print("\n"+"~`"*25)
    print("DP matrix")
    print("            " + "     ".join(seq2))
    for i in range(len(s)) :
        if (i > 0) :
            print(seq1[i-1] + " ",)
        if (i == 0) :
            print("  ",)
        for j in range(len(s[i])) :
            print(" " + "% 2.1f" % s[i,j],)
        print()
 
def needleman_wunsch_matrix(seq1, seq2):
    """
    fill in the DP matrix according to the Needleman-Wunsch algorithm.
    Returns the matrix of scores and the matrix of pointers
    """
 
    match =  1  # match score
    mismatch = -1  # mismatch penalty
    indel = -1 # indel penalty
 
    n = len(seq1)
    m = len(seq2)
    s = np.zeros( (n+1, m+1) ) # DP matrix
    ptr = np.zeros( (n+1, m+1), dtype=int  ) # matrix of pointers
 
    ##### INITIALIZE SCORING MATRIX (base case) #####
 
    for i in range(1, n+1) :
        s[i,0] = indel * i
    for j in range(1, m+1):
        s[0,j] = indel * j
 
    ########## INITIALIZE TRACEBACK MATRIX ##########
 
    # Tag first row by LEFT, indicating initial "-"s
    ptr[0,1:] = LEFT
 
    # Tag first column by UP, indicating initial "-"s
    ptr[1:,0] = UP
 
    #####################################################
 
    for i in range(1,n+1):
        for j in range(1,m+1): 
            # match
            if seq1[i-1] == seq2[j-1]:
                s[i,j] = s[i-1,j-1] + match
                ptr[i,j] = DIAG
            # mismatch
            else :
                s[i,j] = s[i-1,j-1] + mismatch
                ptr[i,j] = DIAG
            # indel penalty
            if s[i-1,j] + indel > s[i,j] :
                s[i,j] = s[i-1,j] + indel
                ptr[i,j] = UP
            # indel penalty
            if s[i, j-1] + indel > s[i,j]:
                s[i,j] = s[i, j-1] + indel
                ptr[i,j] = LEFT
 
    return s, ptr
 
def needleman_wunsch_trace(seq1, seq2, s, ptr) :
 
    #### TRACE BEST PATH TO GET ALIGNMENT ####
    align1 = ""
    align2 = ""
    n, m = (len(seq1), len(seq2))
    i = n
    j = m
    curr = ptr[i, j]
    while (i > 0 or j > 0):        
        ptr[i,j] += 3
        if curr == DIAG :            
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1            
        elif curr == LEFT:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1            
        elif curr == UP:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
 
        curr = ptr[i,j]
 
    return align1, align2
 
def needleman_wunsch(seq1, seq2) :
    """
    computes an optimal global alignment of two sequences using the Needleman-Wunsch
    algorithm
    returns the alignment and its score
    """
    s,ptr = needleman_wunsch_matrix(seq1, seq2)
    alignment = needleman_wunsch_trace(seq1, seq2, s, ptr)
 
    if verbose :
        show_dp_matrix(s, seq1, seq2)
        show_ptr_matrix(ptr, seq1, seq2)
        print("\n"+"~`"*25)
        print("Alignment Score: %f\n" % (s[len(seq1),len(seq2)]))
        print("Alignment:")
        print(alignment[0])
        print(alignment[1])
 
    return alignment, s[len(seq1), len(seq2)]
 
def random_DNA_sequence(length):
    """
    Returns a random DNA of the given length.
    """
    nucleotides = ["A","T","G","C"]
    seq = [ nucleotides[randint(0,3)] for i in range(length) ]
    return "".join(seq)
 
 
if __name__ == "__main__":
 
    seq1 = "ATCGA"
    seq2 = "ATGA"
    needleman_wunsch(seq1, seq2)
 
    seq1 = "ATGTTCGA"
    seq2 = "GGATTACTGA"
    needleman_wunsch(seq1, seq2)
 
    seq1 = "TCCCAGTTATGTCAGGGGACACGAGCATGCAGAGAC"
    seq2 = "AATTGCCGCCGTCGTTTTCAGCAGTTATGTCAGATC"
    needleman_wunsch(seq1, seq2)
 
    seq1 = random_DNA_sequence(20)
    seq2 = random_DNA_sequence(20)
    needleman_wunsch(seq1, seq2)
