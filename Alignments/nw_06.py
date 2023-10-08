# https://studentsxstudents.com/coding-global-sequence-alignment-using-the-needleman-wunsch-algorithm-d47971ebbe5
# https://github.com/karthikm15/needleman-wunsch-and-blast/blob/main/Needleman-Wunsch%20Algorithm.ipynb
# Initializing The Matrix
def nw(x, y, match = 1, mismatch = 1, gap = 2):
    nx = len(x) # length of first base sequence
    ny = len(y) # length of second base sequence
    
    # Initialization process - forming the base matrix
    F = np.zeros((nx + 1, ny + 1)) # an array 
    F[:,0] = np.linspace(0, -gap*nx, nx + 1)
    F[0,:] = np.linspace(0, -gap*ny, ny + 1)
    
    # Pointers to trace through an optimal aligment.
    P = np.zeros((nx + 1, ny + 1))
    P[:,0] = 3
    P[0,:] = 4

# Filling the Matrices
t = np.zeros(3)
for i in range(nx):
  for j in range(ny):
     # Iteration step: take the max (inserting gap in first sequence, inserting gap in second sequence, match or mutation)
       if x[i] == y[j]:
           t[0] = F[i,j] + match
       else:
            t[0] = F[i,j] - mismatch
            
        # Inserting gap in first sequence
        t[1] = F[i,j+1] - gap
        # Inserting gap in second sequence
        t[2] = F[i+1,j] - gap
        tmax = np.max(t)
            
        F[i+1,j+1] = tmax
        if t[0] == tmax:
            P[i+1,j+1] += 2
                
        # Higher weights for inserting gaps rather than matches/mismatches
        if t[1] == tmax:
           P[i+1,j+1] += 3
        if t[2] == tmax:
           P[i+1,j+1] += 4

# Checking Pointers to Find Optimal Arrangement
i = nx
j = ny
rx = []
ry = []
tracer_matrix = np.zeros((nx+1, ny+1))    
while i > 0 or j > 0:
    tracer_matrix[i, j] = -1
    # if there is a match/mismatch    
    if P[i,j] in [2, 5, 6, 9]:
        rx.append(x[i-1])
        ry.append(y[j-1])
            
        i -= 1
        j -= 1
        
    # if there's a gap in the first sequence
    elif P[i,j] in [3, 7]:
        rx.append(x[i-1])
        ry.append('-')
        i -= 1
            
    # if there's a gap in the second sequence
    elif P[i,j] in [4]:
        rx.append('-')
        ry.append(y[j-1])
        j -= 1

# Displaying Results
def nx(x,y,mismatch,match,gap):
    def add_bases (x, y, m):
        matrix = m.tolist()
        matrix.insert(0, [" ",  " "]+list(y))
        for i in range(2, len(x) + 2):
            matrix[i] = [list(x)[i-2]] + matrix[i]
        matrix[1] = [" "] + matrix[1]
        return matrix
    
    def print_matrix(matrix):
        s = [[str(e) for e in row] for row in matrix]
        lens = [max(map(len, col)) for col in zip(*s)]
        fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
        table = [fmt.format(*row) for row in s]
        print('\n'.join(table))
    print("Pointer matrix:")
    pointer_matrix = add_bases(x, y, P)
    print_matrix(pointer_matrix)
    
    print()
    print("Tracer matrix:")
    tracer_matrix = add_bases(x, y, tracer_matrix)
    print_matrix(tracer_matrix)
    
    # Reverse the strings.
    print()
    print("Final result:")
    
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    
    px = "Sequence 1: " + rx
    py = "Sequence 2: " + ry
    return ['\n'.join([px, py]), rx, ry]

# Applications in E. Coli Data
dna = open("e_coli_seq.txt")
# Data structuring off of how the text file is built
seq_arrange = []
for line in dna:
    if line[0] != '>':
        end = len(line)
        seq_arrange.append(line[:-1])
seq_arrange = [i[:13] for i in seq_arrange]
# Choosing random integers for the E. coli sequences
first_index = 0
second_index = 0
while (first_index == second_index):
    first_index = randint(0, len(seq_arrange))
    second_index = randint(0, len(seq_arrange))
first_seq = seq_arrange[first_index]
second_seq = seq_arrange[second_index]
printseq, seq1, seq2 = nw(first_seq, second_seq)
print(printseq)
