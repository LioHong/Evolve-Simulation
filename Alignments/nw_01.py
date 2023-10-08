import numpy as np
import pandas as pd
# https://gist.github.com/slowkow/06c6dba9180d013dfd82bec217d22eb5
def nw(x, y, match = 1, mismatch = 1, gap = 1):
    
    nx = len(x)
    ny = len(y)
    # Optimal score at each possible pair of characters.
    F = np.zeros((nx + 1, ny + 1))
    F[:,0] = np.linspace(0, -nx * gap, nx + 1)
    F[0,:] = np.linspace(0, -ny * gap, ny + 1)
    # Pointers to trace through an optimal aligment.
    P = np.zeros((nx + 1, ny + 1))
    P[:,0] = 3
    P[0,:] = 4
    # Temporary scores.
    t = np.zeros(3)
    for i in range(nx):
        for j in range(ny):
            if x[i] == y[j]:
                t[0] = F[i,j] + match
            else:
                t[0] = F[i,j] - mismatch
            t[1] = F[i,j+1] - gap
            t[2] = F[i+1,j] - gap
            tmax = np.max(t)
            F[i+1,j+1] = tmax
            if t[0] == tmax:
                P[i+1,j+1] += 2
            if t[1] == tmax:
                P[i+1,j+1] += 3
            if t[2] == tmax:
                P[i+1,j+1] += 4
    # Trace through an optimal alignment.
    i = nx
    j = ny
    print(i)
    print(j)
    rx = []
    ry = []
    while i > 0 or j > 0:
        if P[i,j] in [2, 5, 6, 9]:
            rx.append(x[i-1])
            ry.append(y[j-1])
            i -= 1
            j -= 1
        elif P[i,j] in [3, 5, 7, 9]:
            rx.append(x[i-1])
            ry.append('-')
            i -= 1
        elif P[i,j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(y[j-1])
            j -= 1
    # # Reverse the strings.
    # rx = ''.join(rx[::-1])
    # ry = ''.join(ry[::-1])
    # return '\n'.join([rx, ry])


    def pandalign(xSeq,ySeq):
        # aln_df = pd.DataFrame(data=[xSeq,ySeq],columns=['base','new'])
        aln_df = pd.DataFrame(data=[xSeq[::-1],ySeq[::-1]])
        aln_df = aln_df.transpose()
        # aln_df[2] = aln_df[0].equals(aln_df[1])
        aln_df[2] = np.where(aln_df[0] == aln_df[1], '-', '.')
        aln_df = aln_df[[0,2,1]]
        simil = (aln_df[2] == '-').sum() / len(aln_df) * 100
        aln_df = aln_df.rename(columns={2:'d'})
        # Add pipes and crosses somehow.
        print('Similarity: ' + "{:3.1f}".format(simil) + '%')
        print('The globally aligned sequences are:')
        return aln_df

    aln_df = pandalign(rx,ry)
    print(aln_df)

# x = "GATTACA"
# y = "GCATGCU"
# Why isn't the gap placed at the end?
x = ['G','Ax','T','T','A','C','A','T']
y = ['G','C','A','T','G','C','U']
# x = ['Interestingly','he','is','also','doing','a','special','issue','on','Taylor','and','Francis','journal','at','the','same','time']
# y = ['Interestingly','he','is','also','editing','a','special','issue','on','Taylor','and','Francis','journal','at','the','same','time']
print(nw(x, y))
# G-ATTACA
# GCA-TGCU
