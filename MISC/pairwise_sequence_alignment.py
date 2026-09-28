import pandas as pd

match_score = 1
mismatch_score = -1
gap_score = -2

sequence1 = input("Sequence 1: ")
sequence2 = input("Sequence 2: ")

residues1 = ['gap1'] + list(sequence1)
residues2 = ['gap2'] + list(sequence2)

df = pd.DataFrame(columns = residues1, index = residues2)

def cell_score(i, j, Sgh, Sgj, Sih): #where g refers to position i-1 and h refers to position j-1
    if residues2[i] == residues1[j]:
        sij = match_score
    else:
        sij = mismatch_score
    
    score1 = Sgh + sij
    score2 = Sgj - gap_score
    score3 = Sih - gap_score
    
    final_cell_score = max(score1, score2, score3)
    return final_cell_score 


# Fill in cells
for j in range(len(residues2)):
    if j == 0: # initialise first column
        for i in range(len(residues1)):
            df.loc[residues2[i], residues1[j]] = gap_score * i
    # fill in rest of cells
    else:
        for i in range(len(residues1)):
            if i == 0:
                df.loc[residues2[i], residues1[j]] = gap_score * j
            else:
                df.loc[residues2[i], residues1[j]] = cell_score(i - 1, j - 1, df.loc[residues2[i-1], residues1[j-1]], df.loc[residues2[i-1], residues1[j-1]], df.loc[residues2[i], residues1[j]])
    
print(df)