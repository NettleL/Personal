import pandas as pd

# constants - maybe later add the ability for the user to change these constants??
MATCH_SCORE = 1
MISMATCH_SCORE = -1
GAP_SCORE = -2

def cell_score(residue2, residue1, Score_gh, Score_gj, Score_ih):
    #where g refers to position i-1 and h refers to position j-1
    
    trace_arrow = ''
    
    if residue2 == residue1:
        sij = MATCH_SCORE
    else:
        sij = MISMATCH_SCORE

    score1 = Score_gh + sij
    score2 = Score_gj + GAP_SCORE
    score3 = Score_ih + GAP_SCORE

    final_cell_score = max(score1, score2, score3)
    
    if final_cell_score == score1:
        trace_arrow += 'd'
    if final_cell_score == score2:
        trace_arrow += 'l'
    if final_cell_score == score3:
        trace_arrow += 'u'
  
    return final_cell_score, trace_arrow

def form_matrix():
    sequence1 = input("Sequence 1: ")
    sequence2 = input("Sequence 2: ")
    
    residues1 = ['gap1'] + list(sequence1)
    residues2 = ['gap2'] + list(sequence2)
        
    df = pd.DataFrame(index = range(len(residues1)), columns = range(len(residues2))) # so index + column vals are unique
    traceback_df = pd.DataFrame(index = range(len(residues1)), columns = range(len(residues2)))
    
    # Fill in cells
    for j in range(len(residues2)):
        if j == 0: # initialise first column
            for i in range(len(residues1)):
                df.loc[i, j] = GAP_SCORE * i
        else:
            for i in range(len(residues1)):
                # initialise first row
                if i == 0:
                    df.loc[i, j] = GAP_SCORE * j
                # fill in rest of cells
                else:
                    score, trace_value = cell_score(residues2[j],
                                                    residues1[i],
                                                    df.loc[i-1, j-1],
                                                    df.loc[i-1, j],
                                                    df.loc[i, j-1])
                    df.loc[i, j] = score
                    traceback_df.loc[i, j] = trace_value
    return df, traceback_df

def traceback():
    pass

def main():
    alignment_matrix, traceback_matrix = form_matrix()
    
    print(alignment_matrix)
    print(traceback_matrix)
    
    optimal_alignment = alignment_matrix.iloc[-1, -1]

main()

''' returns this - 
Sequence 1: taat
Sequence 2: tcgat
    0   1   2   3   4    5
0   0  -2  -4  -6  -8  -10
1  -2   1   3   5   7    9
2  -4   3   5   7   9   11
3  -6   5   7   9  11   13
4  -8   7   9  11  13   15
     0    1    2    3    4    5
0  NaN  NaN  NaN  NaN  NaN  NaN
1  NaN    d    u    u    u    u
2  NaN    l   lu   lu   lu   lu
3  NaN    l   lu   lu   lu   lu
4  NaN    l   lu   lu   lu   lu
'''
