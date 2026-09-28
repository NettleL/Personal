import pandas as pd

match_score = 1
mismatch_score = -1
gap_score = -2

def cell_score(residue2, residue1, Score_gh, Score_gj, Score_ih):
    #where g refers to position i-1 and h refers to position j-1
    
    trace_arrow = ''
    
    if residue2 == residue1:
        sij = match_score
    else:
        sij = mismatch_score

    score1 = Score_gh + sij
    score2 = Score_gj - gap_score
    score3 = Score_ih - gap_score

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
        
    df = pd.DataFrame(columns = residues1, index = residues2)
    traceback_df = pd.DataFrame(columns = residues1, index = residues2)
    
    # Fill in cells
    for j in range(len(residues2)):
        if j == 0: # initialise first column
            for i in range(len(residues1)):
                df.loc[residues2[i], residues1[j]] = gap_score * i
        else:
            for i in range(len(residues1)):
                # initialise first row
                if i == 0:
                    df.loc[residues2[i], residues1[j]] = gap_score * j
                # fill in rest of cells
                else:
                    score, trace_value = cell_score(residues2[i-1],
                                                    residues1[j-1],
                                                    df.loc[residues2[i-1], residues1[j-1]],
                                                    df.loc[residues2[i-1], residues1[j]],
                                                    df.loc[residues2[i], residues1[j-1]])
                    df.loc[residues2[i], residues1[j]] = score
                    traceback_df.loc[residues2[i], residues1[j]] = trace_value
    return df, traceback_df

def traceback():
    pass

def main():
    alignment_matrix, traceback_matrix = form_matrix()
    
    print(alignment_matrix)
    print(traceback_matrix)
    
    optimal_alignment = alignment_matrix.iloc[-1, -1]

main()