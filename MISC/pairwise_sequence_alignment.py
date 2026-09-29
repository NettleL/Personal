import pandas as pd

# constants - maybe later add the ability for the user to change these constants??
MATCH_SCORE = 1
MISMATCH_SCORE = -1
GAP_SCORE = -2

def cell_score(residue1, residue2, Score_gh, Score_gj, Score_ih):
    #where g refers to position i-1 and h refers to position j-1
    
    trace_arrow = ''
    
    if residue1 == residue2:
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
        trace_arrow += 'u'
    if final_cell_score == score3:
        trace_arrow += 'l'
  
    return final_cell_score, trace_arrow

def form_matrix():
    sequence1 = input("Sequence 1: ")
    sequence2 = input("Sequence 2: ")
    
    residues1 = ['gap1'] + list(sequence1)
    residues2 = ['gap2'] + list(sequence2)
        
    df = pd.DataFrame(index = range(len(residues2)), columns = range(len(residues1))) # so index + column vals are unique
    traceback_df = pd.DataFrame(index = range(len(residues2)), columns = range(len(residues1)))
    
    # Fill in cells
    for j in range(len(residues1)):
        if j == 0: # initialise first column
            for i in range(len(residues2)):
                df.loc[i, j] = GAP_SCORE * i
        else:
            for i in range(len(residues2)):
                # initialise first row
                if i == 0:
                    df.loc[i, j] = GAP_SCORE * j
                # fill in rest of cells
                else:
                    score, trace_value = cell_score(residues1[j],
                                                    residues2[i],
                                                    df.loc[i-1, j-1],
                                                    df.loc[i-1, j],
                                                    df.loc[i, j-1])
                    df.loc[i, j] = score
                    traceback_df.loc[i, j] = trace_value
    return residues1, residues2, df, traceback_df


def total_traceback(traceback_matrix):
    # note - will never be dead ends when you traceback (i think) so don't worry about that
    divergence_points = []
    
    # dont use dictionary for storing tracebacks
    final_paths = []
    
    # make list of divergence points
    for i in range(traceback_matrix.shape[0]): # rows
        for j in range(traceback_matrix.shape[1]): # columns
            if traceback_matrix.iloc[i, j] in ['du', 'dl', 'ul', 'dul']:
                divergence_points.append((i, j))
    
    for i, j in divergence_points:
        for arrow in traceback_matrix.iloc[i,j]:
            matrix_copy = traceback_matrix.copy()
            matrix_copy.iloc[i, j] = arrow
            path = single_traceback(matrix_copy)
            final_paths.append(path)
            
    unique_final_paths = []
    paths_added_to_unique_final_paths = [] # because paths are tuples weird things happen :((
    
    for path in final_paths:
        tupled_path = tuple(path)
        if tupled_path not in paths_added_to_unique_final_paths:
            unique_final_paths.append(path)
            paths_added_to_unique_final_paths.append(tupled_path)
    
    return unique_final_paths
    
def single_traceback(traceback_matrix):
    i = traceback_matrix.shape[0] - 1
    j = traceback_matrix.shape[1] - 1
    
    path = []
    
    while i > 0 or j > 0:
        arrows = traceback_matrix.iloc[i, j]
        
        arrow_used = arrows[0] # default to first arrow if arrow is not specified. do not worry!! each divergence point (i.e. each arrow) will have a pass dedicated to it!! so no divergence points will be missed, even though the first arrow is the default
        
        path.append((i, j, arrow_used))
        if arrow_used == 'd':
            i -= 1
            j -= 1
        elif arrow_used == 'u':
            i -= 1
        elif arrow_used == 'l':
            j -= 1
        
    path.append((0, 0, None)) # add first cell to end of path
    
    chronological_path = path[::-1] # switches path from backwards to forwards
    
    return chronological_path

def main():
    residues1, residues2, alignment_matrix, traceback_matrix = form_matrix()
    
    print(alignment_matrix)
    print(traceback_matrix)
    
    optimal_alignment = alignment_matrix.iloc[-1, -1]
    
    print(total_traceback(traceback_matrix))

main()
