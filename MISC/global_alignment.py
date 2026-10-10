import pandas as pd

# AESTHETICS
from colorama import Fore
import matplotlib.pyplot as plt
import numpy as np

# constants
# note: maybe later add the ability for the user to change these constants??
MATCH_SCORE = 1
MISMATCH_SCORE = -1
GAP_SCORE = -2

def cell_score(residue1, residue2, Score_gh, Score_gj, Score_ih):
    '''
    Determines the score for each cell in the score matrix
    The score of cell_ij is dependent on the scores of the cells preceding it
    
    Note
        - s_ij (lowercase) = score for aligning only aᵢ with bⱼ (match or mismatch)
        - S_ij (uppercase) = score of optimal alignment of a₁..aᵢ with b₁..bⱼ -> i.e. a preceding cell
    
    Note on Parameters:
    (g refers to position i-1 and h refers to position j-1)
        - Score_gh: score of cell i-1, j-1
        - Score_gj: score of cell i-1, j
        - Score_ih: score of cell i, j-1
    
    Three scores are calculated and the maximum score is used
        1. Score_gh + s_ij (match/mismatch):            dependent on the score of the diagonal (upper left) cell (d)
        2. Score_gj + GAP_SCORE (gap in alignment 1):   dependent on the score of the upper cell (u)
        3. Score_ih + GAP_SCORE (gap in alignment 2):   dependent on the score of the left cell (l)
    
    Note that out of the three score calculations, the maximum score may come from multiple scores
    
    For a given cell, it returns the cell score and a list of the preceding cells (d, u, l) that were used to calculate it
    '''
    
    # calculating scores
    if residue1 == residue2:
        sij = MATCH_SCORE
    else:
        sij = MISMATCH_SCORE

    score1 = Score_gh + sij
    score2 = Score_gj + GAP_SCORE
    score3 = Score_ih + GAP_SCORE

    # determines the maximum score
    final_cell_score = max(score1, score2, score3)
    
    # initialises and fills a list the directions of the preceding cells (d, u, l)
    # that were used to calculate the cell score
    trace_arrow = []
    
    if final_cell_score == score1:
        trace_arrow.append('d') # diagonal
    if final_cell_score == score2:
        trace_arrow.append('u') # upper
    if final_cell_score == score3:
        trace_arrow.append('l') # left
  
    return final_cell_score, trace_arrow

def form_matrix(sequence1, sequence2):
    '''
    Given two sequences,
        1. Converts the sequences to a list of residues, adding a gap at the start
        2. Creates two dataframes
            a. a cell alignment score table
            b. a traceback table (which contains the directions of the preceding cells used to determine that cells score)
        3. Fills in the cells of both dataframes
        4. Returns both sequences and both dataframes
    '''
    
    residues1 = ['-'] + list(sequence1) # add a gap
    residues2 = ['-'] + list(sequence2) # add a gap
    
    # initialises a table of cell scores with a default score of 0
    # note: don't use residues - so index + column vals are unique
    df = pd.DataFrame(0, index = range(len(residues2)), columns = range(len(residues1)))
    
    # initialises a table of directions with a default value of []
    # note: magic :)) to initialise every cell as a lil list :(( want to find a better prettier way
    traceback_df = pd.DataFrame([[[] for j in range(len(residues1))] for i in range(len(residues2))],
                                    index = range(len(residues2)), columns = range(len(residues1)))

    # fill in cells
    for j in range(len(residues1)):
        for i in range(len(residues2)):
            # initialise table - if i or j = 0, use initialising values
            # needleman-wunsch - initialisation is dependent on GAP_SCORE
            if i == 0: # initialise first row
                df.loc[i, j] = GAP_SCORE * j
            elif j == 0: # initialise first column
                df.loc[i, j] = GAP_SCORE * i
            else: # determine value of cell ij in both tables using cell_score()
                score, trace_value = cell_score(residues1[j],
                                                residues2[i],
                                                df.loc[i-1, j-1],
                                                df.loc[i-1, j],
                                                df.loc[i, j-1])
                df.loc[i, j] = score
                traceback_df.loc[i, j] = trace_value
    return residues1, residues2, df, traceback_df

def total_traceback(traceback_matrix):
    '''
    Using the traceback matrix
    1. Determines the endpoint which the traceback starts from (last cell in table)
    2. Runs a total traceback (using single_traceback() function)
    and adds the traceback path to the list of final paths
    3. Returns the list of final paths
    '''
    
    # note: dont use dictionary for storing tracebacks
    final_paths = []
    
    # determines endpoint to start traceback from (needleman-wunsch - last cell)
    i = traceback_matrix.shape[0] - 1
    j = traceback_matrix.shape[1] - 1
    
    # runs a traceback
    paths = single_traceback(i, j, traceback_matrix)
    
    # adds reversed path to final_paths
    for path in paths:
        final_paths.extend([path[::-1]])
    
    return final_paths
    
def single_traceback(i, j, traceback_matrix, current_path=None):
    '''
    Recursively implements a traceback through the traceback_matrix to find all possible alignments.
    
    The function uses the values (d, u, l) of cell_ij in the traceback matrix.
    If this is the first call/use of the function (no current_path is provided), the current_path is initialised.
    Otherwise - i.e. when single_traceback() is recursively called within the function - the function uses the given current_path
    
    If the path is complete (needleman-wunsch - i, j  = 0, 0), the function returns the current_path as one finished traceback route.
    
    Otherwise, the function traces back (using directions d, u, l) to a previous cell,
    appends this movement to the current_path, and calls this function again with the new updated path.
    Therefore, if there are multiple arrows, the function recursively calls itself multiple times (once for each arrow),
    making multiple new paths (if these new paths have multiple arrows, the function recursively calls itself, making more paths e.tc.)
    until all possible paths have been traced. Once all the paths are finished (i, j = 0, 0), the function appends them to the list of paths.
    
    Therefore, while the function traces back one step at a time, when it is called by total_traceback(), it recursively calls itself multiple times
    until it and returns a list of all possible paths (where each path is a list of tuples).
    
    Note on variables:
    - current_path contains the traceback path of one single branch of the recursion
    - once the path is complete, it is returned (to single_traceback(), which called itself) and adds to paths
    - finally, the list of paths is returned to total_traceback()
    
    Note: Quite difficult :(( - never wrote this type of recursion before
    '''
    
    arrows = traceback_matrix.iloc[i, j] # retrieves arrows for a single cell (for a single traceback)

    # i.e. when the function is called by total_traceback - starts new empty current_path
    if current_path is None:
        current_path = []
    
    # if current_path is completed, return it to the instance of single_traceback() that called it    
    if i == 0 and j == 0:
        return [current_path + [(i, j, None)]] # note: current_path.append() returns None :((
    
    # note: recursion yippee :))) was a bit hard to do - never done this type of recursion before
    else: 
        paths = []
            
        for arrow in arrows: # for each arrow
            # creates new path in a separate list
            # each arrow in each recursive call has its own path
            new_path = current_path + [(i, j, arrow)]
            
            # function recursively calls itself, continuing the traceback from new_path
            # until finally this collectively returns a list of all complete paths
            # which is all added to paths
            if arrow == 'd' and i > 0 and j > 0:
                paths.extend(single_traceback(i-1, j-1, traceback_matrix, new_path))
            elif arrow == 'u' and i > 0:
                paths.extend(single_traceback(i-1, j, traceback_matrix, new_path))
            elif arrow == 'l' and j > 0:
                paths.extend(single_traceback(i, j-1, traceback_matrix, new_path))
        return paths # full list of paths
        
def form_alignment(residues1, residues2, unique_final_paths):
    alignments = []
    for path in unique_final_paths:
        alignment1 = ''
        alignment2 = ''
                
        for (i, j, step) in path:
            if step == 'd':
                alignment1 += residues1[j]
                alignment2 += residues2[i]
            if step == 'u':
                alignment1 += '-'
                alignment2 += residues2[i]
            if step == 'l':
                alignment1 += residues1[j]
                alignment2 += '-'
        alignments.append((alignment1, alignment2))
    return(alignments)

def deduplicate_alignments(alignments):
    seen = set() # tuples are weird
    unique_paths = []
    for sequence1, sequence2 in alignments:
        if (sequence1, sequence2) not in seen:
            seen.add((sequence1, sequence2))
            unique_paths.append((sequence1, sequence2))
    return unique_paths

def alignment_viewer(alignments):
    all_visualised_alignments = []
    for alignment in alignments:
        alignment1 = alignment[0]
        alignment2 = alignment[1]
        middle_row = ''
        for count in range(len(alignment1)): # for each position in the alignment
            if alignment1[count] == '-' or alignment2[count] == '-':
                middle_row += ' '
            elif alignment1[count] == alignment2[count]: # FOR AESTHETIC
                middle_row += '|'
            else:
                middle_row += Fore.RED + '|' + Fore.RESET
        visualised_alignment = "\n".join([alignment1, middle_row, alignment2])
        all_visualised_alignments.append(visualised_alignment)
    return all_visualised_alignments

# FOR AESTHETIC
def alignment_plot(alignment_matrix, paths):
    alignment_numpy = alignment_matrix.to_numpy(dtype=float)
    plt.matshow(alignment_numpy)
    
    for path in paths:
        x = []
        y = []
        for step in path:
            x.append(step[1]) # to match matrix coordinates - x = j (columns) and y = i (rows)
            y.append(step[0])
        plt.plot(x, y, 'k')
    plt.show()
    
        
def main():
    sequence1 = input("Sequence 1: ")
    sequence2 = input("Sequence 2: ")
    
    residues1, residues2, alignment_matrix, traceback_matrix = form_matrix(sequence1, sequence2)
    optimal_alignment = alignment_matrix.iloc[-1, -1]
    
    print('')
    print('ALIGNMENT MATRIX')
    print(alignment_matrix)
    print("Optimal Alignment Score =", optimal_alignment)
    print('')
    print('TRACEBACK MATRIX')
    print(traceback_matrix)
    print('')
    
    traceback = total_traceback(traceback_matrix)
    alignments = form_alignment(residues1, residues2, traceback)
    unique_alignments = deduplicate_alignments(alignments)
    visualised_alignments = alignment_viewer(unique_alignments)
    
    for alignment in visualised_alignments:
        print('ALIGNMENT')
        print(alignment)
        print('')

    alignment_plot(alignment_matrix, traceback)
        


main()
