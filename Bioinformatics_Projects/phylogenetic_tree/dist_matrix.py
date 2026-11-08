from Bio import SeqIO
import numpy as np
import math
from upgma_global_alignment import *

# ADD ERROR HANDLING
def alignment_score(alignment, row, col):
    arbitrary_large_value = 10
    
    alignment_a = alignment[0]
    alignment_b = alignment[1]
    length = len(alignment_a) # is the same as len(alignment_b)
    raw_score = 0
    
    
    for element in range(length):
        if (alignment_a[element] != alignment_b[element] or
            alignment_a[element] == '-' or
            alignment_b[element] == '-'):
            raw_score += 1
    
    p_val = raw_score/length
    
    if p_val >= 0.75:
        print(f"Saturation detected between sequences {row} and {col} (p={p_val}). Distance capped at 10.0.")
        return arbitrary_large_value
    else:  
        jukes_cantor_score = -(3/4) * math.log(1 - (4/3) * (raw_score/length))
        return jukes_cantor_score

def fill_matrix(records, record_no):
    dist_matrix = np.full((record_no, record_no), -1.0)
    
    for row in range(record_no):
        for col in range(record_no):
            if row > col:
                alignment = main_global_alignment(False, records[row][1], records[col][1])
                score = alignment_score(alignment[0], row, col)
                dist_matrix[row][col] = score
                dist_matrix[col][row] = score
            elif row == col:
                score = 0
                dist_matrix[row][col] = score
                dist_matrix[col][row] = score
            
    return dist_matrix

def convert_matrix(dist_matrix, record_no):
    converted_matrix = [str(record_no)]
    
    for row in range(record_no):
        for col in range(record_no):
            converted_matrix.append(f"{dist_matrix[row][col]}")
            
    return converted_matrix
    
def main_dist_matrix(fasta_file):
    records = {}
    record_no = 0

    for record in SeqIO.parse(fasta_file, "fasta"):
        records[record_no] = (record.id, record.seq)
        record_no += 1
        
    dist_matrix = fill_matrix(records, record_no)
    converted_dist_matrix = convert_matrix(dist_matrix, record_no)
    
    print("------ DISTANCE MATRIX ------")
    print(dist_matrix)
    return converted_dist_matrix, records


        

