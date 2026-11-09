from dist_matrix import *
from upgma_global_alignment import *

# to run upgma.c
import subprocess
import os

# to convert node ids to names in the Newick tree
import re

# for ascii tree visualisation
from Bio import Phylo
from io import StringIO


# TO DO
# implement error handling!!
# comment

def run_upgma(matrix):
    try:
        subprocess.run(["gcc", "upgma.c", "-o", "upgma.exe"], capture_output = True, check = True)
        
        result = subprocess.run(["./upgma.exe"] + matrix, capture_output = True, text = True)
        print(result.stderr)
        
        print("----")
        newick = result.stdout.strip()
        return newick
    except subprocess.CalledProcessError as error:
        print("Failed due to following error: ")
        print(error.stderr)
    
    finally:
        if os.path.exists("upgma.exe"):
            os.remove("upgma.exe")
    

def add_names_to_newick(newick, records):
    named_newick = newick
    
    for index, (name, seq) in records.items():
        named_newick = re.sub(rf'(?<=[(,]){index}(?=:)', name, named_newick)
        # magic!!
        # lookbehind - thing ? before <= the match must be ( or ,
        # index you're looking for
        # lookforward - thing ? after match must be :
        
    return named_newick
    
if __name__ == "__main__":
    fasta = input("Fasta File Path: ")
    print("""
        
        ============================== CREATING DISTANCE MATRIX ============================== 
    
    """)
    matrix, records = main_dist_matrix(fasta)
    print("""
        
        ============================== RUNNING UPGMA ============================== 
    
    """)
    newick = run_upgma(matrix)
    named_newick = add_names_to_newick(newick, records)
    print("")
    print(named_newick)
    print("")
    Phylo.draw_ascii(Phylo.read(StringIO(named_newick), "newick"))