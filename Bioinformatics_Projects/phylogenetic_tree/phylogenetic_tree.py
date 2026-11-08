from dist_matrix import *
from upgma_global_alignment import *
import subprocess
import os

# for ascii tree visualisation
from Bio import Phylo
from io import StringIO

# TO DO
# map node ids to names
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
    

if __name__ == "__main__":
    fasta = input("Fasta File Path: ")
    print("""
        
        ============================== CREATING DISTANCE MATRIX ============================== 
    
    """)
    matrix = main_dist_matrix(fasta)
    print("""
        
        ============================== RUNNING UPGMA ============================== 
    
    """)
    newick = run_upgma(matrix)
    print(newick)
    
    Phylo.draw_ascii(Phylo.read(StringIO(newick), "newick"))