import pandas as pd

sequence1 = input("Sequence 1: ")
sequence2 = input("Sequence 2: ")

residues1 = ['gap1'] + list(sequence1)
residues2 = ['gap2'] + list(sequence2)

df = pd.DataFrame(columns = residues1, index = residues2)

print(df)

def cell_score():
    pass

for j in range(len(residues2)):
    for i in range(len(residues1)):
        pass
    