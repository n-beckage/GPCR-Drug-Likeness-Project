import numpy as np
from rdkit import Chem

mil_smiles=np.load("numpy_objs/selected_1mil_smiles.npy",allow_pickle=True)
print(type(mil_smiles))
print(len(mil_smiles))

mil_mols=[]
for i in np.arange(len(mil_smiles)):
	# Need to add a filter to remove None types from mil_mols
	# B/c Chem.MolFromSmiles is sometimes inadvertently converting some of the smiles to None objects
	# And that's causing the error in the server job
	mil_mols.append(Chem.MolFromSmiles(mil_smiles[i]))

mil_mols=[i for i in mil_mols if i]

print(type(mil_mols))
print(len(mil_mols))
print(mil_mols[0])