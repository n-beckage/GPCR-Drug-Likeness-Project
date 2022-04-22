from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import numpy as np
import pandas as pd

df=pd.read_pickle("numpy_objs/framed_data.zip")

############ Print images of the Lipinksi passes with the highest activity scores ########
# Grouping the data by Lipinski pass or fail
splitLip = df.groupby("Lipinski")

# Getting the grouped dataframes
passLip = splitLip.get_group("Druglike")
failLip = splitLip.get_group("Not Druglike")

# Getting the top 5 active drug-like molecules
top_5_active=passLip.nlargest(5,columns="GPCR_act")
print(top_5_active[['SMILES','GPCR_act']])

# Creating a bond-line diagram of these top five (1=most active, 5=least active)
for i in range(0,len(top_5_active)):
    mol=Chem.MolFromSmiles(top_5_active.iloc[i,0])
    AllChem.Compute2DCoords(mol)
    Draw.MolToFile(mol,"plots/{}_active_LipPass.png".format(i+1))
    print("{}_active_LipPass\n".format(i+1),top_5_active.iloc[i,],"\n")

########### Print Lipinski passes with the lowest activity scores ###############
# Getting the 5 Lip passes with the lowest activity
low_5_active=passLip.nsmallest(5,columns="GPCR_act")
print(low_5_active[['SMILES','GPCR_act']])

for i in range(0,len(low_5_active)):
    mol=Chem.MolFromSmiles(low_5_active.iloc[i,0])
    AllChem.Compute2DCoords(mol)
    Draw.MolToFile(mol,"plots/{}_inactive_LipPass.png".format(i+1))
    print("{}_inactive_LipPass\n".format(i+1),low_5_active.iloc[i,],"\n")

########### Print Lipinski Fails with the highest activity scores ###############
# Getting the 5 Lip fails with the highest activity
top_5_active=failLip.nlargest(5,columns="GPCR_act")
print(top_5_active[['SMILES','GPCR_act']])

for i in range(0,len(top_5_active)):
    mol=Chem.MolFromSmiles(top_5_active.iloc[i,0])
    AllChem.Compute2DCoords(mol)
    Draw.MolToFile(mol,"plots/{}_active_LipFail.png".format(i+1))
    print("{}_active_LipFail\n".format(i+1),top_5_active.iloc[i,],"\n")

########### Print Lipinski Fails with the lowest activity scores ###############
# Getting the 5 Lip fails with the lowest activity
low_5_active=failLip.nsmallest(5,columns="GPCR_act")
print(low_5_active[['SMILES','GPCR_act']])

for i in range(0,len(low_5_active)):
    mol=Chem.MolFromSmiles(low_5_active.iloc[i,0])
    AllChem.Compute2DCoords(mol)
    Draw.MolToFile(mol,"plots/{}_inactive_LipFail.png".format(i+1))
    print("{}_inactive_LipFail\n".format(i+1),low_5_active.iloc[i,],"\n")

############ Finding the molecules Jacob sent me ###################

# # The smile string from Jacob
# first_smil="COC(=O)CC1CCCCN1C(=O)c1ccc(F)cc1F"
# second_smil="CCN(CCC(=O)OC)C(=O)Cc1cn2ccsc2n1"

# # Converting the smile to a mol
# first_mol=Chem.MolFromSmiles(first_smil)
# second_mol=Chem.MolFromSmiles(second_smil)

# # Adding H's to that mol
# first_mol=Chem.AddHs(first_mol)
# second_mol=Chem.AddHs(second_mol)

# first_smil=Chem.MolToSmiles(first_mol)
# second_smil=Chem.MolToSmiles(second_mol)


# # Finding the smile string in the dataframe
# smil1_row=df[df['SMILES']==first_smil]
# print(first_smil)
# print(smil1_row)
# if smil1_row.shape[0] == 0:
#     print("No matching smile string in the dataframe for the first given string")
#     #exit()
# smil2_row=df[df['SMILES']==second_smil]
# if smil1_row.shape[0] == 0:
#     print("No matching smile string in the dataframe for the second given string")
#     #exit()

# # Printing the smile string, activity, and Lipinski results
# print("FIRST SMILE\n",smil1_row[["SMILES", "GPCR_act", "Lipinski"]])
# print("SECOND SMILE\n",smil2_row[["SMILES", "GPCR_act", "Lipinski"]])

# ### It seems that these smiles are not in the dataframe... maybe we could look in the 

# # Drawing the images of these molecules
# ## First we have to remove H's so the drawing looks OK
# first_mol=Chem.RemoveHs(first_mol)
# second_mol=Chem.RemoveHs(second_mol)
# AllChem.Compute2DCoords(first_mol)
# AllChem.Compute2DCoords(second_mol)
# Draw.MolToFile(first_mol, 'plots/first_mol.png')
# Draw.MolToFile(second_mol, 'plots/second_mol.png')

# fpts=np.load("numpy_objs/fingerprints.npy",allow_pickle=True)
# print(np.argwhere(fpts==1).T)