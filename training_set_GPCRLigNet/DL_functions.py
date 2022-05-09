from rdkit.Chem import QED
from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import Crippen
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

############################### Lipinski Ro5 ##############################
# Binary function to see if mol passes Linpinski's rule of 5
# OUTPUT: 1 = druglike
#       : 0 = NOT druglike 
def Ro5(mol):
    strikes=0
    if Lipinski.NumHDonors(mol)>5:
        strikes+=1
    if Lipinski.NumHAcceptors(mol)>10:
        strikes=+1
    if Descriptors.ExactMolWt(mol)>499:
        strikes+=1
    if Descriptors.MolLogP(mol)>4.99:
        strikes+=1
    if strikes>1:
        return "Not Druglike"
    elif strikes<=1:
        return "Druglike"

#################################### Veber #########################################
# OUTPUT: 1 = druglike
#       : 0 = NOT druglike 
def veber(HBD,HBA,ROTB):
    if HBD+HBA<=12:
        if ROTB<=10:
            return "Druglike"
        else: return "Not Druglike"
    else: return "Not Druglike"

##################################### Ghose #########################################
# OUTPUT: 1 = druglike
#       : 0 = NOT druglike 
def ghose(MW, ALOGP, MR, ATOM):
    strike = 0
    if MW <= 160 or MW >= 480:
        strike = strike + 1
    if ALOGP <= -0.4 or ALOGP >= 5.6:
        strike = strike + 1
    if MR <= 40 or MR >= 130:
        strike = strike + 1
    if ATOM <= 20 or ATOM >= 70:
        strike = strike + 1
    if strike == 0:
        return "Druglike"
    else: 
        return "Not Druglike"