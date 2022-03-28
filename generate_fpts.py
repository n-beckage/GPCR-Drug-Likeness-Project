# This script will generate the fingerprints for the mol objects
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

# loading in the array of mol objects
theMols=np.load("theMols.npy",allow_pickle=True)

######################## generate_ftps() ##################################
#INPUT:
# an array or list of mol objects
#OUTPUT:
# a list of fingerprints (each fingerprint is an array)
# a list of corresponding SMILE strings (a list of char strings)
def generate_fpts(molList):
    pos_at=[1,6,7,8,9,11,15,16,17,19,35,53]
    pos_charges=[-1,0,1]
    # the array for the fingerprints
    fpt_arr=[]
    # The corresponding array for the smile strings
    smile_list=[]
    # mol_withH coresponds to moll, also corresponds to moli
    # Looks like you will need to iterate through molList witha for loop and call Chem.addHs(molList[i]), for instance
    for i in np.arange(len(molList)): # iterate through molList
        moli_withH=Chem.AddHs(molList[i])
        smile=Chem.MolToSmiles(moli_withH)
        if moli_withH.GetNumAtoms()<=80 and all([x.GetAtomicNum() in pos_at for x in moli_withH.GetAtoms()]) and all([x.GetFormalCharge() in pos_charges for x in moli_withH.GetAtoms()]):
            # creating an empty numpy array
            arr = np.zeros((0,), dtype=np.float32)
            # Creating the fingerprints and also converting them to a numpy array, storing in arr5
            DataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(moli_withH,4, nBits=1024),arr)
            # adding the fingerprints (as np arrays) to list
            fpt_arr.append(arr)
            # creating a corresponding list of smile strings for the mols
            smile_list.append(smile)
        else: 
            fpt_arr.append(None)
            smile_list.append(None)
    return (fpt_arr,smile_list)


# running the function
fingerprints, smiles = generate_fpts(theMols)

# saving the objects for later use
# np.save("fingerprints.npy",fingerprints,allow_pickle=True)
# np.save("smiles.npy",smiles,allow_pickle=True)
