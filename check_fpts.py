import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

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
        # NOTE: some objects will be removed by this filter, so fpt_array and smile_list will be smaller than molList
        if moli_withH.GetNumAtoms()<=80 and all([x.GetAtomicNum() in pos_at for x in moli_withH.GetAtoms()]) and all([x.GetFormalCharge() in pos_charges for x in moli_withH.GetAtoms()]):
            # creating an empty numpy array
            arr = np.zeros((0,), dtype=np.float32)
            # Creating the fingerprints and also converting them to a numpy array, storing in arr5
            DataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(moli_withH,4, nBits=1024),arr)
            # adding the fingerprints (as np arrays) to list
            fpt_arr.append(arr)
            # creating a corresponding list of smile strings for the mols
            smile_list.append(smile)
    return (fpt_arr,smile_list)

# The smile string from Jacob
first_smil="COC(=O)CC1CCCCN1C(=O)c1ccc(F)cc1F"
second_smil="CCN(CCC(=O)OC)C(=O)Cc1cn2ccsc2n1"

# Converting the smile to a mol
first_mol=Chem.MolFromSmiles(first_smil)
second_mol=Chem.MolFromSmiles(second_smil)

print(np.argwhere(fingerprints[0]==1).T)
print(np.argwhere(fingerprints[1]==1).T)