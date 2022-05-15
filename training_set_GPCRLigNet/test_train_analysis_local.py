# This script should run in the training_set_GPCRLigNet/ folder of Jacob's directory on the server
# The purpose of this script is to run the same analysis and make the same figures for the training/testing data as was done for the PubChem data

import numpy as np
import pandas as pd
from pathlib import Path
from DL_functions import *
from rdkit.Chem import QED
from rdkit import Chem

# defining the generic paths for molact
glass_molact_path="GLASS/molact_4_27_2022_"
dude_molact_path="DUDE/molact_4_27_2022_"

# defining the gneric paths for smiles
glass_smiles_path="GLASS/smiles_4_27_2022_"
dude_smiles_path="DUDE/smiles_4_27_2022_"

GLASS_molact=[]
GLASS_smiles=[]
#* DUDE_molact=[]
DUDE_smiles=[]

# Doing these outside of the loop because they don't fit the naming conventions of the other files
# thus it doesn't really make sense to include them in the loop
#* Now that I think about it, we don't need to read in DUDE mol act b/c they're all inactive
GLASS_molact.extend(np.load(glass_molact_path+"0.npy",allow_pickle=True))
GLASS_smiles.extend(np.load(glass_smiles_path+"0.npy",allow_pickle=True))
#* DUDE_molact.extend(np.load(dude_molact_path+"0.npy",allow_pickle=True))
DUDE_smiles.extend(np.load(dude_smiles_path+"0.npy",allow_pickle=True))

# These will tell us the last file that was read in each DUDE and GLASS folder
dude_count=0
glass_count=0
# start counter at 1 b/c we manually read in the 0 files
i=1
while Path(glass_molact_path+str(i)+"0000.npy").is_file(): # this works b/c GLASS has more files than DUDE
    if Path(dude_molact_path+str(i)+"0000.npy").is_file():
        #* DUDE_molact.extend(np.load(dude_molact_path+str(i)+"0000.npy",allow_pickle=True))
        DUDE_smiles.extend(np.load(dude_smiles_path+str(i)+"0000.npy",allow_pickle=True))
        GLASS_molact.extend(np.load(glass_molact_path+str(i)+"0000.npy",allow_pickle=True))
        GLASS_smiles.extend(np.load(glass_smiles_path+str(i)+"0000.npy",allow_pickle=True))
        dude_count+=10000
        glass_count+=10000
    else:
        GLASS_molact.extend(np.load(glass_molact_path+str(i)+"0000.npy",allow_pickle=True))
        GLASS_smiles.extend(np.load(glass_smiles_path+str(i)+"0000.npy",allow_pickle=True))
        glass_count+=10000
    i+=1

# activities for GLASS need to be transformed into a binary cutoff
# activities are measured in EC50 values with nM (nanomolar) units
# The cutoff is 1 mM (micromolar), which equals 1000 nM
# Molecules with EC50 below the cutoff (1000 nM) are deemed active)

# Will tell us how many molecules are exactly equal to the cutoff
print("There are "+str(GLASS_molact.count(1000))+" molecules with EC50=1000")

# The cutoff will be INCLUSIVE
GLASS_actbin=["Active" if x<=1000 else "Inactive" for x in GLASS_molact]

print("Was the cutoff applied correctly?")
print(len(GLASS_molact)==len(GLASS_actbin))

# Find and replace these samples with full lists before running on server
gl_act_samp=[GLASS_actbin[i] for i in range(1000)]
gl_smiles_samp=[GLASS_smiles[i] for i in range(1000)]
du_smiles_samp=[DUDE_smiles[i] for i in range(1000)]


### Getting Druglikeness metrics for GLASS mols and putting them into a df
qed=[] # we do want QED
lip=[] # We do want LIP
veb=[] 
gho=[]

# Note: order properties in props is:
#       MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS
for i in range(len(gl_act_samp)):
    mol=Chem.MolFromSmiles(gl_smiles_samp[i])
    qed.append(QED.default(mol))
    lip.append(Ro5(mol))
    MR=Crippen.MolMR(mol)
    n_atom=Mol.GetNumAtoms(mol)
    props=QED.properties(mol)
    veb.append(veber(props[3],props[2],props[5]))
    gho.append(ghose(props[0],props[1],MR,n_atom))

temp_dict={"SMILE":gl_smiles_samp,
            "Activity":gl_act_samp,
            "QED":qed,
            "Lipinski":lip,
            "Veber":veb,
            "Ghose":gho}
gldf=pd.DataFrame(temp_dict)

### Repeat the process now for DUDE mols
qed=[]
lip=[]
veb=[] 
gho=[]

# Note: order properties in props is:
#       MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS
for i in range(len(du_smiles_samp)):
    mol=Chem.MolFromSmiles(du_smiles_samp[i])
    qed.append(QED.default(mol))
    lip.append(Ro5(mol))
    MR=Crippen.MolMR(mol)
    n_atom=Mol.GetNumAtoms(mol)
    props=QED.properties(mol)
    veb.append(veber(props[3],props[2],props[5]))
    gho.append(ghose(props[0],props[1],MR,n_atom))

# Note there is no activity column for DUDE, becuase they are all inactive
temp_dict={"SMILE":du_smiles_samp,
            "QED":qed,
            "Lipinski":lip,
            "Veber":veb,
            "Ghose":gho}
dudf=pd.DataFrame(temp_dict)

### SAVE THE DATAFRAMES
gldf.to_pickle("dataframes/glassdf.zip")
dudf.to_pickle("dataframes/dudedf.zip")