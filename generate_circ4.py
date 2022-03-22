import numpy as np
from rdkit import Chem
from rdkit import RDLogger 
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import multiprocessing as mp
import sys
from time import time

RDLogger.DisableLog('rdApp.*')
def generate_fpts(xx):
    smiles=xx[0]
    out_fpt_file=xx[1]
    out_nam_file=xx[2]
    max_nats=80
    pos_at=[1,6,7,8,9,11,15,16,17,19,35,53]
    pos_charges=[-1,0,1]
    pos_aromatic=[True,False]
    pos_bonds=["SINGLE","DOUBLE","TRIPLE","AROMATIC"]
    cicular4=[]
    cicular4_nams=[]
    for iit,smile in enumerate(smiles):
        try:
            #mol=Chem.MolFromSmiles(smile)
            #moll=Chem.AddHs(mol)
            #mol_s=Chem.MolToSmiles(moll)
            mol2=Chem.MolFromSmiles(smile)
            moli=Chem.AddHs(mol2)
            if moli is not None:
                if moli.GetNumAtoms()<=80 and all([x.GetAtomicNum() in pos_at for x in moli.GetAtoms()]) and all([x.GetFormalCharge() in pos_charges for x in moli.GetAtoms()]):
                    #arr5 = np.zeros((0,), dtype=np.float32)
                    #DataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(moli,4, nBits=1024),arr5)
                    #cicular4.append(arr5)
                    cicular4_nams.append(smiles[iit])
        except:
            pass
    #cicular4=np.array(cicular4)
    #np.save(out_fpt_file,cicular4)
    np.save(out_nam_file,cicular4_nams)

    return 1



sdf_file=sys.argv[1]
fpt_file_root=sys.argv[2]


smiles=[]
#for m in Chem.SDMolSupplier(sdf_file):
#with Chem.MultithreadedSDMolSupplier('data/5ht3ligs.sdf') as sdSupl:
with open(sdf_file, 'rb') as reader:
    sdSupl=Chem.ForwardSDMolSupplier(reader) #with Chem.SDMolSupplier(sdf_file) as sdSupl:
    for m in sdSupl:
        if m is not None:
            try:
                mol_s=Chem.MolToSmiles(m)
            except:
                mol_s=None
        if mol_s is not None:
            smiles.append(mol_s)

#fo=open("FDB-17-fragmentset.smi",'r')
#for line in fo:
#    smiles.append(line.strip())
#fo.close()

smiles=np.array(smiles)
n_per_batch=10000

n_batch=int(len(smiles)/n_per_batch)

if n_batch!=0:

    batched_smiles=np.array_split(smiles, int(len(smiles)/n_per_batch))
    
    stuff=[]
    for i in range(n_batch):
        stuff.append([batched_smiles[i],fpt_file_root+"/pc_"+str(i)+"_cicular4",fpt_file_root+"/pc_"+str(i)+"_smiles"])
    
    
    if __name__ == '__main__':
        with mp.Pool(8) as p:
            p.map(generate_fpts, stuff)
    
    
fo=open("completed.txt",'w')
fo.write("\n")
fo.close()





