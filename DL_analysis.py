####################### Make a dataframe of QED values ####################################
from rdkit.Chem import QED
from rdkit import Chem

# we will want QED values 
# need to load in smiles and convert to mol objects again to calculate QED
smiles=np.load("smiles.npy",allow_pickle=True)

# Gerating mol objects from the smile strings

mols=[]
activity=[]
qed=[]
for i in np.arange(len(preds)): # len(preds) is equal to len(mols) FYI
    activity.append(preds[i][0])
    qed.append(QED.default(mols[i]))
    mols.append(Chem.MolFromSmiles(smiles[i]))

len(activity)==len()
len(qed)

######################## Plotting the above results ##############################
import matplotlib.pyplot as plt
#plotting params i like to use:
#the dictionary rcParams has alot of nice things in it and you can look it its keys using .keys() to see what else you can do.
plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams['font.serif'] = ['Calibri']

#predefine the figure object with size
fig = plt.figure(figsize=(3,3))

#see https://matplotlib.org/stable/gallery/color/named_colors.html for colors
plt.hist2d(activity,qed,bins=10,cmap=plt.get_cmap("gist_ncar"))
plt.scatter(activity,qed,marker='o',color='k')
#grab the axes we plotted to
axes=plt.gca()

#set some extra tick options for making it look nice
axes.tick_params(axis='both',direction='in', labelsize=10)
axes.set_xlabel("GPCR Activity Score")
axes.set_ylabel("QED")



#makes a nice layout
plt.tight_layout()

#save it
plt.savefig("activity_vs_qed.png",dpi=300,transparent=True)

plt.show()

############################### Lipinksi Ro5 ##############################
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors

# Returns True if mol passes Linpinski's rule of 5
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
        result=False
    elif strikes<=1:
        result=True
    return result

lip=[]
for i in np.arange(len(mols)): # len(preds) is equal to len(mols) FYI
    lip.append(Ro5(mols[i]))
lip=np.array(lip)

activity_arr=np.array(activity)
lip_pass=activity_arr[lip]
lip_fail=activity_arr[np.invert(lip)]

lip_data=[lip_pass,lip_fail]

# now try to make a boxplot to show these results
#predefine the figure object with size
fig = plt.figure(figsize=(3,3))

#see https://matplotlib.org/stable/gallery/color/named_colors.html for colors
plt.boxplot(lip_data)

#grab the axes we plotted to
axes=plt.gca()

#set some extra tick options for making it look nice
axes.tick_params(axis='both',direction='in', labelsize=10)
axes.set_xlabel("Lipinksi Pass/Fail")
axes.set_ylabel("GPCR Activity Score")
axes.set_xticklabels(['Pass','Fail'])



#makes a nice layout
plt.tight_layout()

#save it
plt.savefig("Lipinksi_boxplot",dpi=300,transparent=True)

plt.show()