from rdkit.Chem import QED
from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors

############################### Lipinksi Ro5 ##############################
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

########################## Program Body ##################################
# need to load in GPCR predictions to harvest activity data
smiles=np.load("numpy_objs/smiles.npy",allow_pickle=True)
preds=np.load("numpy_objs/GPCR_predictions.npy",allow_pickle=True)

# This for loop will be generating 4 seperate lists for each chemical scored in preds: mols, QED scores, Lipinski, and the 1st activity score (not 1-activity score)
mols=[]
activity=[]
qed=[]
lip=[]
for i in np.arange(len(preds)): # len(preds) is equal to len(mols) FYI
    activity.append(preds[i][0])
    mols.append(Chem.MolFromSmiles(smiles[i]))
    qed.append(QED.default(mols[i]))
    lip.append(Ro5(mols[i]))

# should return true
len(activity)==len(mols)==len(qed)==len(preds)==len(lip)

######################## Plotting QED and Activity score - histogram ##############################
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
plt.savefig("plots/activity_vs_qed.png",dpi=300,transparent=True)

plt.show()


########################## Plotting Activity score by Ro5 pass/fail - Boxplot ######################
lip=np.array(lip) # convert to array for easy boolean indexing
activity=np.array(activity) # dido

# subset of activity scores for Ro5 passes
lip_pass=activity[lip]

# subset of activity scores for Ro5 fails
lip_fail=activity[np.invert(lip)]

# combining the data for the box plot
lip_data=[lip_pass,lip_fail]

# now try to make a boxplot to show these results
plt.boxplot(lip_data)

#grab the axes we plotted to
axes=plt.gca()

#set some extra tick options for making it look nice
axes.set_xlabel("Lipinksi Pass/Fail")
axes.set_ylabel("GPCR Activity Score")
axes.set_xticklabels(['Pass','Fail'])

#save it
plt.savefig("Lipinksi_boxplot",dpi=300,transparent=True)

plt.show()

############################# Other Plots to make and analysis to do ###########################
# 1. Graph distribution of QED by activity>0.5 and activity<0.5 
#   - this could be two histograms, boxplots, violin plots, whatever else looks nice
# 2. Repeat #1 but with other drug-likness metrics, such as Lipinski, Veber, Ghose
# 3. Perform the Kolmogorov–Smirnov test (probably on scipy) to see if the two distributions are statistically different