import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

gldf=pd.read_pickle("dataframes/glassdf.zip")
dudf=pd.read_pickle("dataframes/dudedf.zip")

#plotting params i like to use:
#the dictionary rcParams has alot of nice things in it and you can look it its keys using .keys() to see what else you can do.
plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams['font.serif'] = ['Calibri']

#predefine the figure object with size
#fig = plt.figure(figsize=(6,4))

plt.close('All')
### GLASS plots
# glass ghose drug-likeness by activity
glga=gldf.groupby(['Activity','Ghose']).size().unstack()
glga.plot.bar(stacked=True,rot=0,title='GLASS Ghose results by GPCR activity')
plt.tight_layout()
plt.show()

# glass veber drug-likeness by activity
glva=gldf.groupby(['Activity','Veber']).size().unstack()
glva.plot.bar(stacked=True,rot=0,title='GLASS Veber results by GPCR activity')
plt.show()

# glass Lipinski drug-likeness by activity
glla=gldf.groupby(['Activity','Lipinski']).size().unstack()
glla.plot.bar(stacked=True,rot=0,title='GLASS Lipinksi results by GPCR activity')
plt.show()

# glass ACTIVE QED
glqactive=gldf[gldf['Activity']=="Active"]['QED']
# glass INACTIVE QED
glqinactive=gldf[gldf['Activity']=="Inactive"]['QED']
# Plotting
plt.hist(glqactive, bins=100, alpha=0.5, label="Active")
plt.hist(glqinactive, bins=100, alpha=0.5, label="Inactive")
plt.xlabel("QED")
plt.ylabel("")
plt.title("GLASS QED distributions by GPCR activity")
plt.legend(loc='upper right')
plt.show()

### DUDE plots
# 3 binary DL metrics
dudf.groupby("Ghose").size().plot.bar(rot=0,title='DUDE Ghose results')
plt.show()
dudf.groupby("Veber").size().plot.bar(rot=0,title='DUDE Veber results')
plt.show()
dudf.groupby("Lipinski").size().plot.bar(rot=0,title='DUDE Lipinksi results')
plt.show()

# QED
plt.hist(dudf['QED'],bins=100)
plt.xlabel("QED")
plt.ylabel("Count")
plt.title("DUDE QED distribution")
plt.show()

### making sure ends meet
print("Do all the DUDE plots have the same number of observations?")
print(dugo.sum()==duve.sum()==duli.sum()==dudf['QED'].count())
print("For DUDE, there should be "+str(dudf['QED'].count())+" observations")

print("Do all the GLASS plots have the same number of observations?")
print(glga.stack().sum()==glva.stack().sum()==glla.stack().sum()==gldf['QED'].count())
print("For GLASS, there should be "+str(gldf['QED'].count())+" observations")

## Actual VACC output:

# Do all the DUDE plots have the same number of observations?
# True
# For DUDE, there should be 79262 observations
# Do all the GLASS plots have the same number of observations?
# True
# For GLASS, there should be 543475 observations
