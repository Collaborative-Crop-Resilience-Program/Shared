#%% Load Modules
import pandas as pd
import numpy as np
from pgmpy.estimators import MmhcEstimator
from pgmpy.estimators import HillClimbSearch, BicScore
from skbio.stats.composition import clr

#%% Load File
f1 = pd.read_csv('/Users/bkamos/Desktop/DesktopFolders/InRootHackathon/STID12/isolatesXPlantGeno_trunc_noSC.csv', index_col = False)
isolates = f1['Isolate']
f1.set_index('Isolate')

#%%
# batch1Cols = f1.filter(regex='_1')
batch2Cols = f1.filter(regex='_2')
# batch3Cols = f1.filter(regex='_3')

# batch1Cols = batch1Cols.assign(Isolates = isolates)
batch2Cols = batch2Cols.assign(Isolates = isolates)
# batch3Cols = batch3Cols.assign(Isolates = isolates)


# batch1Cols = batch1Cols.set_index('Isolates')
batch2Cols = batch2Cols.set_index('Isolates')
# batch3Cols = batch3Cols.set_index('Isolates')

#%%
row_vars = batch2Cols.var(axis=1)
row_med = row_vars.describe()
rows_to_drop = batch2Cols[row_vars<row_med[4]].index
batch2Cols.drop(rows_to_drop, axis=0, inplace=True)
batch2Cols = batch2Cols.T
# batch2Cols = batch2Cols.iloc[1:,:]
batch2Cols += 1e-6
batch2Cols = batch2Cols.astype(float)

#%%
dataIndex = batch2Cols.index
dataColummns = batch2Cols.columns
data = clr(batch2Cols)
data = pd.DataFrame(data, columns = dataColummns, index= dataIndex)
est = HillClimbSearch(data)
best_model = est.estimate()
# %%
node1 = []
node2 = []
for i in best_model.edges:
    node1.append(i[0])
    node2.append(i[1])

# %%
Network = pd.DataFrame(
    {'Source': node1,
     'Target': node2,
    })
# %%
Network.to_csv('/Users/bkamos/Documents/GitHub/Shared/Programs/Hackathon/Teams/Team1/OutputNetworks/Batch2Network_75x.csv', index = False)
# %%# %%
