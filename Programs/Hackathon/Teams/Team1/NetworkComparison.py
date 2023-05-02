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
# %% drop rows with low variance (microbial samples) and then transpose the data
row_vars = f1.var(axis=1)
row_med = row_vars.describe()
rows_to_drop = f1[row_vars<row_med[6]].index
f1.drop(rows_to_drop, axis=0, inplace=True)
f1transposed = f1.T

#%%
# row_vars = batch1Cols.var(axis=1)
# row_med = row_vars.describe()
# rows_to_drop = f1[row_vars<row_med[5]].index
# batch1Cols.drop(rows_to_drop, axis=0, inplace=True)
# batch1Cols = batch1Cols.T
f1transposed = f1transposed.iloc[1:,:]
f1transposed += 1e-6
f1transposed = f1transposed.astype(float)

# f1transposed1 = (f1transposed1-f1transposed1.min())/(f1transposed1.max()-f1transposed1.min())
#%%
dataIndex = f1transposed.index
dataColummns = f1transposed.columns
data = clr(f1transposed)
data = pd.DataFrame(data, columns = dataColummns, index= dataIndex)
est = HillClimbSearch(data)
best_model = est.estimate()

# est = MmhcEstimator(f1transposed1)
# model = est.estimate()
# print(model.edges())
# %%
# data = pd.DataFrame(np.random.randint(0, 2, size=(2500, 4)), columns=list('XYZW'))
# data['sum'] = data.sum(axis=1)
# est = MmhcEstimator(data)
# model = est.estimate()
# print(model.edges())
# %%
