#%%
import pandas as pd
import numpy as np
from pgmpy.estimators import MmhcEstimator
from pgmpy.estimators import HillClimbSearch, BicScore

#%%
f1 = pd.read_csv('/Users/bkamos/Desktop/DesktopFolders/InRootHackathon/STID12/isolatesXPlantGeno_trunc.csv', index_col = False)
# %%
row_vars = f1.var(axis=1)
row_med = row_vars.describe()
rows_to_drop = f1[row_vars<=row_med[6]].index
f1.drop(rows_to_drop, axis=0, inplace=True)
f1transposed = f1.T
# %%

# colNames = []
# split_cols = f1.columns
# split_cols.to_list
# for i in split_cols:
#     split = i.split('_')
#     if len(split) == 1:
#         colNames.append(split)
#     if len(split) == 2:
#         colNames.append(split[0])
#     else:
#         colNames.append(split[:2])

# colNames.pop(0)
# %%
# f1.columns = colNames
# %%
f1transposed1 = f1transposed.iloc[1:,:]

est = HillClimbSearch(f1transposed1)
best_model = est.estimate(scoring_method=BicScore(f1transposed1))

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
