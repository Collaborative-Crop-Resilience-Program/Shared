#%%
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

#%% loading in data
df = pd.read_csv('/Users/maxgordon/Downloads/Hackathon_2023/STID12/ProcessedData/isolatesXPlantGeno.csv',index_col=0)
df = df.loc[df.index.str.contains('Lj'),:]
df.index = [i.split('=')[0] for i in df.index]
df = df.loc[:,~df.columns.str.contains('NI')]
df = df.loc[:,~df.columns.str.startswith('SC')]

# %%
df['Stick_SC_AGG_1'] = df[['Stick_SC_1','Stick_SC_2']].mean(axis=1)
df['Stick_SC_AGG_2'] = df[['Stick_SC_3','Stick_SC_4']].mean(axis=1)
df['Stick_SC_AGG_3'] = df[['Stick_SC_5','Stick_SC_6']].mean(axis=1)

#%%
#clr_transformation
# kludge to prevent arithmetic errors
df[df==0] = 1e-6
geo_mean = np.exp(np.log(df.values).mean(axis=0))
df.loc[:,:] = np.log(df.values/geo_mean.reshape(1,-1))

ref_df = df.loc[:,df.columns.str.contains('Stick_SC_AGG_')]
comp_df = df.loc[:,~df.columns.str.contains('Stick')]
# %%
for i in range(3):
    