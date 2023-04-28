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

df = df.drop(['Stick_SC_1','Stick_SC_2','Stick_SC_3','Stick_SC_4','Stick_SC_5','Stick_SC_6'],axis='columns')

#%%
#clr_transformation
# kludge to prevent arithmetic errors
df[df==0] = 1e-6
geo_mean = np.exp(np.log(df.values).mean(axis=0))
df.loc[:,:] = np.log(df.values/geo_mean.reshape(1,-1))

ref_df = df.loc[:,df.columns.str.contains('Stick_SC_AGG_')]
comp_df = df.loc[:,~df.columns.str.contains('Stick')]
# %%
results_matrix = pd.DataFrame(index=['score'],columns=comp_df.columns)

for i in range(1,4):
    ref = ref_df.loc[:,ref_df.columns.str.contains('_'+str(i))]
    comp = comp_df.loc[:,comp_df.columns.str.contains('_'+str(i))]
    distances = cdist(ref.T,comp.T)
    results_matrix.loc[:,results_matrix.columns.str.contains('_'+str(i))] = distances
results_matrix = results_matrix.T
# %%
accession_names = [n.split('_')[0] for n in results_matrix.index]
accession_names = pd.unique(accession_names)
results_per_rep = pd.DataFrame(index=accession_names,columns=['score_1','score_2','score_3'])
results_per_rep.loc[:,'score_1'] = results_matrix.loc[results_matrix.index.str.contains('_1'),:].values
results_per_rep.loc[:,'score_2'] = results_matrix.loc[results_matrix.index.str.contains('_2'),:].values
results_per_rep.loc[:,'score_3'] = results_matrix.loc[results_matrix.index.str.contains('_3'),:].values

results_per_rep.to_csv('aitchison_dist.csv')
# %%
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
mds_labels = list(accession_names)
mds_labels.append('Stick')
for i,rep in enumerate(results_per_rep.columns):
    mds = MDS(dissimilarity='euclidean')
    coords = mds.fit_transform(df.loc[:,df.columns.str.contains(f'_{i+1}')].values.T)
    plt.plot(coords[:,0],coords[:,1],'.')
    #for i,label in enumerate(mds_labels):
    #    plt.annotate(label,(coords[i,0],coords[i,1]))

# %%
from sklearn.decomposition import PCA
pca_labels = list(accession_names)
df_stick_subtracted = df.copy()
df_stick_subtracted.loc[:,df_stick_subtracted.columns.str.contains('_1')] = (((df.loc[:,df.columns.str.contains('_1')]).sub((df['Stick_SC_AGG_1']),axis=0)))
df_stick_subtracted.loc[:,df_stick_subtracted.columns.str.contains('_2')] = (((df.loc[:,df.columns.str.contains('_2')]).sub((df['Stick_SC_AGG_2']),axis=0)))
df_stick_subtracted.loc[:,df_stick_subtracted.columns.str.contains('_3')] = (((df.loc[:,df.columns.str.contains('_3')]).sub((df['Stick_SC_AGG_3']),axis=0)))
df_stick_subtracted = df_stick_subtracted.loc[:,~df_stick_subtracted.columns.str.contains('Stick')]
pca = PCA(n_components=2)
coords = pca.fit_transform(df_stick_subtracted.T)
results_matrix = pd.DataFrame(data=coords,index=df_stick_subtracted.columns,columns=['dim0','dim1'])

plt.plot(results_matrix.loc[df_stick_subtracted.columns.str.contains('_1'),'dim0'],results_matrix.loc[df_stick_subtracted.columns.str.contains('_1'),'dim1'],'.')
plt.plot(results_matrix.loc[df_stick_subtracted.columns.str.contains('_2'),'dim0'],results_matrix.loc[df_stick_subtracted.columns.str.contains('_2'),'dim1'],'.')
plt.plot(results_matrix.loc[df_stick_subtracted.columns.str.contains('_3'),'dim0'],results_matrix.loc[df_stick_subtracted.columns.str.contains('_3'),'dim1'],'.')
# %%
# %%
