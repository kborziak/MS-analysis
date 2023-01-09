from __future__ import print_function
import time
import numpy as np
import pandas as pd
#from sklearn.datasets import fetch_mldata
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
##matplotlib inline
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

#Import data
full_matrix = pd.read_csv(r"/data/liver_csc/data/merged_matrix_5_trim_norm.tsv",sep='\t',index_col=0)
sample_class = pd.read_csv(r"/data/liver_csc/data/sample_set_3_trim.tsv",sep='\t',index_col=0)

gene_set = full_matrix.columns
full_matrix['class'] = sample_class['class']
full_matrix['label'] = full_matrix['class'].apply(lambda i: str(i))

#df_subset = full_matrix[gene_set].values

#t-SNE
time_start = time.time()

tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
#tsne_results = tsne.fit_transform(df_subset)
tsne_results = tsne.fit_transform(full_matrix[gene_set].values)

print('t-SNE done! Time elapsed: {} seconds'.format(time.time()-time_start))

#df_subset['tsne-2d-one'] = tsne_results[:,0]
#df_subset['tsne-2d-two'] = tsne_results[:,1]
full_matrix['tsne-2d-one'] = tsne_results[:,0]
full_matrix['tsne-2d-two'] = tsne_results[:,1]

plt.figure(figsize=(30,30))
sns.scatterplot(
    x="tsne-2d-one", y="tsne-2d-two",
    hue="class",
#    palette=sns.color_palette("hls", 14),
    data=full_matrix,
    legend="full",
    alpha=0.3
)
plt.savefig('data/liver CGC t-SNE 5 trim norm.png')

tsne_matrix = full_matrix[['class','tsne-2d-one','tsne-2d-two']]
tsne_matrix.to_csv(r'/data/liver_csc/data/merged_matrix_5_trim_norm_tsne_matrix.tsv', sep="\t")
