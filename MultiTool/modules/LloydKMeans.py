#This script uses Lloyd's Kmeans algorithm to cluster data points and plot clustered data sets. It also calculates the overall
#accuracy of the run.

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.datasets import load_iris
import sklearn.metrics as sm
import numpy as np
import pylab as pl
import pandas as pd

#Because not many people know formal data structures of pandas or pca, this example works by modifying the first two columns
#of 'iris.csv'. This dataset may be found in the sklearn directory under your current python installation directory. If
#uncertain of the location, run python -v and import sklearn.
iris = load_iris()

x = pd.DataFrame(iris.data)
x.columns = ['TsTv','pKaKs']
 
y = pd.DataFrame(iris.target)
y.columns = ['Targets']

kmeans = KMeans(n_clusters=4, random_state=111)
kmeans.fit(iris.data)
pl.figure('K-means with 4 clusters')
pl.scatter(pca_2d[:, 0], pca_2d[:, 1], c=kmeans.labels_)

pl.subplot(1, 2, 1)
for i in range(0, pca_2d.shape[0]):
	if iris.target[i] == 1:
		c1 = pl.scatter(x.TsTv[i],x.pKaKs[i],c='r',marker='o',s=80)
	elif iris.target[i] == 0:
		c2 = pl.scatter(x.TsTv[i],x.pKaKs[i],c='g',marker='o',s=80)
	elif iris.target[i] == 2:
		c3 = pl.scatter(x.TsTv[i],x.pKaKs[i],c='b',marker='o',s=80)
	elif iris.target[i] == 3:
		c4 = pl.scatter(x.TsTv[i],x.pKaKs[i],c='y',marker='o',s=80)


#Title and axis labels of previous plot
pl.title('Real Species Classification')
pl.ylabel('pKaKs')
pl.xlabel('Ts/Tv')

pl.subplot(1, 2, 2)
for i in range(0, pca_2d.shape[0]):
	if kmeans.labels_[i] == 1:
		c1 = pl.scatter(x.TsTv[i],x.pKaKs[i],c='r',marker='o',s=80)
	elif kmeans.labels_[i] == 0:
		c2 = pl.scatter(x.TsTv[i],x.pKaKs[i],c='g',marker='o',s=80)
	elif kmeans.labels_[i] == 2:
		c3 = pl.scatter(x.TsTv[i],x.pKaKs[i],c='b',marker='o',s=80)
	elif kmeans.labels_[i] == 3:
		c4 = pl.scatter(x.TsTv[i],x.pKaKs[i],c='y',marker='o',s=80)

#If needed, add more organisms (remember to add the scatter plotting for them within the for loops as well!)
pl.legend([c1, c2, c3,c4],['Apis mellifera', 'Arabidopsis thaliana','Caenorhabditis briggsae','C. elegans (N2)'], bbox_to_anchor=(1.025, 1), loc=2, borderaxespad=0.)

#Title and axis labels of previous plot
pl.title('K-Means Species Classification | K=4')
pl.xlabel('Ts/Tv')
pl.ylabel('pKaKs')
pl.show()

#Calculate accuracy rate of algorithm
kmeans2 = KMeans(n_clusters=4, random_state=111)
kmeans2.fit(x)
#Convert 1s to 0s and 0s to 1s.
predY = np.choose(kmeans2.labels_, [1, 0, 2]).astype(np.int64)
print(sm.accuracy_score(y, predY)) 

