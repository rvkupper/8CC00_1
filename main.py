"""Main file for PCA Analysis
"""

import AssignmentPCA as apca 
import numpy as np
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

data = []
with open("data/GDSC_RNA_expression.csv") as f:
  for line in f:
    line_n = line.rstrip('\n')
    data.append(line_n.split(',')[1:])
    
data.pop(0) # remove titles

# Convert data
rawData = np.array(data)
normData = StandardScaler().fit_transform(rawData) # normalize data 

pca = apca.AssignmentPCA() # create PCA assignment instance

# Calculate covariance matrix
covData = []
for i in range(len(normData.T)):
    row = []
    for j in range(len(normData.T)):
        val = pca.covariance(normData.T[i], normData.T[j])
        row.append(val)
        
    covData.append(row)

covMatData = np.array(covData)

# Calculate eigenvectors and eigenvalues
eigVals, eigVecs = np.linalg.eig(covMatData)

# Make a list of (eigenvalue, eigenvector) tuples
eigPairs = [(np.abs(eigVals[i]), eigVecs[:,i]) for i in range(len(eigVals))]

# Sort the (eigenvalue, eigenvector) tuples from high to low
eigPairs.sort(key=lambda x: x[0], reverse=True)

# Calculate overall variance
sumEigVals = sum(eigVals)

varianceOverall = [(i / sumEigVals) * 100 for i in sorted(eigVals, reverse=True)]
           # percentage berekenen van elke PC
cumulativeOverallVariance = np.cumsum(varianceOverall) 

# Create bar graph of principal component contribution
plt.figure(figsize=(10,4))
plt.subplot(1, 2, 1)
with plt.style.context('seaborn-whitegrid'):
    # plt.figure(figsize=(6, 4))
    plt.bar(range(30), varianceOverall[:30], alpha=0.5, align='center',
            label='individual explained variance')
    plt.step(range(30), cumulativeOverallVariance[:30], where='mid',
             label='cumulative explained variance')
    plt.ylabel('Explained variance ratio')
    plt.xlabel('Principal components')
    plt.gca().set_title("Overall explained variance")
    plt.legend(loc='lower center')
    plt.tight_layout()

# scree plot 
plt.subplot(1, 2, 2)
plt.plot(eigVals[0:30], 'o-', linewidth=2, label="eigenvalues")
plt.xlabel("Principal components")
plt.ylabel("Eigenvalue")
plt.gca().set_title("PCA Scree plot")
plt.legend()    
plt.tight_layout()
    
plt.show()



    
