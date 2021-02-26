"""Main file for PCA Analysis
"""

import AssignmentPCA as apca 
import CellLineRMAExpression as clre
import numpy as np
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

data = []
cellLines = []
with open("data/GDSC_RNA_expression.csv") as f:
  for line in f:
    line_n = line.rstrip('\n')
    data.append(line_n.split(',')[1:])
    cellLines.append(line_n.split(',')[0])

cellLines.pop(0)    
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

# # Create bar graph of principal component contribution
# plt.figure(figsize=(10,4))
# plt.subplot(1, 2, 1)
# with plt.style.context('seaborn-whitegrid'):
#     # plt.figure(figsize=(6, 4))
#     plt.bar(range(30), varianceOverall[:30], alpha=0.5, align='center',
#             label='individual explained variance')
#     plt.step(range(30), cumulativeOverallVariance[:30], where='mid',
#              label='cumulative explained variance')
#     plt.ylabel('Explained variance ratio')
#     plt.xlabel('Principal components')
#     plt.gca().set_title("Overall explained variance")
#     plt.legend(loc='lower center')
#     plt.tight_layout()
# 
# # scree plot 
# plt.subplot(1, 2, 2)
# plt.plot(eigVals[0:30], 'o-', linewidth=2, label="eigenvalues")
# plt.xlabel("Principal components")
# plt.ylabel("Eigenvalue")
# plt.gca().set_title("PCA Scree plot")
# plt.legend()    
# plt.tight_layout()
# 
# plt.show()

# Reducing dimensions
matrixW = np.hstack((eigPairs[0][1].reshape(len(eigPairs),1),
                      eigPairs[1][1].reshape(len(eigPairs),1),
                      eigPairs[2][1].reshape(len(eigPairs),1)))

newSubspace = normData.dot(matrixW)

# Create CellLineRMAExpression instance for cancertype determination
myCells = clre.CellLineRMAExpression('Cancer')

# Plot principal components 2D
fig = plt.figure(figsize=plt.figaspect(2.))
plt.subplot(1, 2, 1)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Principal Component 1', fontsize=16)
plt.ylabel('Principal Component 2', fontsize=16)
plt.title("2D Principal Component Analysis of Cancer Dataset",fontsize=20)
targets = ['BRCA', 'COAD/READ', 'KIRC', 'NB']
colors = {'BRCA':'r', 'COAD/READ':'g', 'KIRC':'b', 'NB':'k'}
for i in range(len(newSubspace)):
    cancer = myCells.cancerType(cellLines[i])    
    color = colors[cancer]
    x = np.real(newSubspace[i][0])
    y = np.real(newSubspace[i][1])
    plt.scatter(x, y, c = color, s = 50)

# Create legend
markers = [plt.Line2D([0,0],[0,0],color=color, marker='o', linestyle='') for color in colors.values()]
plt.legend(markers, colors.keys(), numpoints=1)

ax = fig.add_subplot(1, 2, 2, projection="3d")
ax.tick_params(labelsize=12)
# plt.yticks(fontsize=12)
# plt.zticks(fontsize=12)
ax.set_xlabel('Principal Component 1', fontsize=16)
ax.set_ylabel('Principal Component 2', fontsize=16)
ax.set_zlabel('Principal Component 3', fontsize=16)
plt.title("3D Principal Component Analysis of Cancer Dataset",fontsize=20)
targets = ['BRCA', 'COAD/READ', 'KIRC', 'NB']
colors = {'BRCA':'r', 'COAD/READ':'g', 'KIRC':'b', 'NB':'k'}
for i in range(len(newSubspace)):
    cancer = myCells.cancerType(cellLines[i])    
    color = colors[cancer]
    x = np.real(newSubspace[i][0])
    y = np.real(newSubspace[i][1])
    z = np.real(newSubspace[i][2])
    ax.scatter(x, y, z, c = color)

# Create legend
markers = [plt.Line2D([0,0],[0,0],color=color, marker='o', linestyle='') for color in colors.values()]
plt.legend(markers, colors.keys(), numpoints=1)
plt.tight_layout()
plt.show()
