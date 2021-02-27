"""Main file for PCA Analysis of cancer cell RNA expression.
"""

import AssignmentPCA as apca 
import CellLineRMAExpression as clre
import numpy as np
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from math import sqrt

# Read the data
data = []
cellLines = []
with open("data/GDSC_RNA_expression.csv") as f:
  for line in f:
    line_n = line.rstrip('\n')
    data.append(line_n.split(',')[1:])
    cellLines.append(line_n.split(',')[0])

cellLines.pop(0)    
genes = data.pop(0) # remove titles


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

# Plot principal components 3D
ax = fig.add_subplot(1, 2, 2, projection="3d")
ax.tick_params(labelsize=12)
ax.set_xlabel('Principal Component 1', fontsize=16)
ax.set_ylabel('Principal Component 2', fontsize=16)
ax.set_zlabel('Principal Component 3', fontsize=16)
plt.title("3D Principal Component Analysis of Cancer Dataset",fontsize=20)
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

# Loading plots of top 50 genes 
def calcLoads(n: int, eigpairs: list, varNames: list) -> list:
    """Calculate the loads of the variables on given PC.
    
    Assumptions: 
    * PC number is in range 
    * eigpairs and varNames have the same length
    * eigpairs and varNames are not empty
    
    :param n: PC number (starting at 1).
    :param eigpairs: Sorted list (high-low) containing tuples of (eigVal, eigVec).
    :param varNames: List containing strings of the variable names in the same order as eigpairs.
    :return: List of (load, varName) tuples, sorted with highest load first.
    """
    k = n - 1
    loadings = eigpairs[k][1] * sqrt(eigpairs[k][0])
    loadingPairs = [(np.abs(np.real(loadings[i])), varNames[i]) for i in range(len(varNames))]
    loadingPairs.sort(key=lambda x: x[0], reverse=True)
    
    return loadingPairs

loadPairsPC1 = calcLoads(1, eigPairs, genes)
loadPairsPC2 = calcLoads(2, eigPairs, genes)
loadPairsPC3 = calcLoads(3, eigPairs, genes)

# Split loads and genes for plotting
loadPC1 = []
loadPC2 = []
loadPC3 = []
genesPC1 = []
genesPC2 = []
genesPC3 = []
for i in range(len(genes)):
    loadPC1.append(loadPairsPC1[i][1])
    loadPC2.append(loadPairsPC2[i][1])
    loadPC3.append(loadPairsPC3[i][1])
    genesPC1.append(loadPairsPC1[i][0])
    genesPC2.append(loadPairsPC2[i][0])
    genesPC3.append(loadPairsPC3[i][0])


# Plot that
fig = plt.figure()

ax1 = fig.add_subplot(3, 1, 1)
ax1.bar(loadPC1[:49], genesPC1[:49])
ax1.set_ylabel('Loads for PC1')
plt.xticks(rotation=90)

ax2 = fig.add_subplot(3, 1, 2)
ax2.bar(loadPC2[:49], genesPC2[:49])
ax2.set_ylabel('Loads for PC2')
plt.xticks(rotation=90)

ax3 = fig.add_subplot(3, 1, 3)
ax3.bar(loadPC3[:49], genesPC3[:49])
ax3.set_ylabel('Loads for PC3')
plt.xticks(rotation=90)

plt.tight_layout()
plt.show()
