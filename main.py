"""Main file for PCA Analysis
"""

import AssignmentPCA as apca 
import numpy as np
from sklearn.preprocessing import StandardScaler

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




    
