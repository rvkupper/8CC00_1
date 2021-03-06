"""Principal Component Analysis tools
"""

import CellLineRMAExpression as clre
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import numpy as np
from math import sqrt

class AssignmentPCA:
    """Class for principal component analysis of the CellLineRMAExpression data.
    """
    listOfCellLineNumbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 
    14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 31, 32, 
    33, 34, 35, 36, 37, 38, 39, 40, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 
    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 
    70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 
    88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 
    105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 
    119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 
    133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 
    147] # copied from email
    
    def __init__(self):
        inf = open('data/GDSC_RNA_expression.csv')
        lines = inf.readlines()
        inf.close()
        
        listOfCellLineNames = []
        for i in self.listOfCellLineNumbers:
            line = lines[i].split(',')
            cellLineName = line[0]
            listOfCellLineNames.append(cellLineName)
            
        self.listOfCellLineNames = listOfCellLineNames
    
    def readRMAExpressionAssigned(self) -> list:
        """Return a list of RMA expression values of assigned cell lines stored 
        in self.listOfCellLineNumbers.
        """
        
        myCells = clre.CellLineRMAExpression('BRCA')
        
        rmaExpressions = []
        for cellLineName in self.listOfCellLineNames:
            rmaExpression = myCells.readRMAExpression(cellLineName)
            rmaExpressions.append(rmaExpression)
                        
        return rmaExpressions
    
    def cumulativeMovingAverage(self, x: list) -> list:
        """Return a list of the cumulative moving average of input parameterlist x.
        
        Assumption: x contains numbers or is empty.
        
        :param x: list of numbers.
        :return: list of float cumulative moving averages.
        
        >>> self.cumulativeMovingAverage([1, 3, 5, 11, 0, 4])
        [1.0, 2.0, 3.0, 5.0, 4.0, 4.0]
        >>> self.cumulativeMovingAverage([])
        []
        """
        N = len(x)
        C = [0] * (N + 1)
        
        for n in range(len(x)):
            C[n + 1] = C[n] + (x[n] - C[n])/(n + 1)
            
        C.pop(0) # remove leading 0
        
        return C   
        
    def plotCumulativeMovingAverage(self, x: list, title: str = "Cumulative moving average") -> None:
        """Plot the cumulative moving average of a list x
        
        Assumption: x contains numbers or is empty.
            
        :param x: List of parameters to be calculated moving average of and plotted.
        :param title: string containing a title for the graph.
        """
        C = self.cumulativeMovingAverage(x)
        plt.plot(C)
        plt.xlabel('index')
        plt.ylabel('cumulative moving average')
        plt.title(title)
        plt.show()
        
    def covariance(self, param1: list, param2: list) -> float:
        """Return the covariance of parameter lists param1 and param2.        
        
        Assumption: param1 and param2 contain numbers and are of equal length.
        
        :param param1: List of parameters to be compared.
        :param param2: List of parameters to compare with .
        :return: covariance of param1 and param2.
        
        >>> self.covariance([1, 3, 5, 11, 0, 4], [2, 6, 2, 78, 1, 4])
        106.4
        >>> self.covariance([1], [1, 2])
        Traceback (most recent call last):
            ...
        AssertionError: Parameter lists must be of the same length.    
        """
        assert len(param1) == len(param2), "Parameter lists must be of the same length."
        
        n = len(param1)
        
        mean1 = np.mean(param1)
        mean2 = np.mean(param2)
                
        arr1 = np.array(param1)
        arr2 = np.array(param2)
                
        arr1_diff = arr1 - mean1
        arr2_diff = arr2 - mean2
                
        multiplied = arr1_diff * arr2_diff
        sumMultiplied = sum(multiplied)
        covar = sumMultiplied/(n - 1.0)
        
        return covar
        
    def calcLoads(self, n: int, eigpairs: list, varNames: list) -> list:
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
