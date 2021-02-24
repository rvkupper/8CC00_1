"""Principal Component Analysis
"""

import CellLineRMAExpression as clre
import matplotlib.pyplot as plt

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
        """return a list of RMA expression values of assigned cell lines stored 
        in self.listOfCellLineNumbers.
        """
        
        myCells = clre.CellLineRMAExpression('BRCA')
        
        rmaExpressions = []
        for cellLineName in self.listOfCellLineNames:
            rmaExpression = myCells.readRMAExpression(cellLineName)
            rmaExpressions.append(rmaExpression)
                        
        return rmaExpressions
    
    def cumulativeMovingAverage(self, x: list) -> list:
        """return a list of the cumulative moving average of input parameterlist x.
        """
        N = len(x)
        C = [0] * (N + 1)
        
        for n in range(len(x)):
            C[n + 1] = C[n] + (x[n] - C[n])/(n + 1)
            
        C.pop(0) # remove leading 0
        
        return C   
        
    def plotCumulativeMovingAverage(self, x: list, title: str = "Cumulative moving average") -> None:
        """plot the cumulative moving average of a list x
        """
        plt.plot(x)
        plt.xlabel('index')
        plt.ylabel('cumulative moving average')
        plt.title(title)
        plt.show()
        
        
    
    
l = [1, 3, 5, 11, 0, 4]    
pca = AssignmentPCA()
pca.plotCumulativeMovingAverage(l)
