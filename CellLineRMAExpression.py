"""Some text on this script (RMA = Rombust Multi-array Averages)
"""

import numpy as np


class CellLineRMAExpression:
    """Class for analysis of cancer cell data.
    """
    
    # instance variables
    def __init__(self, type: str):
        self.type = type
        
    def readRMAExpression(self, cellLine: str) -> list:
        """Read the RMA expression of a single cell line cellLine
        
        parameters:
        cellLine: String, cell line from which the RMA expression is to be read 
        
        returns:
        genes: list of RMA expressions of all 244 genes of cellLine.
        """
        # read data
        inf = open('data/GDSC_RNA_expression.csv')
        lines = inf.readlines()
        inf.close()
        
        nrGenes = len(lines[0])
        
        genes = None
        cellLineIndex = 0
        
        for i, line in enumerate(lines[1:]):
            lineList = line.split(',')
            if lineList[0] == cellLine:
                cellLineIndex = i 
                genes = lineList[1:]
                break
        
        # Convert strings to floats                
        for j in range(len(genes)):
            genes[j] = float(genes[j])
        
        return genes
    
    def cancerType(self, cellLine: str) -> str:
        """Return the name of the cancer type for a given cell line.
        
        parameters: 
        cellLine: String containing the name of a cell line.
        
        returns:
        cancer: String containing the name of the type of cancer with 
        which the cell line is associated.
        """
        data = []
        with open("data/GDSC_metadata.csv") as f:
          for line in f:
            line_n = line.rstrip('\n')
            data.append(line_n.split(',')[1:])
        data.pop(0)
        metadata = np.array(data)
        location = np.where(metadata.T[0] == cellLine)
        cancer = metadata.T[2][location]
        cancer = cancer[0]
        return cancer
        
