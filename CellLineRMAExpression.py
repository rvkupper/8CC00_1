"""Some text on this script (RMA = Rombust Multi-array Averages)
"""

class CellLineRMAExpression:
    """class for analysis of cancer data.
    """
    # class variables
    
    
    # instance variables
    def __init__(self, type: str):
        self.type = type
        
    def readRMAExpression(self, cellLine: str) -> list:
        """read the RMA expression of a single cell line cellLine
        
        parameters:
        cellLine: cell line from which the RMA expression is to be read 
        
        returns:
        genes: list of all 244 genes of cellLine.
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
        
        
someCellLine = CellLineRMAExpression('something')
genen = someCellLine.readRMAExpression('HT-29')
print(len(genen))
