import numpy as np

def parse_blosumfile(blosum_path):
    blosum_matrix = dict()
    with open(blosum_path) as r:
        amino_acids = next(r).strip().split()
        for line in r:
            line = line.strip().split()
            amino_acid,scores = line[0],line[1:]
            for x in range(0,len(scores)):
                subscore = scores[x]
                blosum_matrix[amino_acids[x]+"_"+amino_acid] = int(subscore)
    r.close()
    return blosum_matrix

blosum_matrix = parse_blosumfile("blosum62.txt")

class global_alignment:
    def __init__(self,submatrix,gap,extend):
        self.submatrix = submatrix
        self.gap = gap
        self.extend = extend
    def init_matrices(self):
        self.M = np.zeros((self.len1,self.len2),float)
        self.Ic = np.zeros((self.len1,self.len2),float)
        self.Ir = np.zeros((self.len1,self.len2),float)
        for i in range(1,self.len1):
            self.M[i][0],self.Ic[i][0],self.Ir[i][0] = \
                np.NINF,np.NINF,-(self.gap+self.extend*(i-1))
        for j in range(1,self.len2):
            self.M[0][j],self.Ic[0][j],self.Ir[0][j] = \
                np.NINF,-(self.gap+self.extend*(j-1)),np.NINF
        self.Ic[0][0] = np.NINF
        self.Ir[0][0] = np.NINF
    """
    def score_matrix(self,i,j,):
    """
    def score_alignment(self,peptide1,peptide2):
        self.pep1,self.pep2 = list(peptide1),list(peptide2)
        self.len1,self.len2 = len(self.pep1)+1,len(self.pep2)+1
        self.init_matrices()
        for i in range(1,self.len1):
            amino1 = self.peptide1[i-1]
            for j in range(1,self.len2):
                amino2 = self.peptide2[j-1]
                amino_key = amino1+"_"+amino2
                self.M[i][j],self.Ic[i][j],self.Ir[i][j] = \
                    self.score_matrix(i,j,amino_key)
        print self.Ic
        print self.Ir
        print self.M

test = global_alignment(blosum_matrix,5,1)
test.score_alignment("RRAA","ARTT")
