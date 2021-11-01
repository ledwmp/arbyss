import numpy as np
from parse_cif import list_of_atoms
import sys
from atoms import chain
from calc_rmsd import align_chunks
import matplotlib.pyplot as plt

class global_alignment:
    """Alignment class to implement needleman-wunsch with gap, extend2, and extend2+

    """
    def __init__(self,gap,high_low):
        self.gap = gap
        self.high = high_low
        self.low = -1.*high_low
    def init_matrix(self):
        """Method to initialize the recursion matrices

        """
        self.rawm = np.zeros((self.len1,self.len2),float)
        for i in range(1,self.len1):
            self.rawm[i][0] = self.gap
        for j in range(1,self.len2):
            self.rawm[0][j] = self.gap
    def score_matrix(self,i,j):
        """Method to score the recursion matrices
        Args:
            i: position in peptide1
            j: position in peptide2
            amino_key: amino_acid + "_" + amino_acid
            is_trace: boolean, True if trace, False if construction of matrix
        Returns:
            Tuples of scores from diagonal matrix, insertion peptide 2, insertion peptide1
        """
        #need to replace raw matrix here with other adjust raw matrix
        if self.rawm[i][j] > self.thresh:
            score = self.low
        else:
            score = self.high
        M = ((score+self.scorem[i-1][j-1]),"d")
        X = ((score+self.scorem[i][j-1]),"l")
        Y = ((score+self.scorem[i-1][j]),"u")
        return max(M,X,Y,key = lambda x: x[0])[0],max(M,X,Y,key = lambda x: x[0])[1]

    def init_adjustmatrix(self):
        self.scorem = np.zeros((self.len1,self.len2),float)
        self.tracem = np.zeros((self.len1,self.len2),str)
        for i in range(1,self.len1):
            self.scorem[i][0] = self.low*(i)
            self.tracem[i][0] = "u"
        for j in range(1,self.len2):
            self.scorem[0][j] = self.low*(j)
            self.tracem[0][j] = "l"
        self.tracem[0][0] = "s"
        for i in range(1,self.len1):
            for j in range(1,self.len2):
                self.scorem[i][j],self.tracem[i][j] = self.score_matrix(i,j)
        #add stuff here to adjust the score matrix calculated
    def extract_scores(self):
        self.thresh = 1.5
    def score_alignment(self,chunk1,chunk2):
        """Method to help initialize matrices from peptides
        """
        self.pep1,self.pep2 = chunk1,chunk2
        walk_size = int((len(chunk1[0]) - 1)/2)
        self._onecode1,self._onecode2 = \
            [i[walk_size]._aa_one for i in chunk1],\
            [i[walk_size]._aa_one for i in chunk2]
        self.len1,self.len2 = len(self.pep1)+1,len(self.pep2)+1
        self.init_matrix()
        for i in range(1,self.len1):
            amino_1 = self.pep1[i-1]
            for j in range(1,self.len2):
                amino_2 = self.pep2[j-1]
                self.rawm[i][j] = alignment_score = align_chunks(amino_1,amino_2)
                #self.scorem[i][j],self.tracem[i][j],self.rawm[i][j] = self.score_matrix(amino_1,amino_2,i,j)
        plt.hist(self.rawm.flatten(),bins=25)
        plt.show()
        self.extract_scores()
        self.init_adjustmatrix()
        print(self.rawm)
        print(self.scorem)
        print(self.tracem)
        print(self.tracem.shape)
    def traceback(self):
        """Method to determine which matrix to trace through and create alignment
        """
        i,j = self.len1-1,self.len2-1
        self.traceback_pep1 = list()
        self.traceback_pep2 = list()
        while i > 0 or j > 0:
            pointer = self.tracem[i][j]
            if pointer == "d":
                self.traceback_pep2.append(self._onecode2[j-1])
                self.traceback_pep1.append(self._onecode1[i-1])
                i -= 1
                j -= 1
            elif pointer == "u":
                self.traceback_pep2.append("-")
                self.traceback_pep1.append(self._onecode1[i-1])
                i -= 1
            elif pointer == "l":
                self.traceback_pep2.append(self._onecode2[j-1])
                self.traceback_pep1.append("-")
                j -= 1
            #return "".join(self.traceback_pep1[::-1]),\
        #    "".join(self.traceback_pep2[::-1])
        return "".join(self.traceback_pep1[::-1]),\
                "".join(self.traceback_pep2[::-1])

my_chain_a = list_of_atoms(sys.argv[1],"A")
my_chain_b = list_of_atoms(sys.argv[1],"A")
my_chain_a = chain(my_chain_a)
my_chain_b = chain(my_chain_b)
a = my_chain_a.chunk_out(5)
b = my_chain_b.chunk_out(5)

test = global_alignment(3.8,1.)#3.8)
test.score_alignment(a[10:33],b[10:33])
print(test.traceback())
"""
for i in range(100,105):#len(my_chain_a._residuelist)):
    for j in range(75,125):#len(my_chain_b._residuelist)):
        align_chunks(my_chain_a._chunklist[i],my_chain_b._chunklist[j])
"""
