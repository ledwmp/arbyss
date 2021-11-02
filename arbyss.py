import numpy as np
from parse_cif import list_of_atoms
import sys
from atoms import chain
from calc_rmsd import align_chunks
from needleman_wunsch import needleman_wunsch
from collections import defaultdict


class global_alignment:
    """Alignment class to implement structural alignment with gap

    """
    def __init__(self,gap,inner_chunk):
        self.gap = gap
        self.inner_chunk = inner_chunk
        self.walk = int((inner_chunk-1)/2)
    def init_matrix(self):
        """Method to initialize the raw scoring matrix
        """
        self.rawm = np.zeros((self.len1,self.len2),float)
        self.rawm.fill(self.gap)
        self.rawm[0][0] = 0.
    def score_matrix(self,i,j):
        """Method to score the alignment matrix
        Args:
            i: position in peptide1
            j: position in peptide2
        Returns:
            (minimum score of current box + l or d or u, l or d or u)
        """
        M = ((self.rawm[i][j]+self.scorem[i-1][j-1]),"d") #max score diagonal
        X = ((self.rawm[i][j]+self.scorem[i][j-1]),"l") #max score on left
        Y = ((self.rawm[i][j]+self.scorem[i-1][j]),"u") #max score above
        return min(M,X,Y,key = lambda x: x[0])[0],min(M,X,Y,key = lambda x: x[0])[1]

    def init_adjustmatrix(self):
        """Method to adjust the raw score matrix during alignment
        Returns:
            Filled scoring matrix and tracing matrix
        """
        #create adjusted scoring matrix and tracing matrix
        self.scorem = np.zeros((self.len1,self.len2),float)
        self.tracem = np.zeros((self.len1,self.len2),str)
        #fill adjusted and tracing matrices on edges
        for i in range(1,self.len1):
            self.scorem[i][0] = self.gap*(i)
            self.tracem[i][0] = "u"
        for j in range(1,self.len2):
            self.scorem[0][j] = self.gap*(j)
            self.tracem[0][j] = "l"
        self.tracem[0][0] = "s" #start
        for i in range(1,self.len1):
            for j in range(1,self.len2):
                self.scorem[i][j],self.tracem[i][j] = self.score_matrix(i,j)
    def extract_scores(self):
        self.thresh = 1.5
    def score_alignment(self,chunk1,chunk2):
        """Method to initialize matrices from peptides, use other methods to fill matrices
        Args:
            alignment_block object 1, alignment_block object 2
        Returns:
            Adjusted scoring matrix and tracing matrix

        """
        self.tradition_alignment(chunk1,chunk2)

        def filter_both_chunks(tmp_a,tmp_b,filter_list):
            """Filters two alignment_blocks based on residues that made it into both alignments
            """
            filt_chunk1 = [i._resraw for i in tmp_a._res_list]
            filt_chunk2 = [i._resraw for i in tmp_b._res_list]
            tmp_lista = []
            tmp_listb = []
            for i,j in filter_list:
                if i in filt_chunk1 and j in filt_chunk2:
                    tmp_lista.append(i)
                    tmp_listb.append(j)
            return ([i for i in tmp_a._res_list if i._resraw in tmp_lista],[*tmp_a._res_list]),\
                        ([i for i in tmp_b._res_list if i._resraw in tmp_listb],[*tmp_b._res_list])

        self.pep1,self.pep2 = chunk1._chunklist,chunk2._chunklist
        self.len1,self.len2 = len(self.pep1)+1,len(self.pep2)+1

        self.init_matrix()
        tmp_dict = defaultdict(list)
        for i,j in self.alignment.alignkey:
            tmp_chunk1 = chunk1._chunklist_dict[i]
            tmp_chunk2 = chunk2._chunklist_dict[j]
            tmp_a,tmp_b = filter_both_chunks(tmp_chunk1,tmp_chunk2,self.alignment.alignkey)
            alignment_score = align_chunks(tmp_a,tmp_b,self.walk)
            for k,l in alignment_score.items():
                tmp_dict[k].append(l)
        for i,j in tmp_dict.items():
            self.rawm[i[0]][i[1]] = min(j)
        self.init_adjustmatrix()

    def traceback(self):
        """Method to traverse tracing matrix in reverse and complete alignment
        """
        i,j = self.len1-1,self.len2-1
        self.traceback_pep1 = list()
        self.traceback_pep2 = list()
        max_score = self.tracem[i][j] #starting score
        while i > 0 or j > 0:
            if max_score == "d":
                #previous score was d, align both
                self.traceback_pep2.append(self.pep2[j-1]._aa_one)
                self.traceback_pep1.append(self.pep1[i-1]._aa_one)
                i -= 1
                j -= 1
            elif max_score == "l":
                #previous score was l, insert gap in pep1
                self.traceback_pep2.append(self.pep2[j-1]._aa_one)
                self.traceback_pep1.append("-")
                j -= 1
            elif max_score == "u":
                #previous score was u, insert gap in pep2
                self.traceback_pep2.append("-")
                self.traceback_pep1.append(self.pep1[i-1]._aa_one)
                i -= 1
            max_score = self.tracem[i][j]
        print("".join(self.traceback_pep1[::-1]))
        print("".join(self.traceback_pep2[::-1]))
        return "".join(self.traceback_pep1[::-1]),\
            "".join(self.traceback_pep2[::-1])
    def tradition_alignment(self,a,b):
        """Method to perform needleman-wunsch alignment based on sequence of amino acid
        Args:
            chain object 1, chain object 2
        Returns:
            global_alignment attribute that stores sequence alignment and alignment key
        """
        self.alignment = needleman_wunsch("".join([i._aa_one for i in a._chunklist]),\
                            "".join([i._aa_one for i in b._chunklist]))

def main():
    my_chain_a = list_of_atoms(sys.argv[1],"A")
    my_chain_b = list_of_atoms(sys.argv[2],"A")
    my_chain_a = chain(my_chain_a)
    my_chain_b = chain(my_chain_b)
    my_chain_a.chunk_out(25)
    my_chain_b.chunk_out(25)
    test = global_alignment(3.8,3)
    test.score_alignment(my_chain_a,my_chain_b)
    test.traceback()
if __name__ == "__main__":
    main()
