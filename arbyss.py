import numpy as np
from parse_cif import list_of_atoms
import sys
from atoms import chain
from calc_rmsd import align_chunks
import matplotlib.pyplot as plt
from needleman_wunsch import needleman_wunsch
from collections import defaultdict


class global_alignment:
    """Alignment class to implement needleman-wunsch with gap, extend2, and extend2+

    """
    def __init__(self,gap):
        self.gap = gap
    def init_matrix(self):
        """Method to initialize the recursion matrices

        """
        self.rawm = np.zeros((self.len1,self.len2),float)
        self.rawm.fill(self.gap)
        self.rawm[0][0] = 0.
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
        M = ((self.rawm[i][j]+self.scorem[i-1][j-1]),"d")
        X = ((self.rawm[i][j]+self.scorem[i][j-1]),"l")
        Y = ((self.rawm[i][j]+self.scorem[i-1][j]),"u")
        return min(M,X,Y,key = lambda x: x[0])[0],min(M,X,Y,key = lambda x: x[0])[1]

    def init_adjustmatrix(self):
        self.scorem = np.zeros((self.len1,self.len2),float)
        self.tracem = np.zeros((self.len1,self.len2),str)
        for i in range(1,self.len1):
            self.scorem[i][0] = self.gap*(i)
            self.tracem[i][0] = "u"
        for j in range(1,self.len2):
            self.scorem[0][j] = self.gap*(j)
            self.tracem[0][j] = "l"
        self.tracem[0][0] = "s"
        for i in range(1,self.len1):
            for j in range(1,self.len2):
                self.scorem[i][j],self.tracem[i][j] = self.score_matrix(i,j)
        #add stuff here to adjust the score matrix calculated
    def score_alignment(self,chunk1,chunk2):
        """Method to help initialize matrices from peptides
        """
        self.tradition_alignment(chunk1,chunk2)

        def filter_both_chunks(tmp_a,tmp_b,filter_list):
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
            alignment_score = align_chunks(tmp_a,tmp_b)
            for k,l in alignment_score.items():
                tmp_dict[k].append(l)
        for i,j in tmp_dict.items():
            self.rawm[i[0]][i[1]] = min(j)
        print(self.rawm)
        self.init_adjustmatrix()
        print(self.scorem)
        print(self.tracem)
    def traceback(self):
        """Method to trace through the tracing matrix and complete alignment
        """
        i,j = self.len1-1,self.len2-1
        self.traceback_pep1 = list()
        self.traceback_pep2 = list()
        max_score = self.tracem[i][j] #which matrix to start in
        while i > 0 or j > 0:
            if max_score == "d":
                #previous matrix was M
                self.traceback_pep2.append(self.pep2[j-1]._aa_one)
                self.traceback_pep1.append(self.pep1[i-1]._aa_one)
                i -= 1
                j -= 1
            elif max_score == "l":
                #previous matrix was X
                self.traceback_pep2.append(self.pep2[j-1]._aa_one)
                self.traceback_pep1.append("-")
                j -= 1
            elif max_score == "u":
                #previous matrix was Y
                self.traceback_pep2.append("-")
                self.traceback_pep1.append(self.pep1[i-1]._aa_one)
                i -= 1
            max_score = self.tracem[i][j]
        print("".join(self.traceback_pep1[::-1]))
        print("".join(self.traceback_pep2[::-1]))
        return "".join(self.traceback_pep1[::-1]),\
            "".join(self.traceback_pep2[::-1])

    def tradition_alignment(self,a,b):
        self.alignment = needleman_wunsch("".join([i._aa_one for i in a._chunklist]),"".join([i._aa_one for i in b._chunklist]))


"""
test = global_alignment(blosum_matrix,15,7,1)
test.score_alignment("MATKGTKRSYEQMETDGERQNATEIRASVGKMIDGIGRFYIQMCTELKLSDYEGRLIQNSLTIERMVLSAFDERRNKYLEEHPSAGKDPKKTGGPIYRRVDGKWRRELILYDKEEIRRIWRQANNGDDATAGLTHMMIWHSNLNDATYQRTRALVRTGMDPRMCSLMQGSTLPRRSGAAGAAVKGVGTMVMELIRMIKRGINDRNFWRGENGRRTRIAYERMCNILKGKFQTAAQRTMVDQVRESRNPGNAEFEDLIFLARSALILRGSVAHKSCLPACVYGSAVASGYDFEREGYSLVGIDPFRLLQNSQVYSLIRPNENPAHKSQLVWMACHSAAFEDLRVSSFIRGTKVVPRGKLSTRGVQIASNENMETMESSTLELRSRYWAIRTRSGGNTNQQRASSGQISIQPTFSVQRNLPFDRPTIMAAFTGNTEGRTSDMRTEIIRLMESARPEDVSFQGRGVFELSDEKATSPIVPSFDMSNEGSYFFGDNAEEYDN","MSNMDIDSINTGTIDKTPEELTPGTSGATRPIIKPATLAPPSNKRTRNPSPERTTTSSETDIGRKIQKKQTPTEIKKSVYKMVVKLGEFYNQMMVKAGLNDDMERNLIQNAQAVERILLAATDDKKTEYQKKRNARDVKEGKEEIDHNKTGGTFYKMVRDDKTIYFSPIKITFLKEEVKTMYKTTMGSDGFSGLNHIMIGHSQMNDVCFQRSKGLKRVGLDPSLISTFAGSTLPRRSGTTGVAIKGGGTLVDEAIRFIGRAMADRGLLRDIKAKTAYEKILLNLKNKCSAPQQKALVDQVIGSRNPGIADIEDLTLLARSMVVVRPSVASKVVLPISIYAKIPQLGFNTEEYSMVGYEAMALYNMATPVSILRMGDDAKDKSQLFFMSCFGAAYEDLRVLSALTGTEFKPRSALKCKGFHVPAKEQVEGMGAALMSIKLQFWAPMTRSGGNEVSGEGGSGQISCSPVFAVERPIALSKQAVRRMLSMNVEGRDADVKGNLLKMMNDSMAKKTSGNAFIGKKMFQISDKNKVNPIEIPIKQTIPNFFFGRDTAEDYDDLDY") #A,B
#test.score_alignment("MATKGTKRSYEQMETDGERQNATEIRASVGKMIDGIGRFYIQMCTELKLSDYEGRLIQNSLTIERMVLSAFDERRNKYLEEHPSAGKDPKKTGGPIYRRVDGKWRRELILYDKEEIRRIWRQANNGDDATAGLTHMMIWHSNLNDATYQRTRALVRTGMDPRMCSLMQGSTLPRRSGAAGAAVKGVGTMVMELIRMIKRGINDRNFWRGENGRRTRIAYERMCNILKGKFQTAAQRTMVDQVRESRNPGNAEFEDLIFLARSALILRGSVAHKSCLPACVYGSAVASGYDFEREGYSLVGIDPFRLLQNSQVYSLIRPNENPAHKSQLVWMACHSAAFEDLRVSSFIRGTKVVPRGKLSTRGVQIASNENMETMESSTLELRSRYWAIRTRSGGNTNQQRASSGQISIQPTFSVQRNLPFDRPTIMAAFTGNTEGRTSDMRTEIIRLMESARPEDVSFQGRGVFELSDEKATSPIVPSFDMSNEGSYFFGDNAEEYDN","MAGQGTKRTFEQMETDSKQNTTEIRSAVGRMVKAIGRFYIQMCAELKLDDKEAVLIQNSLTIERMVLSAFDERRNKYLEEHPTVGKDPKKTGGPIYRRKEGKWEREMVLMEKENIRAIWKMANDGEENLSGLSHIMIWHSNLNDTTYQRTRALVRTGMDPRMCSLMQGSTLPRRAGAAGAAIKGVGTLIMELIRMIKRGMNDRNFWKGEQGKRTRAAYERICNNLKNKFQTAPQKAMVDQVKEGKNPGNAEIEDLLFLARSALILRGAVAHKSSLPACVYGLGVSRGFDFEREGYSLVGRDPYMLLQNSQIFSIIRKGENAAHKSQLVWMACHAAAFEDIRVSSFIKGNKIVPRGKLETRGLQIAGSETLDEALVVSLDIKSHYWAIKTRSGGNPQQSRSSAGQIAVQPTFSVQRNIPFEKKTIMAAFSNIEEGRITDMRTEIIKLMENSDPKDKVFLGRGVFEMADEKATNPIVPSLDGNDEGSYFFGDKAEEFDI") #A,BatA
#test.score_alignment("MATKGTKRSYEQMETDGERQNATEIRASVGKMIDGIGRFYIQMCTELKLSDYEGRLIQNSLTIERMVLSAFDERRNKYLEEHPSAGKDPKKTGGPIYRRVDGKWRRELILYDKEEIRRIWRQANNGDDATAGLTHMMIWHSNLNDATYQRTRALVRTGMDPRMCSLMQGSTLPRRSGAAGAAVKGVGTMVMELIRMIKRGINDRNFWRGENGRRTRIAYERMCNILKGKFQTAAQRTMVDQVRESRNPGNAEFEDLIFLARSALILRGSVAHKSCLPACVYGSAVASGYDFEREGYSLVGIDPFRLLQNSQVYSLIRPNENPAHKSQLVWMACHSAAFEDLRVSSFIRGTKVVPRGKLSTRGVQIASNENMETMESSTLELRSRYWAIRTRSGGNTNQQRASSGQISIQPTFSVQRNLPFDRPTIMAAFTGNTEGRTSDMRTEIIRLMESARPEDVSFQGRGVFELSDEKATSPIVPSFDMSNEGSYFFGDNAEEYDN","MASQGTKRSYEQMETGGERQNATEIRASVGRMVGGIGRFYIQMCTELKLSDYEGRLIQNSITIERMVLSAFDERRNKYLEEHPSAGKDPKKTGGPIYRRRDGKWVRELILYDKEEIRRIWRQANNGEDATAGLTHMMIWHSNLNDATYQRTRALVRTGMDPRMCSLMQGSTLPRRSGAAGAAVKGVGTMVMELIRMIKRGINDRNFWRGENGRRTRIAYERMCNILKGKFQTAAQRAMMDQVRESRNPGNAEIEDLIFLARSALILRGSVAHKSCLPACVYGLAVASGYDFEREGYSLVGIDPFRLLQNSQVFSLIRPNENPAHKSQLVWMACHSAAFEDLRVSSFIRGTRVAPRGQLSTRGVQIASNENMETMDSSTLELRSRYWAIRTRSGGNTNQQRASAGQISVQPTFSVQRNLPFERATIMAAFTGNTEGRTSDMRTEIIRMMESSRPEDVSFQGRGVFELSDEKATNPIVPSFDMSNEGSYFFGDNAEEYDN") #A,H5N1A

print test.traceback()
"""
my_chain_a = list_of_atoms(sys.argv[1],"A")
my_chain_b = list_of_atoms(sys.argv[2],"A")
my_chain_a = chain(my_chain_a)
my_chain_b = chain(my_chain_b)
my_chain_a.chunk_out(25,3)
my_chain_b.chunk_out(25,3)


test = global_alignment(3.8)#3.8)
test.score_alignment(my_chain_a,my_chain_b)
test.traceback()
"""
for i in range(100,105):#len(my_chain_a._residuelist)):
    for j in range(75,125):#len(my_chain_b._residuelist)):
        align_chunks(my_chain_a._chunklist[i],my_chain_b._chunklist[j])
"""
