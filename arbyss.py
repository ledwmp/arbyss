import numpy as np
from parse_cif import list_of_atoms
import sys
from atoms import chain
from calc_rmsd import align_chunks
import matplotlib.pyplot as plt
from needleman_wunsch import needleman_wunsch

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
    def score_matrix(self,i,j,low,high):
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
        return max(M,X,Y,key = lambda x: x[0])[0],max(M,X,Y,key = lambda x: x[0])[1]

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
                self.scorem[i][j],self.tracem[i][j] = self.score_matrix(i,j,self.low,self.high)
        #add stuff here to adjust the score matrix calculated
    def score_alignment(self,chunk1,chunk2):
        """Method to help initialize matrices from peptides
        """
        self.tradition_alignment(chunk1,chunk2)
        print(self.alignment.alignkey)
        chunk1.filter_chunks([i[0] for i in self.alignment.alignkey])
        chunk2.filter_chunks([i[1] for i in self.alignment.alignkey])


        ########need to make changes here, must be able to align local alignmens of big chunks, evaluate miniature scoring on small chunks



        print([len(i[0]) for i in chunk2._chunklist])
        self.pep1,self.pep2 = chunk1._chunklist,chunk2._chunklist
        self.len1,self.len2 = len(self.pep1)+1,len(self.pep2)+1
        self.init_matrix()
        for i in range(1,self.len1):
            amino_1 = self.pep1[i-1]
            for j in range(1,self.len2):
                amino_2 = self.pep2[j-1]
                print(amino_1,amino_2)
                self.rawm[i][j] = alignment_score = align_chunks(amino_1,amino_2)
                #self.scorem[i][j],self.tracem[i][j],self.rawm[i][j] = self.score_matrix(amino_1,amino_2,i,j)
        #plt.hist(self.rawm.flatten(),bins=25)
        #plt.show()
        #self.extract_scores()
        self.init_adjustmatrix()
        print(self.scorem)
        print(self.tracem)
    def traceback(self):
        """Method to determine which matrix to trace through and create alignment
        """
        i,j = self.len1-1,self.len2-1
        self.traceback_pep1 = list()
        self.traceback_pep2 = list()
        amino_key = self.pep2[j-1]+"_"+self.pep1[i-1]
        M,X,Y = self.score_matrix(i,j,amino_key,True)
        cell_concat = M+X+Y
        max_score = cell_concat.index(max(cell_concat)) #which matrix to start in
        while i > 0 or j > 0:
            if max_score%3 == 0:
                #previous matrix was M
                self.traceback_pep2.append(self.pep2[j-1])
                self.traceback_pep1.append(self.pep1[i-1])
                i -= 1
                j -= 1
                max_score = M.index(max(M)) #matrix to continue in
            elif max_score%3 == 1:
                #previous matrix was X
                self.traceback_pep2.append(self.pep2[j-1])
                self.traceback_pep1.append("-")
                j -= 1
                max_score = X.index(max(X))
            elif max_score%3 == 2:
                #previous matrix was Y
                self.traceback_pep2.append("-")
                self.traceback_pep1.append(self.pep1[i-1])
                i -= 1
                max_score = Y.index(max(Y))
            amino_key = self.pep2[j-1]+"_"+self.pep1[i-1] #next amino acids
            M,X,Y = self.score_matrix(i,j,amino_key,True) #find scores in next cells
        return "".join(self.traceback_pep1[::-1]),\
            "".join(self.traceback_pep2[::-1])
    def tradition_alignment(self,a,b):
        self.alignment = needleman_wunsch("".join([i[2]._aa_one for i in a._chunklist_inner]),"".join([i[2]._aa_one for i in b._chunklist_inner]))


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
my_chain_a.chunk_out(25,5)
my_chain_b.chunk_out(25,5)


test = global_alignment(3.8,1.)#3.8)
test.score_alignment(my_chain_a,my_chain_b)

"""
for i in range(100,105):#len(my_chain_a._residuelist)):
    for j in range(75,125):#len(my_chain_b._residuelist)):
        align_chunks(my_chain_a._chunklist[i],my_chain_b._chunklist[j])
"""
