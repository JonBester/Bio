################################## IMPORTS ##############################################
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os
from itertools import product
import sys

################################## INSTRUCTIONS ###########################################
'''
When a sequence read is scrabbled due to heteroallelic mutations at a CRISPR-Cas9 target site, 
this program can be used to deconvolute the read and determine the individual genotypes.
This is useful for determining if both mutations are in fact frame-shifting knockouts, or
just point mutations, and can inform which mutants we keep for breeding or experimentation. 

To run the program you need to paste your gRNA sequence (do not include PAM, 5'->3' so that 
the cut site is 3bp from the lefthand site). You must also provide the file names
of the sequencing trace file (this program only accepts ab1 format), and a fasta sequence of
the wild type gene, to be used for comparative purposes. These files should be saved together
in a folder, which you must specify the filepath to under the filepath parameter. 
The gRNA and the wt fasta should be checked to ensure they are in the same orientation 
as the sequencing trace, if they are not then reverse complimentation should be undertaken 
prior to use.

The program makes provisions for low quality sequencing data by having alterable variables
in the advanced parameters section, however for most purposes it should not be necessary to
change any of the default settings. 
'''
def crispr_decon2(gRNA, fasta, ab1, fasta_rv, guide_rv, bracket2, freshhold2, insert_length2, upstream_length2, downstream_length2, length2):
    output_strings = []
    ################################### INPUT PARAMETERS ######################################
    #replace the example below with your gRNA sequence you are using to create the CRISPR in the gene. 
    #note: the cutsite should be 3bp from the left-hand (3') end of your sequence below. 
    GuideRNA = gRNA #"CTACAAGTGGCAGGACCTTA"  #this is an example (Guide1 for my project), delete it and replace with your 
                                        #own gRNA sequence. # #"

    #type the file name of the wild type sequence. This should be saved as a fasta file (.fa), 
    #and can be the gene or entire chromosome. For best results use just the gene, as 
    #larger files will give take longer to process and may give false positives. 
    wt_sequence_filename = fasta

    print("entered answer")
    readme = ab1 #example/testing: "447505101_3H_pR_E05.ab1"
    

    #Store all sequencing data and fasta files of the wt in the same folder, then type 
    #the filepath in the space above. 
    #filepath = "/Users/jonathanbester/20190627_heterosequencingandcallichecker"
    
    ################################## ADVANCED PARAMETERS ##################################
    if fasta_rv == None:
        fasta_rev = False
    else:
        fasta_rev = True
    #Are the wild type fasta sequence and RNA guide sequence in the same orientation as your sequencing trace data?
    #If the wt fasta sequence is in the opposite direction to your trace data change fasta_rev below to True
    if guide_rv == None:
        guide_rev = False
    else:
        guide_rev = True

    #If your guide sequence is in the opposite orientation to the sequencing data, change the below to True
    bracket = int(bracket2)
    #Bracket is the length on either side of the crispr target site that the program 
        #will check your wild type fasta for the new sequence, so it can calculate the 
        #number of deleted basepairs. Unless you have an unusually large deletion, 200bp
        #should be ample. 

    freshhold = int(freshhold2)
    #Freshold is the sensativity of the program when scanning the sequence trace file 
        #for heteroallelic sites. Trial and error has shown about 75 works well, however if 
        #you are getting unexpected results you may want to try increasing/decreasing it, 
        #errors will warn if this is the case. 

    search_length = int(upstream_length2)
    #The search length is the number of bp before the crispr guide cutsite 
        #that the program will search for your indel. Typically CRISPR indels are 
        #<20bp in size, so the program is default to 50bp. You can set this higher if
        #the initial pass does not detect anything (ie 100bp), however this can cause 
        #incorrect identification of the indel if you sequencing data is of low quality.

    insert_length = int(insert_length2) #fix the bug that gives false insertions at high insertion length estimates. 
    #Insert_length is the maximum length of insert the protocol will test for. If set too high
        #this parameter can lead to false positives, and for most organisms 10bp should be ample.
        #some organisms however can have a high frequency of larger insertions due to crispr,
        #if this is the case for your system you may want to increase this variable. 

    downstream_length = int(downstream_length2)
    #The search length is the number of bp before the crispr guide cutsite 
        #that the program will search for your indel. Typically CRISPR indels are 
        #<20bp in size, so the program is default to 50bp. You can set this higher if
        #the initial pass does not detect anything (ie 100bp), however this can cause 
        #incorrect identification of the indel if you sequencing data is of low quality.

    length = int(length2)
    #length is the number of base pairs you want the software to use for its prediction.
        #longer lengths are necessary for longer wt template files as they give higher prediction
        #specificity, however they also increase the probability of peak-read related errors. 
        #for general purposes a length of 14bp is likely to be more than sufficient


    #example guide RNAs:
    Guide1 = "CGCACAAGTGTCCGCCCTGA" 
    Guide2 = "CTACAAGTGGCAGGACCTTA"

    ####################################### CODE ############################################


    read1 = SeqIO.read(readme, "abi")

    GuideS = Seq(GuideRNA, IUPAC.unambiguous_dna)

    if guide_rev == True:
        GuideSeq = GuideS.reverse_complement()
    else:
        GuideSeq = GuideS

    first = ""
    second = ""

    for record in SeqIO.parse(wt_sequence_filename, "fasta"):
        wtstring_ini = record.seq #dirty fix - would be nice if there was a standard way of reading in just a single fasta rather than using the iterator on a single value

    if fasta_rev == True:
        wtstring = wtstring_ini.reverse_complement()
    else:
        wtstring = wtstring_ini

    CutSite = 0

    for i in range(0,len(wtstring)):
        if wtstring[i:i+len(GuideSeq)] == str(GuideSeq):
            CutSite = i+3
            break
    if CutSite == 0:
        output_strings.append("Error: Defined guide sequence does not occur in provided wt sequence. Make sure the guide occurs in the wt, and check that neither needs to be reverse complemented. Please note both must match the orientation of the sequencing trace")
        return output_strings            
    
    trainer = wtstring[CutSite-search_length:CutSite-search_length+15]

            
    for i in range(0,len(read1.seq)):
        Tripper2 = False
        if read1.seq[i:i+len(trainer)] == str(trainer):
            TrainerStart = i
            FirstPeakLoc = TrainerStart
            FinalPeakLoc = TrainerStart + search_length
            Tripper2 = True
            break
    if Tripper2 == False:
        output_strings.append("ERROR: The read1 sequence does not correspond to the training sequence generated from the wild type fasta. Try decreasing the upstream length variable, or changing the training sequence length. It is possible that the indel has occured more than {}bp upstream of the CutSite, the read quality is low, or wt sequence provided is for a different gene to that which you have sequenced.".format(search_length-15))
        return output_strings

    smallprimeseq = wtstring[CutSite-search_length:CutSite+2]
    bigprimeseq = wtstring[CutSite-200:CutSite]

    base_keys = ["A", "C", "G", "T"]
    base_values = {"G": read1.annotations["abif_raw"]["DATA9"], "A": read1.annotations["abif_raw"]["DATA10"], "T": read1.annotations["abif_raw"]["DATA11"], "C": read1.annotations["abif_raw"]["DATA12"]}
    base_locations = read1.annotations["abif_raw"]["PLOC2"]

    initial = 0
    kill_switch = 0

    for i in range(FirstPeakLoc, FinalPeakLoc):
        for nucleotide in base_values:
            if base_values[nucleotide][base_locations[i]] > freshhold:
                remaining = [x for x in base_keys if x != nucleotide]
                for nucleotide2 in remaining:
                    if base_values[nucleotide2][base_locations[i]] > freshhold:
                        if kill_switch == 0:
                            initial = i
                            kill_switch = 1 #dirty fix - would be nice if could get break statement to work.
    if initial == 0:
        #there is no deviation upstream of the cutsite, therefore the is either no deviation, or the deviation must occur at the cutsite (however can be the same nucleotide on both strands)
        #to test for this, the below will examine the downstream sequence downstream_length from the cutsite for double peaking, if there is evidence of two sequences then it will set the 
        #CutSite as the point of deviation. 
        Tripped = False
        DevFreqChecker = []
        for i in range (FinalPeakLoc, FinalPeakLoc+downstream_length):
            for nucleotide in base_values:
                if base_values[nucleotide][base_locations[i]] > freshhold:
                    remaining = [x for x in base_keys if x != nucleotide]
                for nucleotide2 in remaining:
                    if base_values[nucleotide2][base_locations[i]] > freshhold:
                        Tripped = True
                        DevFreqChecker.append(i)
        if Tripped == True:
            initial = FinalPeakLoc
            if len(DevFreqChecker) < 3:
                print("WARNING: the sequence downstream of the CutSite is not highly deviated, it is possible that the results found here are due to low sequencing read quality rather than actual deviation. Please manually examine the sequencing trace to validate the results generated here")
        else:
            print("ERROR: No indels could be detected {} upstream or {} downstream of the cutsite. Try reducing the freshhold value in advanced settings, however it is possible that this sequence is wt and was not edited by the Cas9.".format(search_length, downstream_length))
            sys.exit()

    terminal = initial + length
    print("up to tie breaker")
    def tiebreaker(basenum, freshhold):
        diction = {}
        for nucleotide in base_values:
            if base_values[nucleotide][base_locations[basenum]] >= freshhold:
                diction[base_values[nucleotide][base_locations[basenum]]] = str(nucleotide)
        diction2 = dict(diction)
        first_value = max(diction.keys())
        first_value_nuc = diction[first_value]
        del diction2[first_value]
        if all(key < freshhold for key in diction2.keys()):
            second_value = first_value
        else:
            second_value = (max(diction2.keys()))
        second_value_nuc = diction[second_value]
        return first_value_nuc, second_value_nuc
            

    for i in range(initial,terminal):
        values = tiebreaker(i,freshhold)
        first += values[0]
        second += values[1]

    poss1 = [''.join(s) for s in product(*zip(first, second))]
    print(poss1)
    
    def scanner(poss, template, length):
        matches = []
        for i in range(0, len(poss)):
            if poss[i] in template:
                matches += [poss[i]]
        f_matches = set(matches)
        #print(f_matches)
            #set removes duplicates, but also destroys order. If need order use: https://stackoverflow.com/questions/7961363/removing-duplicates-in-lists
        while len(f_matches) <= 1:
            for j in range(1, insert_length):
                for i in range(0, len(poss)):
                    if poss[i][j:] in template:
                        matches += [poss[i]]
                        f_matches = set(matches)
        #print(f_matches)
        return f_matches

    wtstring_trun = wtstring[CutSite-bracket:CutSite+bracket]
    print("about to generate leader seqs")
    Leader_Seqs = list(scanner(poss1, wtstring_trun, length))
    print("leader seqs generated")
    SeqPos = []
    Deletions = []
    Insertions = []
    Net = []
    for Sequ in Leader_Seqs:
        no_ins = False
        for i in range(0,len(wtstring)):
            #first check deletions only
            if wtstring[i:i+len(Sequ)] == Sequ:
                no_ins = True
                Pos = i
                SeqPos.append(Pos)
                for m in range(0,len(wtstring)):
                    if str(read1.seq[initial-20:initial]) == str(wtstring[m-20:m]):
                        for n in range(0, 50):
                            if str(read1.seq[initial-20:initial]) + Sequ == str(wtstring[m-20:m]) + str(wtstring[m+n:m+n+len(Sequ)]):
                                DeletedBp = str(wtstring[m:m+n])
                                if DeletedBp == "":
                                    DeletedBp = 0
                                Deletions.append(DeletedBp)
                                InsertedBp = 0
                                Insertions.append(InsertedBp)
                                if DeletedBp == 0:
                                    length_DeletedBp = 0
                                else:
                                    length_DeletedBp = len(DeletedBp)
                                NetBp = ((-1)*(length_DeletedBp))
                                Net.append(NetBp)
            #now for the insertions
            if no_ins == False:
                for k in range(1, insert_length):
                    if wtstring[i:i+len(Sequ)-k] == Sequ[k:len(Sequ)]:
                        Pos = i
                        SeqPos.append(Pos)
                        for m in range(0,len(wtstring)):
                            if str(read1.seq[initial-20:initial]) == str(wtstring[m-20:m]):
                                for n in range(0, 50):
                                    if str(read1.seq[initial-20:initial]) + Sequ[k:] == str(wtstring[m-20:m]) + str(wtstring[m+n:m+n+len(Sequ[k:])]):
                                        DeletedBp = str(wtstring[m:m+n])
                                        #if DeletedBp == "":
                                            #DeletedBp = 0
                                        Deletions.append(DeletedBp)
                                        InsertedBp = Sequ[0:k]
                                        Insertions.append(InsertedBp)
                                        NetBp = len(InsertedBp) - len(DeletedBp)
                                        Net.append(NetBp)
    #print(Deletions)
    if len(Leader_Seqs) > 2:
        output_strings.append("Warning. 3 or more indel candidates have been identified. One or more of these is likely to be inaccurate. Try increasing the length variable to increase test stringency. ")
        output_strings.append("Leader sequence candidates are: {}".format(Leader_Seqs))
    if len(Leader_Seqs) == 2:
        for i in range(0,2):
            output_strings.append("Leader Sequence {} = {} ".format(i+1, Leader_Seqs[i]))
            output_strings.append("\tDeleted bases: {} ".format(Deletions[i]))
            output_strings.append("\tInserted bases: {} ".format(Insertions[i]))
            output_strings.append("\tThe position of the sequence on the origonal is: {} ".format(SeqPos[i]))
            if Deletions[i] == 0 and Insertions[i] == 0:
                output_strings.append("It appears that this strand is the wild type sequence, as there is no loss or gain of mutations ")
            else:
                output_strings.append("\tThe net loss/gain of  basepairs is: {} ".format(Net[i]))
            if ((Net[i])%3) == 0:
                output_strings.append("\tALERT: This mutation is non-frameshifting")
            else:
                output_strings.append("\tThis is a frameshift mutation, your gene should be knocked out")
            output_strings.append("\n ")
        output_strings.append("For your reference: the position of the crispr guide cut site in the fasta file is {} ".format(CutSite))
    if len(Leader_Seqs) == 1:
        output_strings.append("Only one indel detected. Mutation may be homoallelic. Confirm by lowering detection freshhold by reducing freshhold variable value ")
        output_strings.append("The leader Sequence is: ", Leader_Seqs[0])
    if len(Leader_Seqs) < 1:
        output_strings.append("Error: No leader sequences detected. Examine raw data, possible no Cas9 event has occured.")
    if output_strings == []:
        output_strings.append("something went wrong")
    print(len(output_strings))
    return output_strings

