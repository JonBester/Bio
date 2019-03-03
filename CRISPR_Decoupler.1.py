#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 21:52:30 2019

@author: jonathanbester
"""

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


################################### INPUT PARAMETERS ######################################
#replace the example below with your gRNA sequence you are using to create the CRISPR in the gene. 
#note: the cutsite should be 3bp from the left-hand (3') end of your sequence below. 
GuideRNA = "CGCACAAGTGTCCGCCCTGA" #this is an example (Guide1 for my project), delete it and replace with your 
                                    #own gRNA sequence. 

#type the file name of the wild type sequence. This should be saved as a fasta file (.fa), 
#and can be the gene or entire chromosome. For best results use just the gene, as 
#larger files will give take longer to process and may give false positives. 
wt_sequence_filename = "RiceC7BACFASTA.fa"


readme = "Paste the ab1 trace file filename here. remember to include .ab1 on the end"

#Store all sequencing data and fasta files of the wt in the same folder, then type 
#the filepath in the space above. 
filepath = "/Users/jonathanbester/20190119TestSequencingData"

################################## ADVANCED PARAMETERS ##################################

fasta_rev = False
#Are the wild type fasta sequence and RNA guide sequence in the same orientation as your sequencing trace data?
#If the wt fasta sequence is in the opposite direction to your trace data change fasta_rev below to True

guide_rev = False

#If your guide sequence is in the opposite orientation to the sequencing data, change the below to True

bracket = 200
#Bracket is the length on either side of the crispr target site that the program 
    #will check your wild type fasta for the new sequence, so it can calculate the 
    #number of deleted basepairs. Unless you have an unusually large deletion, 200bp
    #should be ample. 

freshhold = 75
#Freshold is the sensativity of the program when scanning the sequence trace file 
    #for heteroallelic sites. Trial and error has shown about 75 works well, however if 
    #you are getting unexpected results you may want to try increasing/decreasing it, 
    #errors will warn if this is the case. 

insert_length = 3 #fix the bug that gives false insertions at high insertion length estimates. 
#Insert_length is the maximum length of insert the protocol will test for. If set too high
    #this parameter can lead to false positives, and for most organisms 10bp should be ample.
    #some organisms however can have a high frequency of larger insertions due to crispr,
    #if this is the case for your system you may want to increase this variable. 

search_length = 30
#The search length is the number of bp before the crispr guide cutsite 
    #that the program will search for your indel. Typically CRISPR indels are 
    #<20bp in size, so the program is default to 50bp. You can set this higher if
    #the initial pass does not detect anything (ie 100bp), however this can cause 
    #incorrect identification of the indel if you sequencing data is of low quality.

downstream_length = 30
#The search length is the number of bp before the crispr guide cutsite 
    #that the program will search for your indel. Typically CRISPR indels are 
    #<20bp in size, so the program is default to 50bp. You can set this higher if
    #the initial pass does not detect anything (ie 100bp), however this can cause 
    #incorrect identification of the indel if you sequencing data is of low quality.

length = 14
#length is the number of base pairs you want the software to use for its prediction.
    #longer lengths are necessary for longer wt template files as they give higher prediction
    #specificity, however they also increase the probability of peak-read related errors. 
    #for general purposes a length of 14bp is likely to be more than sufficient



################################## IMPORTS ###########################################
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os
from itertools import product
import sys

#example guide RNAs:
Guide1 = "CGCACAAGTGTCCGCCCTGA"
Guide2 = "CTACAAGTGGCAGGACCTTA"

####################################### CODE ############################################

os.chdir(filepath)

#SequenceFiles, used in testing/debugging the program:
readwt2 = "448050101_C_PR_F03.ab1"
readwt = "446625701_WT_SP1-FW_G01.ab1"
read18_1 = "448385501_18-1_PRR_A08.ab1"
read18_2a = "448385501_18-2a_PRR_B08.ab1"
read18_2b = "448385501_18-2b_PRR_C08.ab1"
read18_3 = "448385501_18-3_PRR_D08.ab1"
read35B = "448385501_35B_PRR_E08.ab1"
read35C = "448385501_35C_PRR_F08.ab1"
read59C = "448385501_59C_PRR_G08.ab1"
read3_rv = "447287101_3H_PS_G11.ab1"
read3_fw = "447505101_3H_pR_E05.ab1"

read1 = SeqIO.read(readme, "abi")



GuideS = Seq(GuideRNA, IUPAC.unambiguous_dna)

if guide_rev == True:
    GuideSeq = GuideS.reverse_complement()
else:
    GuideSeq = GuideS

first = ""
second = ""

for record in SeqIO.parse(open(wt_sequence_filename), "fasta"):
    wtstring_ini = record.seq #dirty fix - would be nice if there was a standard way of reading in just a single fasta rather than using the iterator on a single value

if fasta_rev == True:
    wtstring = wtstring_ini.reverse_complement()
else:
    wtstring = wtstring_ini

CutSite = 0

for i in range(0,len(wtstring)):
    if wtstring[i:i+len(GuideSeq)] == str(GuideSeq):
        CutSite = i+3
        #print(CutSite)
        break
if CutSite == 0:
    print("Error: Defined guide sequence does not occur in provided wt sequence. Make sure the guide occurs in the wt, and check that neither needs to be reverse complemented. Please note both must match the orientation of the sequencing trace")
    sys.exit()
    
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
    print("ERROR: The read1 sequence does not correspond to the training sequence generated from the wild type fasta. Try changing the training sequence length. It is possible that the indel has occured more than {}bp upstream of the CutSite, the read quality is low, or wt sequence provided is for a different gene to that which you have sequenced.".format(search_length-15))
    sys.exit()

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

Leader_Seqs = list(scanner(poss1, wtstring_trun, length))
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
                                    Deletions.append(DeletedBp)
                                    InsertedBp = Sequ[0:k]
                                    Insertions.append(InsertedBp)
                                    NetBp = len(InsertedBp) - len(DeletedBp)
                                    Net.append(NetBp)
    
if len(Leader_Seqs) > 2:
    print("Warning. 3 or more indel candidates have been identified. One or more of these is likely to be inaccurate. Try increasing the length variable to increase test stringency")
    print("Leader sequence candidates are: {}".format(Leader_Seqs))
if len(Leader_Seqs) == 2:
    for i in range(0,2):
        print("Leader Sequence {} = {}".format(i+1, Leader_Seqs[i]))
        print("\tDeleted bases: {}".format(Deletions[i]))
        print("\tInserted bases: {}".format(Insertions[i]))
        print("\tThe position of the sequence on the origonal is: {}".format(SeqPos[i]))
        if Deletions[i] == 0 and Insertions[i] == 0:
            print("It appears that this strand is the wild type sequence, as there is no loss or gain of mutations")
        else:
            print("\tThe net loss/gain of  basepairs is: {}".format(Net[i]))
        if ((Net[i])%3) == 0:
            print("\tALERT: This mutation is non-frameshifting")
        else:
            print("\tThis is a frameshift mutation, your gene should be knocked out")
    print("\n")
    print("For your reference, the position of the crispr guide cut site is: {}".format(CutSite))
if len(Leader_Seqs) == 1:
    print("Only one indel detected. Mutation may be homoallelic. Confirm by lowering detection freshhold by reducing freshhold variable value")
    print("The leader Sequence is:", Leader_Seqs[0])
if len(Leader_Seqs) < 1:
    print("Error: No leader sequences detected. Examine raw data, possible no Cas9 event has occured.")
    

