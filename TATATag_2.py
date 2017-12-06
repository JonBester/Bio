#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 12:39:51 2017

@author: jonathanbester
"""

##   -   TATATAG: Pribnow Box Detector (E.coli)   -   ##

#   This program identifies candidate Pribnow boxes and possible sites of bacterial 
#   transcription. If you are searching for a Pribnow box in your sequence you 
#   are welcome to use this code by following the instructions below, however there 
#   are many factors that contribute to a site being a promoter - this program may 
#   suggest sites that are not promoters, or it may miss sites that are. 

##      FUNCTION:
##      Pribnow boxes are genetic motifs that enable the binding of sigma factors 
##      to sites upstream of transcription regions, which in turn recruit RNA Polymerase
##      to express the gene. They consist of a -35 element and a -10 element of
##      of six base pairs seperated by 15-19bp of interluding sequence. The architypal 
##      Pribnow box is TTCAGA (-35) and TATAAT (-10) in E. coli. However this sequence
##      never actually occurs in E. coli (you can check yourself by running the 
##      with both cut parameters set to 1.0), apparently this is because an actual 
##      TTCAGA TATAAT sequence would result in expression levels over 9000 and 
##      kill the cell. Instead the motifs typically differ from TTCAGA and TATAAT 
##      by a few base pairs. Based on the frequencies of base pair occurence in 
##      Pribnow boxes reported by Lisser and Margalit (1993), TATATAG predicts 
##      possible transcription start sites in an input sequence.
##      This program does not account for possible interaction effects of individual
##      base pairs, or UP sequences etc that also affect transcription regulation.


###       HOW TO USE:
###  1.   Either copy paste in a string on the first line behind "inputseq =" or
###       have a .txt file saved in your current directory that represents the 
###       nucleotide sequence you would like to analyse.
     
###  2.   If you use a text file fill in the file name/directory into the read_seq 
###       function where indicated, and delete "inputseq =" on the first line. 
###       Alternatively if you copy pasted a string, then delete "inputseq = 
###       read_seq(#"filename.txt here")"

###  3.   Set cut parameters. A cut parameter of 1 will only return the consensus
###       sequence, while 0 will cause your computer to explode. A cutoff of 0.1
###       is recommended, although this will miss a few promoters.

###  *    Any characters that are not T, A, C, or G will exclude the motifs in 
###       which they are present from analysis.
###  *    Big sequences take a long time, especially if you set the cut parameters
###       low. 


##Promoter sequence probability data derived from Lisser and Margalit (1993)
##does not account for UP sequences etc.
#inputseq = ##Paste your string here or write the file below##
cut35 = 0.1
cut10 = 0.1
#inputseq = 
import pandas as pd

def read_seq(inputfile):
    """takes input text file and outputs as a string. All characters that are not A,T,G or C are replaced with N, this subsequently excludes any motifs containing them from possible inclusion as a pribnow box."""
    with open(inputfile, "r") as f:
        seq = f.read()
        seq = seq.replace("\n", "")
        seq = seq.replace("r", "")
        for i in range(0,len(seq)):
            if seq[i] != "A" and seq[i] != "T" and seq[i] != "G" and seq[i] != "C":
                seq = seq.replace(seq[i], "N") 
        return seq

inputseq = read_seq(#"filename.txt" here)


##Create panda dataframes for each probability
d16n10 = {1 : [0.06,0.06,0.08,0.81,0], 2 : [0.75,0.12,0.00,0.13,0], 3 : [0.13,0.10,0.15,0.62,0], 4 : [0.69,0.12,0.12,0.08,0], 5 : [0.62,0.17,0.10,0.12,0], 6 : [0.06,0.04,0.06,0.85,0]}
df16n10 = pd.DataFrame(d16n10, index = ["A","C","G","T","N"])
0.1369612395


d17n10 = {1 : [0.06,0.15,0.08,0.71,0], 2 : [0.82,0.05,0.04,0.09,0], 3 : [0.19,0.09,0.16,0.57,0], 4 : [0.53,0.18,0.17,0.12,0], 5 : [0.55,0.18,0.09,0.18,0], 6 : [0.04,0.09,0.03,0.84,0]}
df17n10 = pd.DataFrame(d17n10, index = ["A","C","G","T","N"])
0.08125777043999999


d18n10 = {1 : [0,0.08,0.06,0.86,0], 2 : [0.71,0.04,0.08,0.18,0], 3 : [0.14,0.16,0.18,0.53,0], 4 : [0.61,0.12,0.14,0.14,0], 5 : [0.61,0.20,0.08,0.12,0], 6 : [0.12,0.02,0.02,0.84,0]}
df18n10 = pd.DataFrame(d18n10, index = ["A","C","G","T","N"])
0.10115133655199998


d16n35 = {1 : [0.19,0.08,0.04,0.69,0], 2 : [0.04,0.06,0.12,0.79,0], 3 : [0.06,0.13,0.60,0.21,0], 4 : [0.52,0.17,0.06,0.25,0], 5 : [0.19,0.56,0.13,0.12,0], 6 : [0.60,0.06,0.15,0.19,0]}
df16n35 = pd.DataFrame(d16n35, index = ["A","C","G","T","N"])
0.05714392320000001


d17n35 = {1 : [0.09,0.09,0.08,0.74,0], 2 : [0.09,0.05,0.09,0.77,0], 3 : [0.09,0.14,0.57,0.19,0], 4 : [0.56,0.20,0.11,0.13,0], 5 : [0.23,0.50,0.10,0.16,0], 6 : [0.50,0.13,0.18,0.19,0]}
df17n35 = pd.DataFrame(d17n35, index = ["A","C","G","T","N"])
0.045470039999999996


d18n35 = {1 : [0.04,0.22,0.08,0.67,0], 2 : [0.06,0.04,0.06,0.84,0], 3 : [0.10,0.06,0.73,0.12,0], 4 : [0.59,0.12,0.18,0.12,0], 5 : [0.20,0.51,0.06,0.24,0], 6 : [0.57,0.14,0.16,0.14,0]}
df18n35 = pd.DataFrame(d18n35, index = ["A","C","G","T","N"])
0.07046508697199999


#Calculate probabilities and append where appropriate
def TATATAG(seqtest):
    output = pd.DataFrame(columns = ["Trans_Start_Index", "Estimated_Start_Seq", "Distance_Between_-10and-30", "-35_Prob_Coefficient", "-35_Sequence", "-10_Prob_Coefficient", "-10_Sequence"])
    j = 0
    minus35 = df16n35
    n = 0.05714392320000001
    minus10 = df16n10
    m = 0.1369612395
    L = 16
    for i in range(35, len(seqtest)):
        prob_coefficient1 = minus35.loc[seqtest[i-35],1] * minus35.loc[seqtest[i-34],2] * minus35.loc[seqtest[i-33],3] * minus35.loc[seqtest[i-32],4] * minus35.loc[seqtest[i-31],5] * minus35.loc[seqtest[i-30],6]
        prob1 = prob_coefficient1/n
        if prob1 >= cut35:
            prob_coefficient2 = minus10.loc[seqtest[i-13],1] * minus10.loc[seqtest[i-12],2] * minus10.loc[seqtest[i-11],3] * minus10.loc[seqtest[i-10],4] * minus10.loc[seqtest[i-9],5] * minus10.loc[seqtest[i-8],6]
            prob2 = prob_coefficient2/m
            if prob2 >= cut10:
                start = seqtest[i:i+10]
                seq1 = seqtest[i-35:i-29]
                seq2 = seqtest[i-13:i-7]
                output.loc[j] = [i, start, L, prob1, seq1, prob2, seq2]
                j += 1
    minus35 = df17n35
    n = 0.045470039999999996
    minus10 = df17n10
    m = 0.08125777043999999
    L = 17
    for i in range(35,len(seqtest)):
        prob_coefficient1 = minus35.loc[seqtest[i-35],1] * minus35.loc[seqtest[i-34],2] * minus35.loc[seqtest[i-33],3] * minus35.loc[seqtest[i-32],4] * minus35.loc[seqtest[i-31],5] * minus35.loc[seqtest[i-30],6]
        prob1 = prob_coefficient1/n
        if prob1 >= cut35:
            prob_coefficient2 = minus10.loc[seqtest[i-12],1] * minus10.loc[seqtest[i-11],2] * minus10.loc[seqtest[i-10],3] * minus10.loc[seqtest[i-9],4] * minus10.loc[seqtest[i-8],5] * minus10.loc[seqtest[i-7],6]
            prob2 = prob_coefficient2/m
            if prob2 >= cut10:
                start = seqtest[i:i+10]
                seq1 = seqtest[i-35:i-29]
                seq2 = seqtest[i-12:i-6]
                output.loc[j] = [i, start, L, prob1, seq1, prob2, seq2]
                j += 1
    minus35 = df18n35
    n = 0.07046508697199999
    minus10 = df18n10
    m = 0.10115133655199998
    L = 18
    for i in range(35,len(seqtest)):
        prob_coefficient1 = minus35.loc[seqtest[i-35],1] * minus35.loc[seqtest[i-34],2] * minus35.loc[seqtest[i-33],3] * minus35.loc[seqtest[i-32],4] * minus35.loc[seqtest[i-31],5] * minus35.loc[seqtest[i-30],6]
        prob1 = prob_coefficient1/n
        if prob1 >= cut35:
            prob_coefficient2 = minus10.loc[seqtest[i-11],1] * minus10.loc[seqtest[i-10],2] * minus10.loc[seqtest[i-9],3] * minus10.loc[seqtest[i-8],4] * minus10.loc[seqtest[i-7],5] * minus10.loc[seqtest[i-6],6]
            prob2 = prob_coefficient2/m
            if prob2 >= cut10:
                start = seqtest[i:i+10]
                seq1 = seqtest[i-35:i-29]
                seq2 = seqtest[i-11:i-5]
                output.loc[j] = [i, start, L, prob1, seq1, prob2, seq2]
                j += 1
    return output

 
            
##run and print results
print(TATATAG(inputseq))

