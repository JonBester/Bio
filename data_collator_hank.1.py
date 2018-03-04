#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 17:34:02 2018

@author: jonathanbester
"""

# This is a python script for taking 96 well plate data from per-plate excel 
# documents and combining it into one spreadsheet where there is 96 collumns 
# representing each cell, and each row represents a different plate reading tagged
# for time and reading number. Fill in the variables below and set your working 
# directory to where you would like the combined output file to be saved, then run
# the program.



#Type the directory of the folder that your excel files are stored in between the "" below
cwd = "/Users/jonathanbester/JonB" #example
#Type a name for the output excel file between the brackets below
name = "Collated_Tecan_1"

#Define the bounds of the excel sheet using the variables below
Xstart = 1 #This is the x-axis coordinate of the start of the 96 well grid, ie if your data starts in "collumn B" then put 1 here, as it is indexed from 0. 
Xend = Xstart + 12 #This is the end of the grid, unless otherwise specified it will be Xstart+12
Ystart = 27 #This is the y-axis start of the grid, ie if your data starts in row 28 then this will be specified to 27
Yend = Ystart + 8 #This is the end of the grid, unless otherwise specified it should default to Ystart + 8

####################################### - Library Imports - #######################################
import pandas as pd
import numpy as np 
import os
import time
import xlrd
import datetime

#######################################  - Start of Code -  #########################################
now = datetime.datetime.now()
Width = []
Width.append("filename")
Width.append("Time_Created")
for a in range(1,97):
    Width.append(a)
l = os.listdir(cwd)
L = len(l)
Length = list(range(1,L+1))
df = pd.DataFrame(np.nan, index= Length, columns= Width)

count = 1
Alist = []


for filename in os.listdir(cwd):
    if filename.endswith(".xls"):
        print(count)
        workbook = xlrd.open_workbook(cwd + "/" + str(filename))
        worksheet = workbook.sheet_by_name("Sheet2")
        Alist.append(str(filename))
        T = time.ctime(os.path.getctime(str(cwd + "/" + filename)))
        Alist.append(T)
        print(Alist)
        print(worksheet.col(1)[30])
        if worksheet.col(1)[30].value != "":
            for i in range(Xstart, Xend):
                for j in range(Ystart, Yend):
                    Alist.append(worksheet.col(i)[j].value)
            df.loc[count] = Alist
            Alist = []
            print(filename)
            count += 1
        else:
            print("WARNING: " + filename + " contains empty cells, ie Tecan measurements were not made or were made atypically, all values for this reading are included as NA in the final table")
            count += 1 
            Alist = []
            continue
    elif filename.endswith(".xlsx"):
        print("WARNING: This folder contains .xlsx files - the new excel format from post-2007 versions of excel. For a quick fix set excel to save as in the old pre-2007 xls format. Alternatively, rework the above code using the Openpyxl module instead of xlrd.")
    else:
        print("WARNING: There are non-excel files in your folder, or excel files that are not of xls or xlsx format")
print(df)
name2 = name + ".xlsx"
writer = pd.ExcelWriter(name2)
df.to_excel(writer, "Sheet1")
writer.save()