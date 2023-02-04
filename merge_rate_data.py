#! /usr/bin/env python3

import os
import sys
import pandas as pd
import glob
  
#input arguments should be director/of/csvs name_protein  
# use glob to get all the csv files 
# in the folder
protein = sys.argv[1]
path = sys.argv[0]
csv_files = glob.glob(os.path.join(path, "*.csv"))
  
files = []  
# loop over the list of csv files
for f in csv_files:
    # read the csv file
    df = pd.read_csv(f)
    files.append(df) #create list of csvs
    
new_df = pd.DataFrame()

for i in len(files): #for each csv in unmerged data list
  if i == 0:
    new_df = files[i] #save first file as new df
  if i > 0:
    new_df = pd.merge(new_df, file[i], on = ["taxon","taxon order"], how = "outer") 
    #Merge first file with next file in list by taxa and taxon order.
    #Use outer join so that files are joined by union of taxa and order, because not every aa position
    #have the same group of taxa.
    #Taxa and taxa order only appears one time at the start of the mastersheet.
    
new_df.to_csv(protein+".csv")
    
    
