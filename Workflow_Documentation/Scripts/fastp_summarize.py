#!/usr/bin/env python

import json
import pandas as pd
import glob
import sys
import os


folder_name = sys.argv[1]

folder_path = str(os.path.abspath(folder_name))

#change directory to the folder
os.chdir(folder_path)
print('Directory:')
print(os.getcwd())


##make a list of the files in the directory
files = glob.glob("*.json")

print('File list:')
print(files)


##add files to dictionary
file_dictionary_before = {}
file_dictionary_after = {}
file_dictionary_result = {}
for file in files:
     file_json = open(file)
     json_obj = json.load(file_json)    
     json_summary = json_obj["summary"]
     #convert each piece to a data frame
     json_before = json_summary['before_filtering']
     json_before_df = pd.DataFrame.from_dict(json_before, orient = 'index')
     json_after = json_summary['after_filtering']
     json_after_df = pd.DataFrame.from_dict(json_after, orient = 'index')
     json_result = json_obj['filtering_result']
     json_result_df = pd.DataFrame.from_dict(json_result, orient = 'index')
     #add the data frames to a dictionary, with key as sample name
     file_dictionary_before[file] = json_before_df
     file_dictionary_after[file] = json_after_df
     file_dictionary_result[file] = json_result_df
    
##make data frame of all sample values
big_df_result = pd.concat(file_dictionary_result.values(), axis = 1, keys = file_dictionary_result.keys())
big_df_before = pd.concat(file_dictionary_before.values(), axis = 1, keys = file_dictionary_before.keys())
big_df_after = pd.concat(file_dictionary_after.values(), axis = 1, keys = file_dictionary_after.keys())


#calculate means
big_df_result_mean = big_df_result.copy()
big_df_result_mean['average'] = big_df_result_mean.mean(numeric_only=True, axis=1)

big_df_before_mean = big_df_before.copy()
big_df_before_mean['average'] = big_df_before_mean.mean(numeric_only=True, axis=1)

big_df_after_mean = big_df_after.copy()
big_df_after_mean['average'] = big_df_after_mean.mean(numeric_only=True, axis=1)

#to csv
big_df_result_mean.to_csv('fastp_results_summary.csv')
big_df_before_mean.to_csv('fastp_before_filtering_summary.csv')
big_df_after_mean.to_csv('fastp_after_filtering_summary.csv')


print('done')
