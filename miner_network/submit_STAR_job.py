#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated on June 23 12:05 PM 2021

@author: sturkars
"""
import glob
import sys
import os
import re

# data and results directories
run_dir = "/proj/omics4tb2/ishitaM/GBM_R01_Swedish/"

data_folders = glob.glob('%s/DrugRNAseq/*' %(run_dir))
#data_folders = ['%s/TL-20-C7EDD0/' %(run_dir)]
#data_folders = list(filter(lambda x:'TL-' in x, data_folders))
#data_folders = [element for element in data_folders if element in ('%s/DrugRNAseq') %(run_dir)]
print

print(data_folders)
print('Total Samples to process: %s' %(len(data_folders)))
#sys.exit()
 
jobscripts_dir = "%s/STAR_jobscripts" %(run_dir)
jobscripts_logs = "%s/STAR_job_logs" %(jobscripts_dir)


# function to create job script
def create_job_scripts(data_folder):
    
    print('Data Folder: %s' %data_folder)
    folder_name = data_folder.split('/')[-1]
    print('Folder name: %s' %folder_name)

    jobscript_file = '%s/%s_STAR.sh' %(jobscripts_dir, folder_name)
    jobscripts_logs_dir = '%s/jobscripts_logs' %jobscripts_dir
    
    # create jobscripts_dir
    if not os.path.exists(jobscripts_dir):
        os.makedirs(jobscripts_dir)

    # create jobscripts_logs_dir
    if not os.path.exists(jobscripts_logs_dir):
        os.makedirs(jobscripts_logs_dir)
    

    # Run STAR_RSEM pipeline_command
    #python run_STAR_Salmon.py $genome_dir $data_root $data_folder $out_folder
    star_cmd = 'python /proj/omics4tb2/ishitaM/GBM_R01_Swedish/code/run_STAR_RSEM.py --starPrefix %s_STAR %s %s %s/%s' %(folder_name, data_dir, folder_name,results_dir,folder_name)  
    print('     STAR CMD: ', star_cmd, '\n')           
    
    # Open jobscript for writing
    with open(jobscript_file,'w') as g:
      g.write('#!/bin/bash\n\n')
      g.write('#\n')
      g.write('#SBATCH -J %s\n'%(folder_name))
      g.write('#SBATCH -o %s/%s_outlog.txt\n' %(jobscripts_logs_dir,folder_name))
      g.write('#SBATCH -e %s/%s_outlog.txt\n' %(jobscripts_logs_dir,folder_name))

      g.write('%s\n' %star_cmd)

    g.close()

# data and results directories
data_dir = "/proj/omics4tb2/ishitaM/GBM_R01_Swedish/DrugRNAseq"
results_dir = "/proj/omics4tb2/ishitaM/GBM_R01_Swedish/process_dat"

folderCount = 1
for data_folder in data_folders:
    # get folder name
    folder_name = data_folder.split('/')[-1]
    # folder count
    folderCount = folderCount + 1
    # create jobscript
    create_job_scripts(data_folder)
    # submit jobscript
    submit_cmd = 'sbatch %s/%s_STAR.sh' %(jobscripts_dir, folder_name)
    print("Submit cmd: ", submit_cmd)
    #os.system(submit_cmd)
#    sys.exit()

