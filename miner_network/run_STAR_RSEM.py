#!/usr/bin/env python3

#############################################################
##### RNASeq Analysis Pipeline with STAR, RSEM & Salmon #####
##### Last update: 01/11/2023 Serdar Turkarslan         #####
##### Institute for Systems Biology                     #####
############################################################
import glob, sys, os, string, datetime, re
import argparse

DESCRIPTION = """run_STAR_SALMON.py - run STAR and Salmon"""

####################### Create results directories ######################################
def create_dirs(data_trimmed_dir, fastqc_dir, results_dir, rsem_results_dir):
    dirs = [data_trimmed_dir, fastqc_dir, results_dir, rsem_results_dir]
    print()
    for dir in dirs:
        # create results folder
        #print(dir)
        if not os.path.exists('%s' %(dir)):
            print('\033[31m %s directory DOES NOT exist. Creating. \033[0m' %(dir))
            os.makedirs('%s' %(dir))
        else:
            print('\033[31m %s directory exists. Not creating. \033[0m' %(dir))
    print()

####################### Trimgalore for quality and trimming ###############################
def trim_galore(first_pair_file, second_pair_file, folder_name, sample_id, file_ext, data_trimmed_dir,
                fastqc_dir):
    #print("1stpair:%s, 2ndpair:%s, folder_name:%s, sample_name:%s")%(first_pair_file,second_pair_file,folder_name,sample_name)
    print()
    print ("\033[33m*** Running TrimGalore \033[0m")
    
    # create sample spepcific trimmed directory
    if not os.path.exists('%s' %(data_trimmed_dir)):
        os.makedirs('%s' %(data_trimmed_dir))  
    # create sample spepcific fastqcdirectory
    if not os.path.exists('%s' %(fastqc_dir)):
        os.makedirs('%s' %(fastqc_dir))
   
   # TrimGalore Run Command
    cmd = 'trim_galore --fastqc_args "--outdir %s/" --paired --output_dir %s/ %s %s' %(fastqc_dir,data_trimmed_dir,first_pair_file, second_pair_file)
    print( '    ++++++ Trimgalore Command:', cmd)
    os.system(cmd)


####################### Collect trimmed data files ###############################
def collect_trimmed_data(data_trimmed_dir, file_ext):
    # define result files
    if file_ext == "gz":
        first_pair_trimmed = glob.glob('%s/*_val_1.fq.gz'%(data_trimmed_dir))
        second_pair_trimmed = glob.glob('%s/*_val_2.fq.gz'%(data_trimmed_dir))
    else:
        first_pair_trimmed = glob.glob('%s/*_val_1.fq'%(data_trimmed_dir))
        second_pair_trimmed = glob.glob('%s/*_val_2.fq'%(data_trimmed_dir))
    print()
    print('     Trimmed Files:\n 1st:%s \n 2nd:%s' %(first_pair_trimmed,second_pair_trimmed))
    
    first_pair_group = ','.join(first_pair_trimmed)
    second_pair_group = ','.join(second_pair_trimmed)
    pair_files = []
    for file in first_pair_trimmed:
        mate_file = file.replace('_1_val_1.fq','2_val_2.fq')
        paired_mates = file + ' ' + mate_file
        pair_files.append(paired_mates)

    star_input_files = ' '.join(pair_files)

    return first_pair_group,second_pair_group, star_input_files


####################### Run STAR #####################################
def run_star(first_pair_group, second_pair_group, results_dir, star_input_files,
             folder_name, genome_dir):
    print()
    print('\033[33m*** Running STAR! \033[0m')

    outfile_prefix = '%s/%s' %(results_dir, args.starPrefix)

    # ## Create tmp directory
    # outfile_tmp_prefix = outfile_prefix + '_STARtmp'
    # print("Temo file" + outfile_tmp_prefix)
    # if not os.path.exists('%s' %(outfile_tmp_prefix)):
    #     print('\033[31m %s directory DOES NOT exist. Creating. \033[0m' %(dir))
    #     os.makedirs('%s' %(outfile_tmp_prefix))
    # else:
    #     print('\033[31m %s directory exists. Not creating. \033[0m' %(dir))

    star_options = """--runThreadN 32 \
    --outSAMattributes All \
    --genomeLoad LoadAndKeep \
    --outFilterType Normal \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonical \
    --outSAMtype BAM Unsorted \
    --limitBAMsortRAM 5784458574 \
    --readFilesCommand zcat \
    --outReadsUnmapped Fastx \
    --outFilterMatchNminOverLread 0.33 \
    --outFilterScoreMinOverLread 0.33 \
    --outFilterMultimapScoreRange 1 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 10 \
    --alignIntronMax 500000 \
    --alignMatesGapMax 1000000 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 2 \
    --quantMode TranscriptomeSAM GeneCounts \
    """ 

    cmd = 'STAR --genomeDir %s/STAR_reference %s --readFilesIn %s %s --outFileNamePrefix %s' % (genome_dir, star_options,first_pair_group, second_pair_group, outfile_prefix)
    print(' STAR run command:%s' %cmd)
    os.system(cmd)

####################### Run Salmon Count ###############################
def run_salmon_quant(results_dir, folder_name, genome_fasta):
    star_outfile_prefix = '%s/%s_' %(results_dir, args.starPrefix)
    print()
    print('\033[33m*** Running salmon-quant! \033[0m')
    salmon_input = '%sAligned.out.bam' % (star_outfile_prefix)

    cmd = 'salmon quant -t %s -l A -a %s -o %s/salmon_quant_%s' % (genome_fasta, salmon_input, results_dir,args.salmonPrefix)
    print('salmon-count run command:%s' %cmd)
    #os.system(cmd)


####################### Run HTSEq Count ###############################
# def run_htseq(htseq_dir, results_dir, folder_name, genome_gff):
#     print()
#     print('\033[33m*** Running htseq-count! \033[0m')
#     htseq_input = '%s/%s_Aligned.sortedByCoord.out.bam' %(results_dir, folder_name)
#     cmd = 'htseq-count -s "reverse" -t "exon" -i "Parent" -r pos --max-reads-in-buffer 60000000 -f bam %s %s > %s/%s_htseqcounts.txt' %(htseq_input,genome_gff,htseq_dir,folder_name)
#     print('htseq-count run command:%s' %cmd)
#     #os.system(cmd)

####################### Run RSEM Calculate Expression ###############################
# We are using prebuilt index for human genome
def run_RSEM(results_dir,rsem_results_dir,folder_name):
    # build RSEM reference index
    print()
    print ("\033[33m*** Building RSEM index... \033[0m")
    rsem_ref_cmd = 'rsem-prepare-reference --gtf %s %s %s/RSEM_reference/human_ref' %(genome_gff, genome_fasta, genome_dir)
    print(rsem_ref_cmd)
    #os.system(rsem_ref_cmd)

    print()
    print ("\033[33m*** Running RSEM... \033[0m")
    rsem_cmd = 'rsem-calculate-expression -p 32 --paired-end --estimate-rspd --append-names --bam %s/%s_STARAligned.toTranscriptome.out.bam %s/RSEM_reference/human_ref %s/%s' %(results_dir, folder_name, genome_dir,rsem_results_dir,folder_name)
    print(rsem_cmd)
    os.system(rsem_cmd)


####################### Create STAR index ###############################
def create_genome_index(genome_dir, genome_fasta):
    index_cmd = 'STAR --runMode genomeGenerate --runThreadN 32 --genomeDir %s/STAR_reference --genomeFastaFiles %s --genomeSAindexNbases 12   --genomeSAsparseD 3   --limitGenomeGenerateRAM 15000000000  --limitSjdbInsertNsj 383200   --sjdbGTFfile %s  --sjdbOverhang 149' %(genome_dir, genome_fasta, genome_gff)
    print()
    print ("\033[33m*** Indexing genome for STAR... \033[0m")
    print(index_cmd)
    print()
    if os.path.exists('%s/STAR_reference/SAindex' % (genome_dir)):
        print ('*** Genome indexes exist. Not creating!')
    else:
        print ('*** Creating genome indexes')
        #os.system(index_cmd)


####################### Running the Pipeline ###############################

def run_pipeline(data_folder, results_folder, genome_dir, genome_fasta, genome_gff):
    folder_count = 1
    error_file = "run_errors.txt"

    # Loop through each data folder
    #for data_folder in data_folders:
    folder_name = data_folder.split('/')[-1]
    print()
    print('\033[34m Processing Folder: %s\033[0m' % (folder_name))

    # Get the list of first file names in paired end sequences
    DATA_SEARCH1 = '%s/*R1*.fastq*' % data_folder # test
    #DATA_SEARCH1 = '%s/*_1.fq*' % data_folder
    print("SEARCHING FIRST PAIRS IN: ", DATA_SEARCH1)
    first_pair_files1 = glob.glob("%s/*[_.R]1*.fastq.gz" %(data_folder), recursive=True) # Search first set of folders
    first_pair_files2 = glob.glob("%s/*/*[_.R]1*.fastq.gz" %(data_folder), recursive=True) # search additional folders if there are different levels
    # combine multiple searches
    first_pair_files = first_pair_files1 + first_pair_files2
    
    if len(first_pair_files) == 0:
        print("Error: I didnt find any file with '*fastq*' extension in %s\n" %(folder_name))
        with open(error_file,'w') as f:
            f.write("Error: I didnt find any file with '*fastq*' extension in %s" %(folder_name))
            f.write(str(first_pair_files))
        f.close()    
    
    #concat
    #first_pair_files = glob.glob('%s/*_1.fq*' % (data_folder))
    #second_pair_files = glob.glob('%s/_R2*.fastq*' %(data_folder))

    # Program specific results directories
    data_trimmed_dir = "%s/trimmed" % (data_folder)
    fastqc_dir = "%s/fastqc_results" % (data_folder)

    results_dir = "%s/results_STAR" %(results_folder)
    rsem_results_dir = "%s/results_RSEM" %(results_folder) 
    #htseq_dir = "%s/htseqcounts" % (results_dir)
    

    # Run create directories function to create directory structure
    create_dirs(data_trimmed_dir, fastqc_dir, results_dir, rsem_results_dir)

    print("FIRST_PAIR_FILES: ", first_pair_files)

    # Loop through each file and create filenames
    file_count = 1
    for first_pair_file in first_pair_files:
        first_file_name_full = first_pair_file.split('/')[-1]

        if re.search("_R1_001.fastq.gz", first_pair_file):
            print("Found a match for: '_R1_001.fastq.gz'\n")
            second_pair_file = re.sub("_R1_001.fastq.gz", "_R2_001.fastq.gz", first_pair_file)
        elif re.search(".R1.fastq.gz", first_pair_file):
            print("Found a match for: '.R1.fastq.gz'\n")
            second_pair_file = re.sub(".R1.fastq.gz", ".R2.fastq.gz", first_pair_file)
        elif re.search("_1.fastq.gz", first_pair_file):
            print("Found a match for: '_1.fastq.gz'\n")
            second_pair_file = re.sub("_1.fastq.gz", "_2.fastq.gz", first_pair_file)        
        elif re.search(".1.fastq.gz", first_pair_file):
            print("Found a match for: '.1.fastq.gz'\n")
            second_pair_file = re.sub(".1.fastq.gz", ".2.fastq.gz", first_pair_file) 
        else:
            print("No match found for any defined types\n")
            with open(error_file,'w') as f:
                f.write("No match found for any defined types\n")
                f.write(str(first_pair_file))
            f.close()   

        
        print(first_pair_file)
        print(second_pair_file)

        second_file_name_full = second_pair_file.split('/')[-1]
        file_ext = first_pair_file.split('.')[-1]
        print()
        print ('\033[32m Processing File: %s of %s (%s)\033[0m' %(file_count, len(first_pair_files), first_file_name_full ))

        first_file_name = re.split('.fq|.fq.gz',first_file_name_full)[0]
        second_file_name = re.split('.fq|.fq.gz',second_file_name_full)[0]
        print('first_file_name:%s, second_file_name:%s' %(first_file_name,second_file_name))

        # Collect Sample attributes
        exp_name = folder_name
        print("exp_name: %s" %(exp_name))
        lane = first_file_name.split("_")[-1]
        print("Lane: %s" %(lane))
        sample_id = re.split('.fq|.fq.gz', first_file_name)[0]
        print("sample_id: %s" %(sample_id))
        #sys.exit()

        # Run TrimGalore
        trim_galore(first_pair_file,second_pair_file,folder_name,sample_id,file_ext,data_trimmed_dir,fastqc_dir)
        file_count = file_count + 1
        #sys.exit()

    # Collect Trimmed data for input into STAR
    first_pair_group,second_pair_group,star_input_files = collect_trimmed_data(data_trimmed_dir,file_ext)

    # Run STAR
    run_star(first_pair_group,second_pair_group,results_dir,star_input_files, folder_name, genome_dir)

    # Run Salmon Quant
    #run_salmon_quant(results_dir, folder_name, genome_fasta)

    # Run RSEM Quantification
    run_RSEM(results_dir,rsem_results_dir,folder_name)

    # Run HTSeq count
    #run_htseq(htseq_dir, results_dir, folder_name, genome_gff)

    folder_count += 1

    return data_trimmed_dir, fastqc_dir, results_dir


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    #parser.add_argument('genomedir', help='genome directory')
    parser.add_argument('dataroot', help="parent of input directory")
    parser.add_argument('indir', help="input directory (R<somenumber>)")
    parser.add_argument('outdir', help='output directory')
    parser.add_argument('--starPrefix', help="STAR output file name prefix")
    parser.add_argument('--salmonPrefix', help="Salmon output folder name prefix")

    args = parser.parse_args()

    now = datetime.datetime.now()
    timeprint = now.strftime("%Y-%m-%d %H:%M")
    data_folder = "%s/%s" % (args.dataroot, args.indir)

    genome_dir = "/proj/omics4tb2/SYGNAL/reference_genome/homo_sapiens"
    genome_fasta = "%s/STAR_reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa" %genome_dir
    genome_gff = "%s/STAR_reference/Homo_sapiens.GRCh38.99.gtf" %genome_dir
    #tt = glob.glob('/proj/omics4tb2/SYGNAL/XCures/patients_raw_data_SNO/*/*/RNA/*[_.R]1*.fastq.gz',recursive=True)
    #print(tt)
    #sys.exit()

    create_genome_index(genome_dir, genome_fasta)
    data_trimmed_dir,fastqc_dir,results_dir = run_pipeline(data_folder, args.outdir, genome_dir, genome_fasta, genome_gff)
