#importing necessary packages
from collections import defaultdict
import gzip
from Bio.Align import PairwiseAligner
import numpy as np
#import matplotlib.pyplot as plt
import Bio.Seq
import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns
import Bio.Seq
import Bio.SeqIO
from Bio.SeqRecord import SeqRecord
import os
import sys
from pathlib import Path

#-------Loading in system args
# Parse command-line arguments to read in files of interest
if len(sys.argv) != 7:
    print("Usage: python3 sensor_extraction.py <input_df> <R1_FILE> <R2_FILE> splitby -o <folder_name>")
    sys.exit(1)

input_df = pd.read_csv(Path(sys.argv[1]))
R1_FILE= Path(sys.argv[2])
R2_FILE= Path(sys.argv[3])
splitby = str(sys.argv[4])
folder_name = str(sys.argv[6])

#-----GLOBAL VARIABLES------
MIN_QUALITY = 30 #originally set to 30
GZ=False

#THIS FUNCTION MAY NEED TO BE UPDATED DEPENDING ON HOW THE DATA IS STRUCTURED
def fastq_reader(fname, gz=False):
    _open = gzip.open if gz else open
    proc_read = (lambda line: line.strip().decode()) if gz else (lambda line: line.strip())
    proc_qual = (lambda line: line.strip()) if gz else (lambda line: line.strip().encode('utf-8'))
    with _open(fname, 'r') as f:
        for i, line in enumerate(f):
            if i % 4 == 0:
                identifier = proc_read(line)
            elif i % 4 == 1:
                read = proc_read(line)
            elif i % 4 == 3:
                quality = proc_qual(line)
                yield identifier[-20:], read, quality #protospacer (index2), read, quality score


def quality_checker(q1, q2):
    """ 
    Checks quality score of the reads
    """

    qual1 = (np.frombuffer(q1, dtype=np.uint8) - 33).mean()
    qual2 = (np.frombuffer(q2, dtype=np.uint8) - 33).mean()

    low_qual = []

    if qual1 < MIN_QUALITY:
        low_qual.append('1')
    if qual2 < MIN_QUALITY:
        low_qual.append('2')
    
    if len(low_qual)==0:
        quality_out = 'good quality'
    else:
        quality_out = 'low_qual_r' + ''.join(low_qual)

    return quality_out

'''Prepare fastq records.'''
def to_IOSeq_rec(sensor_seq, i, q1):
    """ 
    Modified to only extract 42 nt sensor
    After 15 nt barcode
    """
    qname = 'read_' + str(i)
    record = SeqRecord(Bio.Seq.Seq(sensor_seq[15:15+42]), id=qname, name=qname, description='', dbxrefs=[])

    #add quality score to the record
    qual = (np.frombuffer(q1, dtype=np.uint8) - 33)
    record.letter_annotations["phred_quality"] = qual[15:15+42]
    
    return record

def extraction_filtration(folder_name, input_df,R1_FILE,R2_FILE, splitby='barcode', breakpoint=False):

    """
    Takes in reads and returns dataframe containing pegRNA counts AND sensor outcomes.

    Also returns summary of outcomes including:
    - low quality (and which of the reads are low quality (or if all are low quality))
    - no extension match
    - no protospacer match
    - decoupled extension-protospacer
    - correct identification
    """

    os.mkdir(folder_name)#make the new directory (needs to be a non-existent folder)

    #write empty fastq files to a new folder for EACH gRNA
    #need to modify based on where the guide_ID column is/what its name is
    for i in list(input_df['gRNA_id']):
        f_name = i + '.fastq'
        with open(folder_name + '/' + f_name, 'w') as fp:
            pass

    #also an array of the extension sequences
    #ext_sequences = np.array(list(input_df['PBS_RTT_5to3']))
    #sensor_handles = np.unique(input_df['sensor_handle'])

    #lists and dicts for matching up reads
    proto_list = [i[1:] for i in input_df['Protospacer']]
    bc_list = [str(Bio.Seq.Seq(i).reverse_complement()) for i in input_df['Hamming_BC']]

    proto_bc_dict = dict(zip(bc_list, proto_list))
    confusion_mat = np.zeros((len(proto_list), len(bc_list)))

    
    #------initialize a dataframe for holding the pegRNA counts and sensor outcomes
    d1 = pd.DataFrame(dict(zip(['Guide_ID', 'sgRNA_no_Gstart', 'unique_BC'], [list(input_df['gRNA_id']), proto_list, bc_list])))

    cols = ['guide_count', 'matched_sensor_count', 'bc_count']
    z = np.zeros((len(cols), len(input_df)))
    d2 = pd.DataFrame(dict(zip(cols, z)))
    df1 = pd.concat((d1, d2), axis=1).set_index('Guide_ID')

    #and annotating the duplicate guides
    i, v = np.unique(d1['sgRNA_no_Gstart'], return_counts=True)

    duplicate_dict = dict(zip(i, v>1))
    df1['duplicate_sgRNA'] = [duplicate_dict[i] for i in df1['sgRNA_no_Gstart']]
    df1[df1['duplicate_sgRNA']==True]

    #------initialize a dataframe for holding metadata about the identification of sensors
    outcomes = ['good quality', 'low_qual_r1', 'low_qual_r2', 'low_qual_r12', 
                'proto_identified',
                'bc_identified',
                'no_match_proto',
                'no_match_bc',
                'mismatch',
                'correct_match'
                ]
    
    outcomes_count = np.zeros(len(outcomes))
    class_df = pd.DataFrame(dict(zip(['classification', 'count'],[outcomes, outcomes_count]))).set_index('classification')

    #iterating through the reads...
    for i, ((proto1, r1, q1), (proto2, r2, q2)) in enumerate(zip(fastq_reader(R1_FILE, gz=GZ), fastq_reader(R2_FILE, gz=GZ)), 1):
        #r1 = sensor read
        #r2 = extension read
        #proto1/2 = protospacer read (index2)

        #first check the quality
        quality_out = quality_checker(q1,q2)
        class_df.loc[quality_out, 'count']+=1

        #if quality only good for r2, use it to identify the protospacer
        if quality_out == 'low_qual_r1':
            proto = r2[:20]
            
            #not including these in the metadata count...
            #just the guide count (separate thing)
            if proto in proto_list:
                df1.loc[df1['sgRNA_no_Gstart']==proto, 'guide_count']+=1
                
            else:
                continue
            
        elif quality_out=='good quality':
            proto = r2[:20]
            bc = r1[:15]

            if proto in proto_list:
                #guide counts
                df1.loc[df1['sgRNA_no_Gstart']==proto, 'guide_count']+=1
                class_df.loc['proto_identified', 'count']+=1

            elif proto not in proto_list:
                class_df.loc['no_match_proto', 'count']+=1

            if bc in bc_list:
                df1.loc[df1['unique_BC']==bc, 'bc_count']+=1
                class_df.loc['bc_identified', 'count']+=1

            elif bc not in bc_list:
                class_df.loc['no_match_bc', 'count']+=1


            if proto in proto_list:

                if bc in bc_list:
                    bc_idx = bc_list.index(bc)
                    proto_match = proto_bc_dict[bc]
                    
                    if proto == proto_match:
                        #correct_match+=1
                        class_df.loc['correct_match', 'count']+=1

                        #and add in the matched_sensor_count
                        df1.loc[df1['unique_BC']==bc, 'matched_sensor_count']+=1

                        #and annotate the confusion matrix
                        confusion_mat[bc_idx][bc_idx]+=1

                        #-----ADD IN CODE FOR PUTTING THINGS IN SEPARATE FASTQ FILE HERE------
                        #add the sensor read to the appropriate fastq file (numbered according to the pegRNA index)
                        guide_id = df1[df1['unique_BC']==bc].index[0]

                        out_file = folder_name + '/' +  guide_id + '.fastq'

                        record = to_IOSeq_rec(r1, i, q1)

                        with open(out_file, 'a') as fq:
                            #and write it to the appropriate file
                            Bio.SeqIO.write(record, fq, 'fastq')

                    else:
                        #mismatch+=1
                        class_df.loc['mismatch', 'count']+=1

                        proto_idx = proto_list.index(proto)
                        confusion_mat[proto_idx][bc_idx]+=1

            #if splitting up reads by barcode only, add additional reads THAT AREN"T MISMATCHED to the fastq file for sensor analysis
            if splitby=='barcode':
                if proto not in proto_list:

                    if bc in bc_list:
                        
                        guide_id = df1[df1['unique_BC']==bc].index[0]

                        out_file = folder_name + '/' +  guide_id + '.fastq'

                        record = to_IOSeq_rec(r1, i, q1)

                        with open(out_file, 'a') as fq:
                            #and write it to the appropriate file
                            Bio.SeqIO.write(record, fq, 'fastq')

                    else:
                        continue
            else:
                continue

        if breakpoint != False:
            if i>breakpoint:
                break

    
    #and then prune the duplicates...
    #both the files and the counts table...

   
    return df1, class_df, confusion_mat

#---RUNNING THE SCRIPT------
count_df, class_df, mat = extraction_filtration(folder_name, input_df, R1_FILE,R2_FILE,splitby, breakpoint=False)

#save the counts and classification dataframes
count_df.to_csv('counts/' + folder_name + '_count_df.csv')
class_df.to_csv('classification/' + folder_name + '_classification_df.csv')

confusion_mat = pd.DataFrame(mat)
confusion_mat.to_csv('confusion_mats/' + folder_name + '_confusion_matrix.csv')