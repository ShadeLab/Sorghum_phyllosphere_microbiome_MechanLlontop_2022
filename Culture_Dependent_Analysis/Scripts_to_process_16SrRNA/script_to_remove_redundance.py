import sys
import os
import csv
import subprocess
import shutil
from Bio import SeqIO

fn = sys.argv[1]
if os.path.exists(fn):
    print ("working directory: ",os.path.dirname(fn))
else:
    print("Directory",fn, "do not exist\n\n ")

current_direc = os.getcwd()
sequences_nr_d = '16S_sequences_nr'
sequences_nr_d_path = os.path.join(current_direc,sequences_nr_d)
if not os.path.isdir(sequences_nr_d_path):
    os.mkdir(sequences_nr_d_path)

direct = fn
os.chdir(direct)

general_table = []
for file in os.listdir(direct):
    if file.endswith(".fasta"):
        os.chdir(direct)
        in_file  = str(file)
        total_16_before = 0
        for rec in SeqIO.parse(in_file,'fasta'):
            total_16_before += 1

        out_file = str(in_file.replace('.fasta','-nr.fasta'))
        to_comand_line = str('cd-hit -i '+in_file+' -o '+out_file+' -aL 0.9 -c 0.9')
        subprocess.call(to_comand_line, shell=True)

        total_16_after = 0
        for rec in SeqIO.parse(out_file,'fasta'):
            total_16_after += 1

        source = os.path.join(direct,out_file)
        destin = os.path.join(sequences_nr_d_path,out_file)
        shutil.copyfile(source,destin)
        os.remove(out_file)
        os.remove(out_file+'.clstr')
        general_table.append([file , total_16_before , total_16_after])
general_table.insert(0,['genome','total_16_before' , 'total_16_after'])

filename_table = os.path.join(current_direc,'General_table_cdhit_analysis.csv')
open_table_table =  open(filename_table,'w')
writ_table_table = csv.writer(open_table_table,delimiter='\t')
writ_table_table.writerows(general_table)
open_table_table.close()
