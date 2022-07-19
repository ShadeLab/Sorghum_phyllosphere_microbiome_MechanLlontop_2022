import sys
import os
import subprocess
import shlex
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

fn = sys.argv[1]
if os.path.exists(fn):
    print ("Directory of annotated genomes: ",os.path.dirname(fn))
else:
    print("Directory",fn, "do not exist\n\n ")

current_direc = os.getcwd()
sequences_16S_d = '16S_sequences'
sequences_16S_d_path = os.path.join(current_direc,sequences_16S_d)
if not sequences_16S_d in os.listdir(current_direc):
    os.mkdir(sequences_16S_d_path)

direct = fn
os.chdir(direct)

for d in os.listdir(direct):
    if os.path.isdir(d) and d != '2585428185':
        print(d)
        d_path = os.path.join(direct,d)
        for f in os.listdir(d_path):
            if f.endswith('.gbk'):
                os.chdir(d_path)
                fasta_records = []
                name_fasta_16S = d+'_16S.fasta'
                for rec in SeqIO.parse(f,'gb'):
                    print(len(rec.id))
                    for feat in rec.features:
                        if feat.type == 'rRNA' and feat.qualifiers['product'][0] == '16S ribosomal RNA':
                            ffn_id  = feat.qualifiers["locus_tag"][0]
                            ffn_des = feat.qualifiers["product"][0]
                            ffn_seq  = feat.location.extract(rec).seq
                            ffn_rec = SeqRecord(Seq(str(ffn_seq)),\
                                                id = str(ffn_id),\
                                                description = str(ffn_des))
                            fasta_records.append(ffn_rec)
                os.chdir(sequences_16S_d_path)
                SeqIO.write(fasta_records,name_fasta_16S,'fasta')
                os.chdir(direct)
