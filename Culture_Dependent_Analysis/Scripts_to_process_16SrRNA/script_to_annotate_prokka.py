#!usr/bin/env python3
import sys
import os
import subprocess
import shlex

fn = sys.argv[1]
if os.path.exists(fn):
    print ("Directory of genomes to anotate: ",os.path.dirname(fn))
else:
    print("Directory",fn, "do not exist\n\n ")

direct = fn
os.chdir(direct)

for file in os.listdir(direct):
    if file.endswith(".fna"):
        print(file)
        in_file  = str(file)
        out_file = str(in_file[:-4])
        
        outdir = fn+out_file+"/"
        strain_tag = in_file.split(".")[0]
        outdir = fn+out_file+"/"
        print('#####',outdir)

        print("annotating genome: ",strain_tag)

        to_comand_line = str('prokka '+in_file+' --outdir '+outdir+' --locustag '+strain_tag)
        comand_line = to_comand_line
        subprocess.call(comand_line, shell=True)
        
        print(strain_tag, " ... Finished")
    else:
        print("ERROR ... ",file,"\t could not identify annotation format")
