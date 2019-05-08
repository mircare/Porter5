""" split fasta file with multiple sequences into multiple fasta files with 1 sequence each"""
import os
import sys

# check arguments
if len(sys.argv) < 2:
    print("\n---->>You need to specify where the big fasta file is")
    exit()


# catch big fasta
fasta = open(sys.argv[1], "r").readlines()

# create new directory to save outcome
os.mkdir("Fastas")
os.chdir("Fastas")

# fix formatting
i = 0
while i < len(fasta):
    pid = fasta[i].replace(">", "").split()[0]
    if os.path.isfile(pid):
        j = 0
        while os.path.isfile(pid+str(j)):
            j += 1
        pid += str(j) 
    i+=1
    aa = ">"+pid+"\n"
    while i < len(fasta) and fasta[i][0] != ">":
        aa = aa + fasta[i].strip()
        i+=1
    f = open(pid,"w")
    f.write(aa+"\n")

f.close()

print("Done. You can find the splitted fasta files inside Fastas/. Please run multiple_fasta.py on it to start predicting.")