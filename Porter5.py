import argparse, configparser
import os, sys
import time

## check Python version
if sys.version_info[0] < 3:
   print("Python2 detected, please use Python3.")
   exit()
   
### set argparse
parser = argparse.ArgumentParser(description="This is the standalone of Porter5. It is sufficient to run it on a FASTA file to start the prediction of its Secondary Structure in 3- and 8-classes.\
    Please run multiple_fasta.py if you have multiple protein sequences to predict, or split_fasta.py if you have one fasta file with multiple protein sequences.", 
epilog="E.g., to run Porter on 4 cores: python3 Porter5.py -i example/2FLGA.fasta --cpu 4")
parser.add_argument("-i", type=str, nargs=1, help="Indicate the FASTA file containing the protein to predict.")
parser.add_argument("--cpu", type=int, default=1, help="Specify how many cores to assign to this prediction.")
parser.add_argument("--fast", help="Use only HHblits (skipping PSI-BLAST) to perform a faster prediction.", action="store_true")
parser.add_argument("--tmp", help="Leave output files of HHblits and PSI-BLAST, i.e. log, hhr, psi, chk, and blastpgp files.", action="store_true")
parser.add_argument("--setup", help="Initialize Porter5 from scratch. Run it when there has been any change involving PSI-BLAST, HHblits, Porter itself, etc.", action="store_true")
args = parser.parse_args()

## check arguments
if not args.i:
    print("Usage: python3 "+sys.argv[0]+" -i <fasta_file> [--cpu CPU_number] [--fast]\n--help for the full list of commands")
    exit()

# save protein path and name, and current PATH
filename = "".join(args.i)
path = os.path.abspath(sys.argv[0].replace("Porter5.py","scripts/"))
predict = path+"/Predict_BRNN/Predict"
models = path+"/Predict_BRNN/models/"

## PSI-BLAST and HHblits variables and paths.
config = configparser.ConfigParser()
if not os.path.exists(path+"/config.ini") or args.setup:
    psiblast = input("Please insert the absolute path to psiblast (e.g., /home/username/psiblast): ")
    uniref90 = input("Please insert the absolute path to uniref90 (e.g., /home/username/UniProt/uniref90.fasta): ")
    hhblits = input("Please insert the call to HHblits (e.g., hhblits): ")
    uniprot20 = input("Please insert the absolute path to uniprot20 - DATABASE NAME INCLUDED (e.g., /home/username/uniprot20_2016_02/uniprot20_2016_02): ")
    
    config['DEFAULT'] = {'psiblast': psiblast,
                    'uniref90': uniref90,
                    'hhblits': hhblits,
                    'uniprot20' : uniprot20}

    with open(path+"/config.ini", 'w') as configfile:
        config.write(configfile)

    # compile predict and set absolute paths for all model files
    os.system('cd %s/Predict_BRNN; make -B; bash set_models.sh; cd %s' % (path, path))
    
    print("\n>>>>> Setup completed successfully. If you encounter any problems in future, please run \"python3 Porter5.py --setup\". <<<<<\n")
else:
    config.read(path+"/config.ini")


print("~~~~~~~~~ Prediction of "+filename+" started ~~~~~~~~~")

time0 = time.time()
if not args.fast:
    ### run PSI-BLAST and process output
    os.system('%s -query %s -out_pssm %s.chk -num_threads %d -dbsize 0 -num_alignments 300000 -num_iterations 2 -evalue 0.001 -inclusion_ethresh 1e-3 -pseudocount 10 -comp_based_stats 1 -db %s >> %s.log' % (config['DEFAULT']['psiblast'], filename, filename, args.cpu, config['DEFAULT']['uniref90'], filename))
    os.system('%s -in_pssm %s.chk -out %s.blastpgp -num_threads %d -dbsize 0 -num_alignments 300000 -num_iterations 1 -evalue 0.001 -inclusion_ethresh 1e-3 -comp_based_stats 1 -db %s >> %s.log 2>&1' % (config['DEFAULT']['psiblast'], filename, filename, args.cpu, config['DEFAULT']['uniref90'], filename))
    os.system('%s/process-blast.pl %s.blastpgp %s.flatblast %s' % (path, filename, filename, filename))

    time1 = time.time()
    print('PSI-BLAST executed in %.2fs' % (time1-time0))
else:
    time1 = time.time()

#### run HHblits and process output
os.system('%s -d %s -i %s -opsi %s.psi -cpu %d -n 3 -maxfilt 150000 -maxmem 27 -v 2 2>> %s.log >> %s.log' % (config['DEFAULT']['hhblits'], config['DEFAULT']['uniprot20'], filename, filename, args.cpu, filename, filename))
os.system('%s/process-psi.sh %s.psi' % (path, filename))

time2 = time.time()
print('HHblits executed in %.2fs' % (time2-time1))

### encode alignments made with HHblits or PSI-BLAST
os.system('python3 %s/process-alignment.py %s.flatpsi flatpsi %d' % (path, filename, args.cpu)) # generated with HHblits
aa = list("".join(line.strip() for line in open(filename, "r").readlines()[1:])) # get plain list of AA from FASTA
length = len(aa)

flatpsi_ann = open(filename+".flatpsi.ann", "r").readlines()
if not args.fast:
    os.system('python3 %s/process-alignment.py %s.flatblast flatblast %d' % (path, filename, args.cpu)) # generated with PSI-BLAST

    ## concatenate the previous 2
    flatblast_ann = open(filename+".flatblast.ann", "r").readlines()

    # write header, protein name, and length
    flatblastpsi_ann = open(filename+".flatblastpsi.ann", "w")
    flatblastpsi_ann.write("1\n44 3\n"+ filename +"\n"+ str(length) +"\n")

    # concatenate
    tmp = flatblast_ann[4].strip().split(" ")
    tmp0 = flatpsi_ann[4].split(" ")

    for j in range(length):
        x = j*22
        tmp.insert(x + 22 + j, " ".join(tmp0[x:x+22]))
    flatblastpsi_ann.write(" ".join(tmp))
    flatblastpsi_ann.close()

time3 = time.time()
print('Alignments encoded in %.2fs' % (time3-time2))

### predict SS in 3 classes
os.system('%s %smodelv8_ss3 %s.flatpsi.ann > /dev/null' % (predict, models, filename))
if not args.fast:
    os.system('%s %smodelv7_ss3 %s.flatblast.ann > /dev/null' % (predict, models, filename))
    os.system('%s %smodelv78_ss3 %s.flatblastpsi.ann > /dev/null' % (predict, models, filename))

time4 = time.time()
print('Prediction in 3 classes made in %.2fs' % (time4-time3))


### ensemble predictions and process output
secstruc = {0 : "H", 1 : "E", 2 : "C"}
SS = [0] * 3
prediction = open(filename+".ss3", "w")
prediction.write("#\tAA\tSS\tHelix\tSheet\tCoil\n")

prob_hh = list(map(float, open(filename+".flatpsi.ann.probsF", "r").readlines()[3].split()))
if not args.fast:
    prob_psi = list(map(float, open(filename+".flatblast.ann.probsF", "r").readlines()[3].split()))
    prob_psihh = list(map(float, open(filename+".flatblastpsi.ann.probsF", "r").readlines()[3].split()))

    for i in range(length):
        for j in range(3):
            SS[j] = round((3*prob_psi[i*3+j]+3*prob_hh[i*3+j]+prob_psihh[i*3+j])/7, 4)
        index = SS.index(max(SS))
        prediction.write(str(i+1)+"\t"+aa[i]+"\t"+secstruc[index]+"\t"+str(SS[0])+"\t"+str(SS[1])+"\t"+str(SS[2])+"\n")
else:
    for i in range(length):
        for j in range(3):
            SS[j] = round(prob_hh[i*3+j], 4)
        index = SS.index(max(SS))
        prediction.write(str(i+1)+"\t"+aa[i]+"\t"+secstruc[index]+"\t"+str(SS[0])+"\t"+str(SS[1])+"\t"+str(SS[2])+"\n")
prediction.close()


######## eight-state prediction ########
def generate8statesANN(extension, prob, ann):
    input_size = int(ann[1].split()[0])

    ss3 = open(filename+"."+extension+".ann+ss3", "w")
    ss3.write("1\n"+str(input_size+3)+" 8\n"+ filename +"\n"+ str(length) +"\n")

    prob = list(map(str, prob))
    tmp = ann[4].strip().split(" ")

    for j in range(length):
        ss3.write(" ".join(tmp[j*input_size:(j+1)*input_size])+" "+" ".join(prob[j*3:j*3+3])+" ")
    ss3.close()

### generate inputs
generate8statesANN("flatpsi", prob_hh, flatpsi_ann)
if not args.fast:
    generate8statesANN("flatblast", prob_psi, flatblast_ann)

    flatblastpsi_ann = open(filename+".flatblastpsi.ann", "r").readlines()
    generate8statesANN("flatblastpsi", prob_psihh, flatblastpsi_ann)


### predict in 8 classes
os.system('%s %smodelv8_ss8 %s.flatpsi.ann+ss3 > /dev/null' % (predict, models, filename))
if not args.fast:
    os.system('%s %smodelv7_ss8 %s.flatblast.ann+ss3 > /dev/null' % (predict, models, filename))
    os.system('%s %smodelv78_ss8 %s.flatblastpsi.ann+ss3 > /dev/null' % (predict, models, filename))
time5 = time.time()
print('Prediction in 8 classes made in %.2fs' % (time5-time4))


### ensemble predictions and process output
secstruc = {0 : "G", 1 : "H", 2 : "I", 3 : "E", 4 : "B", 5 : "C", 6 : "S", 7 : "T"}
SS = [0] * 8

prediction = open(filename+".ss8", "w")
prediction.write("#\tAA\tSS\tG\tH\tI\tE\tB\tC\tS\tT\n")

prob_hh = list(map(float, open(filename+".flatpsi.ann+ss3.probsF", "r").readlines()[3].split()))
if not args.fast:
    prob_psi = list(map(float, open(filename+".flatblast.ann+ss3.probsF", "r").readlines()[3].split()))
    prob_psihh = list(map(float, open(filename+".flatblastpsi.ann+ss3.probsF", "r").readlines()[3].split()))

    for i in range(length):
        for j in range(8):
            SS[j] = round((3*prob_psi[i*8+j]+3*prob_hh[i*8+j]+prob_psihh[i*8+j])/7, 4)
        index = SS.index(max(SS))
        prediction.write(str(i+1)+"\t"+aa[i]+"\t"+secstruc[index]+"\t"+str(SS[0])+"\t"+str(SS[1])+"\t"+str(SS[2])+"\t"+str(SS[3])+"\t"+str(SS[4])+"\t"+str(SS[5])+"\t"+str(SS[6])+"\t"+str(SS[7])+"\n")
else:
    for i in range(length):
        for j in range(8):
            SS[j] = round(prob_hh[i*8+j], 4)
        index = SS.index(max(SS))
        prediction.write(str(i+1)+"\t"+aa[i]+"\t"+secstruc[index]+"\t"+str(SS[0])+"\t"+str(SS[1])+"\t"+str(SS[2])+"\t"+str(SS[3])+"\t"+str(SS[4])+"\t"+str(SS[5])+"\t"+str(SS[6])+"\t"+str(SS[7])+"\n")
prediction.close()

timeEND = time.time()
print('Porter5 executed on %s in %.2fs (TOTAL)' % (filename, timeEND-time0))

### remove temporary files
os.system('rm -f %s.flatblast.ann+ss3.probs %s.flatpsi.ann.probs %s.flatblast.ann.probs %s.flatblastpsi.ann+ss3.probsF %s.flatblast.ann+ss3 %s.flatblastpsi.ann.probsF %s.flatblastpsi.ann %s.flatpsi.ann+ss3.probsF %s.flatpsi.ann %s.flatblast %s.flatblastpsi.ann+ss3.probs %s.flatblastpsi.ann+ss3 %s.flatblastpsi.ann.probs %s.flatpsi %s.flatpsi.ann+ss3.probs %s.flatblast.ann %s.flatblast.ann+ss3.probsF %s.flatpsi.ann+ss3 %s.flatpsi.ann.probsF %s.flatblast.app %s.flatblast.ann.probsF 2> /dev/null' % (filename, filename, filename, filename, filename, filename, filename, filename, filename, filename, filename, filename, filename, filename, filename, filename, filename, filename, filename, filename, filename))
if not args.tmp:
    os.system('rm -f %s.blastpgp, %s.chk, %s.hhr, %s.psi, %s.log 2> /dev/null' % (filename, filename, filename, filename, filename))