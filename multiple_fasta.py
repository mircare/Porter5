# Usage: python3 multiple_fasta.py -i Fastas/ --cpu 4 --parallel 2
# (run 2 parallel instances of Porter5 on 4 cores - total of 8 cores)

import os, sys
import argparse
from multiprocessing import Pool

### parallel code ##
def loop(line):
    if args.fast and args.tmp:
        os.system('python3 %s -i %s --cpu %d --fast --tmp' % (executable, line, args.cpu))
    elif args.fast:
        os.system('python3 %s -i %s --cpu %d --fast' % (executable, line, args.cpu))
    elif args.tmp:
        os.system('python3 %s -i %s --cpu %d --tmp' % (executable, line, args.cpu))
    else:
        os.system('python3 %s -i %s --cpu %d' % (executable, line, args.cpu))

### set argparse
parser = argparse.ArgumentParser(description="This is the standalone of Porter5 for multiple inputs. It is sufficient to specify a directory containing FASTA files to start the prediction of their Secondary Structure in 3- and 8-classes. It is also possible to run multiple predictions in parallel (TOTAL cpu = --cpu x --parallel). Please run Porter5.py if you have only 1 protein sequence to predict.", 
epilog="E.g., to run 2 instances of Porter5 on 4 cores (total of 8 cores): python3 multiple_fasta.py -i Fastas/ --cpu 4 --parallel 2")
parser.add_argument("-i", type=str, nargs=1, help="Indicate the directory containing the FASTA files.")
parser.add_argument("--cpu", type=int, default=1, help="Specify how many cores to assign to each prediction.")
parser.add_argument("--parallel", type=int, default=1, help="Specify how many instances to run in parallel.")
parser.add_argument("--fast", help="Use only HHblits (skipping PSI-BLAST) to perform a faster prediction.", action="store_true")
parser.add_argument("--tmp", help="Leave output files of HHblits and PSI-BLAST, i.e. log, hhr, psi, chk, and blastpgp files.", action="store_true")
parser.add_argument("--setup", help="Initialize Porter5 from scratch. Run it when there has been any change involving PSI-BLAST, HHblits, Porter itself, etc.", action="store_true")
args = parser.parse_args()

## check arguments
if not args.i:
    print("Usage: python3 "+sys.argv[0]+" -i <fasta_dir> [--cpu CPU_number] [--parallel instances] [--fast]\n--help for the full list of commands")
    exit()

#initialization variables
executable = os.path.abspath(os.path.dirname(sys.argv[0]))+"/Porter5.py"
if not os.path.isfile(executable):
    print("\n---->>No executable retrieved at", executable)
    exit()

if not os.path.isdir("".join(args.i)):
    print("\n---->>","".join(args.i),"isn't a directory! Please consider running split_fasta.py.")
    exit()

if not os.path.isfile(os.path.abspath(os.path.dirname(sys.argv[0]))+"/scripts/config.ini") or args.setup:
    os.system("python3 %s --setup" % executable)

# fetch all the inputs from the passed directory, and sort them by size
os.chdir("".join(args.i))
sorted_files = sorted(os.listdir(os.getcwd()), key = os.path.getsize, reverse=True)

#ligth the bomb // launch the parallel code
with Pool(args.parallel) as p:
    p.map(loop, sorted_files, 1)
