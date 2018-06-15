# Porter 5 
### Light, fast and high quality prediction of protein secondary structure in 3 and 8 classes

The web server of Porter 5 is available at http://distilldeep.ucd.ie/porter/.  
The train and test sets are available at http://distilldeep.ucd.ie/porter/data/.


### Reference
Porter 5: fast, state-of-the-art ab initio prediction of protein secondary structure in 3 and 8 classes<br>
Mirko Torrisi, Manaz Kaleel and Gianluca Pollastri; bioRxiv 289033; doi: https://doi.org/10.1101/289033.


## Setup
```
$ git clone https://github.com/mircare/Porter5/
```
or
```
$ wget http://distilldeep.ucd.ie/porter/data/Porter5.tgz
$ tar zxvf Porter5.tgz
```

### Requirements
1. Python3 (https://www.python.org/downloads/);
1. NumPy (https://www.scipy.org/scipylib/download.html);
1. HHblits (https://github.com/soedinglab/hh-suite/);
1. uniprot20 (http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/old-releases/uniprot20_2016_02.tgz).

#### Optionally (for more accurate predictions):
1. PSI-BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/); 
1. UniRef90 (ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz).


## How to run Porter 5 with/without PSI-BLAST
```
# To exploit HHblits only (for fast and accurate predictions)
$ python3 Porter5/Porter5.py -i Porter5/example/2FLGA.fasta --cpu 4 --fast

# To exploit both PSI-BLAST and HHblits (for very accurate predictions)
$ python3 Porter5/Porter5.py -i Porter5/example/2FLGA.fasta --cpu 4
```

## Performances on large independent test set
| Method | Q3 per AA | SOV per AA | Q3 per protein | SOV per protein |
| :--- | :---: | :---: | :---: | :---: |
| **Porter 5** | **83.81%** | **83.29%** | **84.32%** | **84.57%** |
| SPIDER 3 | 83.15% | 82.04% | 83.42% | 83.17% |
| **Porter 5 *HHblits only*** | **83.06%** | **82.21%** | **83.68%** | **83.71%** |
| SSpro 5 *with templates* | 82.58% | 80.13% | 83.94% | 82.49% |
| PSIPRED 4.01 | 81.88% | 80.17% | 82.48% | 81.70% |
| RaptorX-Property | 81.86% | 81.86% | 82.57% | 83.13% |
| Porter 4 | 81.66% | 81.11% | 82.29% | 82.51% | 
| SSpro5 *ab initio* | 81.17% | 78.54% | 81.10% | 79.45% |
| DeepCNF | 81.04% | 80.84% | 81.16% | 81.46% |

Reference: Table 1 in https://doi.org/10.1101/289033.

## License
This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

Email us at gianluca[dot]pollastri[at]ucd[dot]ie if you wish to use it for purposes not permitted by the CC BY-NC-SA 4.0.

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a>
