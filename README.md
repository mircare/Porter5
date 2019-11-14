[![PWC](https://img.shields.io/endpoint.svg?url=https://paperswithcode.com/badge/deeper-profiles-and-cascaded-recurrent-and/protein-secondary-structure-prediction-on-2)](https://paperswithcode.com/sota/protein-secondary-structure-prediction-on-2?p=deeper-profiles-and-cascaded-recurrent-and)

# Porter 5 
### Light, fast and high quality prediction of protein secondary structure in 3 and 8 classes

The web server of Porter 5 is available at http://distilldeep.ucd.ie/porter/.   
The train and test sets are available at http://distilldeep.ucd.ie/porter/data/.

More protein structure annotations predicted at http://distilldeep.ucd.ie/brewery/.

### References
Deeper Profiles and Cascaded Recurrent and Convolutional Neural Networks for state-of-the-art Protein Secondary Structure Prediction, Scientific Reports, Nature Publishing Group<br>
Mirko Torrisi, Manaz Kaleel and Gianluca Pollastri; doi: https://doi.org/10.1038/s41598-019-48786-x.

Porter 5: fast, state-of-the-art ab initio prediction of protein secondary structure in 3 and 8 classes<br>
Mirko Torrisi, Manaz Kaleel and Gianluca Pollastri; bioRxiv 289033; doi: https://doi.org/10.1101/289033.

Protein Structure Annotations; Essentials of Bioinformatics, Volume I. Springer Nature<br>
Mirko Torrisi and Gianluca Pollastri; doi: https://doi.org/10.1007/978-3-030-02634-9_10.


## Setup
```
$ git clone https://github.com/mircare/Porter5/ --depth 1 && rm -rf Porter5/.git
```

### Requirements
1. Python3 (https://www.python.org/downloads/);
1. NumPy (https://www.scipy.org/scipylib/download.html);
1. HHblits (https://github.com/soedinglab/hh-suite/);
1. uniprot20 (http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/old-releases/uniprot20_2016_02.tgz).

#### Optionally (for more accurate predictions):
1. PSI-BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/); 
1. UniRef90 (ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz).


## How to run Porter 5
```
# For fast and accurate predictions (exploiting HHblits only)
$ python3 Porter5/Porter5.py -i Porter5/example/2FLGA.fasta --cpu 4 --fast

# For very accurate predictions (exploiting both HHblits and PSI-BLAST)
$ python3 Porter5/Porter5.py -i Porter5/example/2FLGA.fasta --cpu 4
```

### How to run Porter 5 on multiple sequences
```
# To split a FASTA file with multiple sequences (Optional)
$ python3 Porter5/split_fasta.py many_sequences.fasta

# To predict all the fasta files in a given directory (Fastas)
$ python3 Porter5/multiple_fasta.py -i Fastas/ --cpu 4 --fast

# To run multiple predictions in parallel (using a total of 8 cores)
$ python3 Porter5/multiple_fasta.py -i Fastas/ --cpu 4 --parallel 2 --fast
```

### Use the docker image
```
# Set-up docker image
$ docker pull mircare/porter5

# set the absolute PATHs for databases and query sequences (stored locally)
$ docker run --name porter5 -v /**PATH_to_uniprot20_2016_02**:/uniprot20 \
-v /**PATH_to_UniRef90_optional**:/uniref90 -v /**PATH_to_fasta_to_predict**:/Porter5/query \
--cap-add IPC_LOCK mircare/porter5 sleep infinity &

# How to run a prediction using 5 cores and HHblits only
$ docker exec porter5 python3 Porter5.py -i query/2FLGA.fasta --cpu 5 --fast
```


## Performances in 3 states on large independent test set
| Method | Q3 per AA | SOV'99 per AA | Q3 per protein | SOV'99 per protein |
| :--- | :---: | :---: | :---: | :---: |
| **Porter 5** | **83.81%** | **80.41%** | **84.32%** | **81.05%** |
| SPIDER 3 | 83.15% | 79.43% | 83.42% | 79.79% |
| **Porter 5 *HHblits only*** | **83.06%** | **79.49%** | **83.68%** | **80.26%** |
| SSpro 5.1 *with templates* | 82.58% | 78.54% | 83.94% | 80.29% |
| PSIPRED 4.01 | 81.88% | 77.36% | 82.48% | 78.22% |
| RaptorX-Property | 81.86% | 78.08% | 82.57% | 78.99% |
| Porter 4 | 81.66% | 78.05% | 82.29% | 78.61% | 
| SSpro 5.1 *ab initio* | 81.17% | 76.87% | 81.10% | 76.92% |
| DeepCNF | 81.04% | 76.74% | 81.16% | 76.99% |

Calculated with http://dna.cs.miami.edu/SOV/.

## Performances in 8 states on large independent test set in 
| Method | Q8 per AA | SOV8'99 per AA | Q8 per protein | SOV8'99 per protein |
| :--- | :---: | :---: | :---: | :---: |
| **Porter 5** | **73.02%** | **69.91%** | **73.92%** | **70.76%** |
| SSpro 5.1 *with templates* | 71.91% | 68.68% | 74.46% | 71.74% |
| **Porter 5 *HHblits only*** | **71.8%** | **68.87%** | **72.83%** | **69.79%** |
| RaptorX-Property | 70.74% | 67.59% | 71.78% | 68.36% |
| DeepCNF | 69.76% | 66.42% | 70.14% | 66.44% |
| SSpro 5.1 *ab initio* | 68.85% | 65.33% | 69.27% | 65.97% |

Calculated with http://dna.cs.miami.edu/SOV/.


## License
This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

Email us at gianluca[dot]pollastri[at]ucd[dot]ie if you wish to use it for purposes not permitted by the CC BY-NC-SA 4.0.

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a>
