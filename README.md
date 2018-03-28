# Porter5

The web server of Porter5 is available at http://distilldeep.ucd.ie/porter/.
Datasets used for training and testing purposes available at http://distilldeep.ucd.ie/porter/data/.

## Reference
Porter 5: state-of-the-art ab initio prediction of protein secondary structure in 3 and 8 classes
Mirko Torrisi, Manaz Kaleel, Gianluca Pollastri

bioRxiv 289033; doi: https://doi.org/10.1101/289033

## Setup
>$ git clone git@github.com:mircare/Porter5.git
<br/>$ python3 Porter5/Porter5.py -i Porter5/example/2FLGA.fasta


or
>$ wget http://distilldeep.ucd.ie/porter/data/Porter5.tgz<br/>$ tar zxvf Porter5.tgz<br/>$ python3 Porter5/Porter5.py -i Porter5/example/2FLGA.fasta

### Requirements
1. Python3 (https://www.python.org/downloads/);
1. HHblits (https://github.com/soedinglab/hh-suite/);
1. uniprot20 (http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/old-releases/uniprot20_2016_02.tgz).

Optionally:
1. PSI-BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/); 
1. UniRef90 (ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz).


### License

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.
