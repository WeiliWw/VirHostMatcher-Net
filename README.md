# README

VHiNet is a Python tool for predicting virus-host interactions. Current version predicts hosts of given viruses from a database of 31,986 bacteria candidates. VHiNet has two modes: predicting for complete genomes and predicting for short viral contigs.

### Dependencies
* [Biopython](http://biopython.org/wiki/Download)
* [pandas](https://pandas.pydata.org/) 
* [numpy](https://www.scipy.org/scipylib/download.html)
* [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279671/) 

### Installation
VHiNet requires Python 3.4+. To begin with, install dependencies and make sure the command `blastn` from `BLAST` is added to a directory in `$PATH`.

#### Building local dependent modules
##### Linux: 
```
git ?????
cd VhiNet
python setup.py build_ext 
python setup.py install --install-platlib=.
```
##### MacOS
```
brew install gcc@6
export CC=gcc-6
export CXX=g++-6

git ?????
cd VHiNet
python setup.py build_ext 
python setup.py install --install-platlib=.
```

### Data preparation
The prediction model of VHiNet depends on a large amount of data: BLAST index files of all bacteria and their CRISPRs, WIsH models(short viral contig mode) and hash files for calculating s~2~^*^, etc.

#### Downloading
##### Complete genome mode alone
At the directory of VHiNet, run
```
wget http://www-rcf.usc.edu/~weiliw/vhinet/data_vhinet.tar.gz    
tar xf data_vhinet.tar.gz
```
In MacOS, use `curl` instead to download the data.
##### Complete genome mode and short viral contig mode
At the directory of VHiNet, run
```
wget http://www-rcf.usc.edu/~weiliw/vhinet/data_vhinet_short_viral_contig.tar.gz    
tar xf data_vhinet_short_viral_contig.tar.gz 
```
In MacOS, use `curl` instead to download the data.

##### Required format of query sequences
VHiNet accepts files in FASTA format(i.e., files with .fasta or .fa suffix.)


### Usage 
    VHiNet.py [-h] -q QUERY_VIRUS_DIR [-t NUM_THREADS] [--short-contig] -o
                 OUTPUT_DIR [-n topN]
#### Options
      -h, --help          show this help message and exit
      -q QUERY_VIRUS_DIR  Directory containing query virus genomes with .fasta or
                          .fa suffix
      -t NUM_THREADS      Number of threads (CPUs) to use in the BLAST search.
                          Default = 1
      --short-contig      Predict hosts for short viral contigs. WIsH model files
                          are required in this mode
      -o OUTPUT_DIR       Output directory
      -n topN             Number of top predictions written to the output files.
                          All predictions will be output if there is a tie at the
                          highest score. Default = 1


### examples
#### Predict hosts of virus genomes
```
mkdir tmp
python VHiNet.py -q test_query -t 2 -o tmp -n 5
```

#### Predict hosts of viral contigs
```
mkdir tmp2
python VHiNet.py -q test_query2 -t 2 --short-contig -o tmp2 -n 5
```

