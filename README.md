# README

VHiNet is a network-based computational tool for predicting virus-host interactions. Current version predicts hosts of given viruses from a database of 31,986 bacteria candidates. VHiNet has two modes: predicting for complete genomes and predicting for short viral contigs.

### Dependencies

VHiNet requires Python 3.4+ together with the following packages and `BLAST`.

* Python packages
    + [Biopython](http://biopython.org/wiki/Download)
    + [pandas](https://pandas.pydata.org/) 
    + [numpy](https://www.scipy.org/scipylib/download.html)
* [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK52640/) 

We recommend to use [Miniconda](https://conda.io/miniconda.html) to install all dependencies. After installing Miniconda, simply run
```
conda install numpy pandas Biopython
conda install -c bioconda blast
``` 

Alternatively, users may want to install Python3 and dependencies manually. In this case, please make sure Python3.X is the default Python version and the command `blastn` from `BLAST` is added to a directory in `$PATH`. To check directories in `$PATH`, run `echo $PATH`.

## Installation
In addition to dependencies, VHiNet requires to build local modules and download data needed.

### Building local dependent modules
##### Linux: 
```
git clone https://github.com/WeiliWw/VHiNet.git 
cd VhiNet
CC=gcc python setup.py build_ext 
python setup.py install --install-platlib=./src/
```
##### MacOS
For MacOS, please use [brew](https://brew.sh/) to install `gcc-6` first.
```
brew install gcc@6
git clone https://github.com/WeiliWw/VHiNet.git
cd VHiNet
export CC=gcc-6
export CXX=g++-6
python setup.py build_ext 
python setup.py install --install-platlib=./src/
```

### Data preparation
The prediction model of VHiNet depends on a large amount of data: BLAST index files of all bacteria and their CRISPRs, WIsH models(short viral contig mode) and hash files for calculating s<sup>*</sup><sub>2</sub>, etc.

#### Downloading
##### Complete genome mode alone
At the directory of VHiNet, run
```
wget -c http://www-rcf.usc.edu/~weiliw/vhinet/data_vhinet.tar.gz    
tar xf data_vhinet.tar.gz
```

##### Complete genome mode and short viral contig mode
At the directory of VHiNet, run
```
wget -c http://www-rcf.usc.edu/~weiliw/vhinet/data_vhinet_short_viral_contig.tar.gz    
tar xf data_vhinet_short_viral_contig.tar.gz 
```
In MacOS, use command `curl` instead to download the data.

#### Required format of query sequences
VHiNet accepts files in FASTA format(i.e., files with .fasta or .fa suffix.)


## Usage 
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
                          All predictions will be output if there is a tie in 
                          score. Default = 1

### Examples

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

In both modes, VHiNet outputs a prediction file for each query virus to the specified directory. A prediction file is in .csv format where each row represents one candidate host with detailed taxanomic information and prediction score. Matrices of feature values are stored in a subdirectory `feature_values` under the output directory.
