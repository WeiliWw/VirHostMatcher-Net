# README

VirHostMatcher-Net is a network-based computational tool for predicting virus-host interactions. Current version predicts hosts of given viruses from a database of 31,986 bacteria candidates. VirHostMatcher-Net has two modes: predicting for complete genomes and predicting for short viral contigs.

### Dependencies

VirHostMatcher-Net requires Python 3.4+ together with the following packages and `BLAST`.

* Python packages
    + [Biopython](http://biopython.org/wiki/Download)
    + [pandas](https://pandas.pydata.org/) 
    + [numpy](https://www.scipy.org/scipylib/download.html)
* [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK52640/) 

We recommend to use [Miniconda](https://conda.io/miniconda.html) to install all dependencies. After installing Miniconda, simply run
```
conda install numpy pandas Biopython pytables
conda install -c bioconda blast
``` 

Alternatively, users may want to install Python3.X and dependencies manually. In this case, please make sure Python3.X is the default Python version and the command `blastn` from `BLAST` is added to a directory in `$PATH`. To check directories in `$PATH`, run `echo $PATH`.

## Installation
In addition to dependencies, VirHostMatcher-Net requires to build local modules and download data needed.

### Building local dependent modules
##### Linux: 
```
git clone https://github.com/WeiliWw/VirHostMatcher-Net.git 
cd VirHostMatcher-Net
CC=g++ python setup.py install --install-platlib=./src/
```
##### MacOS
```
git clone https://github.com/WeiliWw/VirHostMatcher-Net.git
cd VirHostMatcher-Net
MACOSX_DEPLOYMENT_TARGET=10.9 CC=g++ python setup.py install --install-platlib=./src/
```

### Data preparation
The prediction model of VirHostMatcher-Net depends on a large amount of data: BLAST index files of all bacteria and their CRISPRs, WIsH models(short viral contig mode) and hash files for calculating s<sup>*</sup><sub>2</sub>, etc.

#### Downloading
(In MacOS, use command `curl -C` instead to download the data.)

##### Complete genome mode alone
At the directory of VirHostMatcher-Net, run
```
wget -c http://www-rcf.usc.edu/~weiliw/VirHostMatcher-Net/data_VirHostMatcher-Net_complete_genome_mode_alone.tar.gz    
tar xf data_VirHostMatcher-Net_complete_genome_mode_alone.tar.gz
```

##### Complete genome mode and short viral contig mode
At the directory of VirHostMatcher-Net, run
```
wget -c http://www-rcf.usc.edu/~weiliw/VirHostMatcher-Net/data_VirHostMatcher-Net_both_modes.tar.gz    
tar xf data_VirHostMatcher-Net_both_modes.tar.gz
```

#### Required format of query sequences
VirHostMatcher-Net accepts files in FASTA format.


## Usage 
    python VirHostMatcher-Net.py [-h] -q QUERY_VIRUS_DIR [-t NUM_THREADS] [--short-contig] -o
                 OUTPUT_DIR [-n topN] [-i INTERMEDIATE_DIR] [-l GENOME_LIST]
#### Options
      -h, --help          show this help message and exit
      -q QUERY_VIRUS_DIR  Directory containing query virus genomes with .fasta or
                          .fa suffix
      -t NUM_THREADS      Number of threads to use.  Default = 1
      --short-contig      Predict hosts for short viral contigs. WIsH model files
                          are required in this mode
      -o OUTPUT_DIR       Output directory
      -n topN             Number of top predictions written to the output files.
                          All predictions will be output if there is a tie in 
                          score. Default = 1
      -i INTERMEDIATE_DIR  Directory storing intermediate result. Default =
                          ./intermediate_res                   
      -l GENOME_LIST       Location of the file containing host NCBI genome names of
                           interest

### Examples

#### Predict hosts of virus genomes
```
mkdir output
python VirHostMatcher-Net.py -q ./test/VGs -o output -i tmp -n 3 -t 8
```

#### Predict hosts of viral contigs
(Including 418 contigs, host prediction for this test set may take about one hour.)
```
mkdir output2
python VirHostMatcher-Net.py -q ./test/mVCs --short-contig -o output2 -n 3 -t 8 -l genome_list/hmp359.txt
```



In both modes, VirHostMatcher-Net outputs a prediction file for each query virus to the specified directory. A prediction file is in .csv format where each row represents one candidate host with detailed taxanomic information, a prediction score, values (*_val) and percentiles (*_pct) of each feature. The feature percentile of a virus-host pair is defined as the percentile of this feature score among all scores between that virus and all the candidate hosts. A very high percentile (i.e. >95%) suggests significance of the feature on contributing to the prediction. In the output, the percentile of SV<sub>-</sub>, with a negative coefficient, is reversed to keep consistent with other feature percentiles. Tables of feature values are stored in a subdirectory `feature_values` under the output directory.

Users can use a subset of candidate hosts for prediction by the option `-l` to specify the NCBI genome names. Two example lists of 3529 marine hosts and 359 HMP hosts can be found in the directory `genome_list`. The current entire host database can be found in the table `data_info/hostData.csv`.

Other related information about the training data and benchmark test data can also be found in the directory `data_info`. 

### Bug reports
Please open a Github issue or contact Weili Wang weiliw@usc.edu


