# README

VirHostMatcher-Net is a network-based computational tool for predicting virus-host interactions. Current version predicts hosts of given viruses from a database of 62,493 Bacteria and Archaea candidates. VirHostMatcher-Net has two modes: predicting for complete genomes and predicting for short viral contigs.

### Update on 2021/04/10: Please check the message below before usage!
***Important:*** We corrected a critical bug in the software: the coefficients for SV+ and SV- (see Equation 6 in the paper) were mistakenly swapped when the software was published. Note this does not impact the results in the paper: we computed the feature values separately without the software (because the computation/simulation used in the study was too time-consuming for a single run by the software, especially for the blast feature that we evaluated in the study), and directly applied the coefficients to calculate the final score. 

## Citation
Wang *et al.* "A network-based integrated framework for predicting virus–prokaryote interactions" NAR Genomics and Bioinformatics, Volume 2, Issue 2, June 2020, lqaa044, https://doi.org/10.1093/nargab/lqaa044.


### Dependencies

VirHostMatcher-Net requires Python 3.4+ together with the following packages and `BLAST`.

* Python packages
    + [Biopython](http://biopython.org/wiki/Download)
    + [pandas](https://pandas.pydata.org/) 
    + [numpy](https://www.scipy.org/scipylib/download.html)

We recommend to use [Miniconda](https://conda.io/miniconda.html) to install all dependencies. After installing Miniconda, simply run
```
conda install numpy pandas Biopython 
conda install -c bioconda blast    # skip this line if you've already installed BLAST
``` 

Alternatively, users may want to install Python3.X and dependencies manually. 

## Installation
In addition to dependencies, VirHostMatcher-Net requires to build local modules and download database.

### Building local dependent modules
##### Linux: 
```
git clone https://github.com/WeiliWw/VirHostMatcher-Net.git 
cd VirHostMatcher-Net
CC=g++ python setup.py install --install-platlib=./src/

## Optional
# Include the VirHostmatcher-Net directory in your $PATH
# The main script is executable
export PATH=/path/to/VirHostMatcher-Net/:$PATH
VirHostMatcher.py -h
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
(In MacOS, use command `curl -C -` instead to download the data.)

##### Complete genome mode alone
At the directory of VirHostMatcher-Net, run
```
wget -c https://virhostmatcher-net.s3.us-east-2.amazonaws.com/data_VirHostMatcher-Net_complete_genome_mode_alone.tar.gz    
tar xf data_VirHostMatcher-Net_complete_genome_mode_alone.tar.gz
```

##### Complete genome mode and short viral contig mode
At the directory of VirHostMatcher-Net, run
```
wget -c https://virhostmatcher-net.s3.us-east-2.amazonaws.com/data_VirHostMatcher-Net_both_modes.tar.gz    
tar xf data_VirHostMatcher-Net_both_modes.tar.gz
```

> Be advised
>
> The extracted `data` folder takes up 125G of disk space.

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
```
mkdir output2
python VirHostMatcher-Net.py -q ./test/mVCs --short-contig -o output2 -n 3 -t 8 -l genome_list/marine_host_list.txt
```



In both modes, VirHostMatcher-Net outputs a prediction file for each query virus to the specified directory. A prediction file is in .csv format where each row represents one candidate host with detailed taxanomic information, a prediction score, values (*_val) and percentiles (*_pct) of each feature. The feature percentile of a virus-host pair is defined as the percentile of this feature score among all scores between that virus and all the candidate hosts. A very high percentile (i.e. >95%) suggests significance of the feature on contributing to the prediction. In the output, the percentile of SV<sub>-</sub>, with a negative coefficient, is reversed to keep consistent with other feature percentiles. Tables of feature values are stored in a subdirectory `feature_values` under the output directory.

Users can use a subset of candidate hosts for prediction by the option `-l` to specify a list of NCBI genome names. Two example lists of 4034 marine hosts and 9097 human-associated hosts can be found in the directory `genome_list`.  

### Training and validation data
Training and validation data are shared in [Google Drive](https://drive.google.com/drive/folders/1ilhe-xPQa89jZL8C33NgHNcymz-hkTEC?usp=sharing).

### Bug reports
Please open a Github issue or contact Weili Wang weiliw@usc.edu

-----------------------------------------------------------------------------------------------
## Copyright and License Information
Copyright (C) 2019 University of Southern California, Weili Wang and Fengzhu Sun

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
