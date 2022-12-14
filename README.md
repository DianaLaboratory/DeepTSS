## DIANA-DeepTSS - Deep learning computational framework for TSS identification

### Prerequisites
- Unix system

## Installation

##### Install Anaconda python:
Download and install anaconda python [(installation guide)](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)

##### Install conda enviroment using the DeepTSS.yml file provided.
```
conda env create -f DeepTSS.yml
conda activate DeepTSS
```
### 
## Manual installation

##### Create and activate a conda environment with python 3.7
```
conda create -n py37 python=3.7
conda activate py37
```
##### Install python libraries
```
pip intall pandas
pip install numpy
pip install sklearn
pip install pyBigWig
pip install tensorflow-gpu
pip install keras
conda install -c bioconda bedtools
conda install -c bioconda samtools
```

## Download input files

#### PhyloP bigWig file can be downloaded from UCSC 
```
http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw (for -cons tag)
```
 #### Download fasta file from UCSC
```
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz (for -hg tag)
```

## Running the algorithm
```
Examle usage of tool:
python main.py -i H9.bam -hg hg38.fa -cons hg38.phyloP100way.bw

  -h, --help            show this help message and exit
  -i INPUT_FILE         Input file                                           | bam alignments or bed TSS representatives file | Mandatory
  -hg HG38_FA           Specify human genome file (requires fasta file)      | Mandatory
  -cons CONS            Specify PhyloP bigWig file (.bw)                     | If empty conservation branch won't be used
  -g FAI                Specify genome chrom sizes (hg38, hg19, or hg18)     | Deafult hg38
  -Q QUALITY            Skip alignments with map quality lower than <-Q>     | Default: 10
  -mapd MAP_DISTANCE    Set map distance                                     | Default: 50
  -t THREADS            Set number of threads                                | Default half of total threads
  -tpm TPM              Specify tpm threshold                                | Default 1
  -min MIN_CLUST_SIZE   Set minimum cluster size                             | Default 1
  -max MAX_CLUST_SIZE   Set maximum cluster size                             | Default 1000
  -out OUT_DIR          Specify output directory                             | Optional
  -dir_name DIR_NAME    Specify output dir name                              | Optional
```

### Example use

##### Without conservation branch
* `python main.py -i <bam or bed file> -hg <human genome file> `

##### Using conservation branch
* `python main.py -i <bam or bed file> -hg <human genome file> -cons <PhyloP.bw> `

### Mandatory arguments

* `-i <specify bam or bed file>`
   Input bam is a CAGE-seq file.
   Clusters and representatives will be calculated based on bam input.
   If a bed file is provided (instead of bam), the algorithm assumes it contains precalculated CAGE tag cluster representatives.
      
* `-hg <specify genome>`
   User must provide fasta file of desired genome. e.g. hg38.fa

* `-cons <specify PhyloP bigWig file>`
   Specify PhyloP.bw file for evolutionary conservation. If empty CNN conservation branch will not be used.

### Optional arguments

* ` -g <specify genome> `

   Choose between hg38, hg19 or hg18 or provide external chrom.zises file.
   ###### Examples:
   ``` 
   python main.py -i <bam or bed file> -hg <human genome file> -cons <PhyloP.bw>  -g hg19
   python main.py -i <bam or bed file> -hg <human genome file> -cons <PhyloP.bw>  -g /user/path/hg19.chrom.sizes
   ```
<br>

* ` -Q <INT> `

   Reads having quality below <-Q> will be discarded
   
* `-mapd <INT>`

   Maximum required distance between discovered CAGE CTSSs to be consider as a single peak. Default = 50.
    
* `-t <INT>`

    Maximum threads to be used. Default half of system's availiable.

* `-tpm <FLOAT>`

   Minimum tpm value of accepted clusters (>=tpm). Default = 1.
   
* `-min <INT>`

   Minimum cluster size (bp). Default = 1.
    
* `-max <INT>`

   Maximum cluster size (bp). Default = 1000.
   
* `-out <STR>`

   Specify output directory. If empty an output directory will be created in current terminal path.
   ###### Example:
  `python main.py -i <bam or bed file> -hg <human genome file> -cons <PhyloP.bw>  -out /user/mydir`

* `-dir_name <STR>`

   Specify output directory name. If empty an output directory will be created with unique name (based on date & time)
   
   ###### Example:
  `python main.py -i <bam or bed file> -hg <human genome file> -cons <PhyloP.bw>  -dir_name output_dir_name`

### Output explained
DeepTSS returns 3 bed files.
- Clusters file          (clusters around cage peaks)
- Representatives TSS    (Representative TSS for each clusters)
- Scored Representatives (Scored TSS representatives based on the CNN model)



### Authors

Dimitris Grigoriadis, Nikos Perdikopanis, Georgios K. Georgakilas and Artemis Hatzigeorgiou.


### Please cite
Grigoriadis, D., Perdikopanis, N., Georgakilas, G.K. et al. DeepTSS: multi-branch convolutional neural network for transcription start site identification from CAGE data. BMC Bioinformatics 23 (Suppl 2), 395 (2022). https://doi.org/10.1186/s12859-022-04945-y


### License
MIT LICENSE

#### Please address any problems or comments to: 

jim.grigor@gmail.com

