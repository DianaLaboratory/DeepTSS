#!/usr/bin/env python3.6
import argparse
from src.clustering import *
from src.All_Predict import *

start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, help='Input file | bam or bed', dest='input_file', required=True)
parser.add_argument('-hg', type=str, help='Specify human genome file  | Default hg38', dest='hg38_fa',required=True)
parser.add_argument('-cons', type=str, help='Specify PhyloP ".bw" file | If empty conservation won\'t be used', dest='cons')
parser.add_argument('-g', type=str, help='Specify genome chrom sizes file "e.g. hg38.fa.fai" | Deafult hg38', dest='fai')
parser.add_argument('-Q', type=int, help='Skip alignments with map quality lower than <-Q> | Default: 10', dest='quality')
parser.add_argument('-mapd', type=int, help='Set map distance | Default: 50', dest='map_distance')
parser.add_argument('-t', type=int, help='Set number of threads | Default half of total threads', dest='threads')
parser.add_argument('-tpm', type=float, help='Specify tpm threshold | Default 0.5', dest='tpm')
parser.add_argument('-min', type=int, help='Set minimum cluster size | Default 1', dest='min_clust_size')
parser.add_argument('-max', type=int, help='Set maximum cluster size | Default 1000', dest='max_clust_size')
parser.add_argument('-out', type=str, help='Specify output directory | Optional', dest='out_dir')
parser.add_argument('-dir_name', type=str, help='Specify output dir name | Optional', dest='dir_name')


args = parser.parse_args()

a = my_dictionary()
for arg in vars(args):
    a.add(arg, getattr(args, arg))

a = arg_inputs(a)

if (os.path.splitext(a['input_file'])[1] == ".bam") | (os.path.splitext(a['input_file'])[1] == ".bed"):
    work_dir, representatives = main_clustering(a)
    predictions = main_predict(a, representatives, work_dir)
else:
    sys.exit(f'ERROR: Unknown file type \" {os.path.splitext(a["input_file"])[1]} \"')

print(f'Output TSS file: {predictions}')

print(f'\nTotal time: {round(((time.time() - start_time) / 60), 2)} mins ---')
