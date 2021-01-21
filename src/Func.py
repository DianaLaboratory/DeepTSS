import os
import sys
import pandas as pd
import time
import multiprocessing

class my_dictionary(dict):
    def __init__(self):
        self = dict()

    def add(self, key, value):
        self[key] = value

def chromosomes():
    chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
       'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    return(chr)

def create_dir(path):
	if not os.path.isdir(path):
		os.system(f'mkdir {path}')

def time_date_id():
    import time
    cur_time = time.asctime()
    cur_time = cur_time.replace(":", "-")
    cur_time = cur_time.replace(" ", "_")
    return(cur_time)

def fix_clusters_represent(clust, representatives, tpm_threshold, ttl_reads):
    df1 = pd.DataFrame.from_records(clust)
    df1 = df1.rename(columns={0 : 'chr', 1:'start',2:'stop',3:'id',4:'reads',5:'strand'})
    df2 = pd.DataFrame.from_records(representatives)

    clusters = pd.concat([df1, df2], axis=1)

    clusters.sort_values(by=['chr', 'start', 'stop'])

    # Scale for cluster lenght
    # clusters['tmp'] = clusters['reads'] / (clusters['stop'] - clusters['start'])
    # Don't scale for cluster length
    clusters['tmp'] = clusters['reads']

    # Count total reads mapped to clusters
    # total_reads = clusters['tmp'].sum()

    # Scale per million
    total_reads = ttl_reads / 1_000_000

    # Final tpm
    clusters['reads'] = clusters['tmp'] / total_reads
    clusters = clusters.rename(columns={'reads' : 'tpm'})

    # Filter tpm >= tpm_threshold
    clusters = clusters[clusters['tpm'] >= tpm_threshold]

    # Round tpm %f.2
    clusters['tpm'] = clusters['tpm'].round(2)

    # Fix id for clusters Dataframe
    clusters['id'] = [f'CTSS_{x}' for x in range(1, len(clusters)+1)]

    represent = pd.DataFrame()
    represent['chr'] = clusters[0]
    represent['start'] = clusters[1]
    represent['stop'] = clusters[1].astype(int) + 1
    represent['id'] = clusters['id']
    represent['tpm'] = clusters['tpm']
    represent['strand'] = clusters[5]

    # Drop columns
    clusters = clusters.drop(columns=['tmp', 0, 1, 2, 3, 4, 5])

    return(clusters, represent)

def arg_inputs(dict):
    if dict['quality'] is None:
        dict['quality'] = 10
    if dict['map_distance'] is None:
        dict['map_distance'] = 50
    if dict['threads'] is None:
        dict['threads'] = int(multiprocessing.cpu_count()/2)
    elif dict['threads'] <= 0:
        sys.exit(f"\nNumber of threads can't be <= 0.\n")

    if dict['fai'] is None:
        dict['fai'] = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), 'files', 'genomes', 'hg38.chrom.sizes')
    elif dict['fai'] == 'hg19':
        dict['fai'] = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), 'files', 'genomes', 'hg19.chrom.sizes')
    elif dict['fai'] == 'hg18':
        dict['fai'] = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), 'files', 'genomes', 'hg18.chrom.sizes')
    elif dict['fai'] == 'hg38':
        dict['fai'] = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), 'files', 'genomes', 'hg38.chrom.sizes')

    if dict['tpm'] is None:
        dict['tpm'] = 1
    elif dict['tpm'] < 0:
        sys.exit(f"\nTpm parameter must be a positive number.\n")

    if dict['min_clust_size'] is None:
        dict['min_clust_size'] = 1
    if dict['max_clust_size'] is None:
        dict['max_clust_size'] = 1000
    # if dict['min_clust_tags'] is None:
    #     dict['min_clust_tags'] = 1

    if dict['out_dir'] is None:
        terminal_path = os.path.abspath(os.getcwd())
        dict['out_dir'] = terminal_path
    elif dict['out_dir'][-1] == '/':
        dict['out_dir'] = dict['out_dir'][:-1]
    # Unique dir with time & date
    if dict['dir_name'] is None:
        dict['dir_name'] = f'DeepTSS_OUT_{time_date_id()}'
    return(dict)