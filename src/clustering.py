import pandas as pd
import multiprocessing as mp
import time
import os
from src.FindClusters import *
from src.Func import *

def _args_(input_file):
    filename = os.path.basename(input_file)
    name = os.path.splitext(filename)[0]
    return(filename, name)

def findAllRepresentatives(line):
    line = line.split('\t')
    if line[5].rstrip() == '+':
        out = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(line[0], line[1], int(line[1].rstrip()) + 1, line[3], line[4], line[5])
    elif line[5].rstrip() == '-':
        out = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(line[0], int(line[2].rstrip()) - 1, line[2], line[3], line[4], line[5])
    return(out)

def genome_cov(row, current, hg38_fa_fai):
    strand = row[0]
    out = row[1]
    os.system(f"bedtools genomecov -i {current} -bg -strand {strand} -g  {hg38_fa_fai} > {out}")


def main_clustering(dict):
    start_time = time.time()

    input_file = dict['input_file']
    hg38_fa_fai = dict['fai']
    quality = dict['quality']
    map_distance = dict['map_distance']
    threads = dict['threads']
    tpm_threshold = dict['tpm']
    min_clust_size = dict['min_clust_size']
    max_clust_size = dict['max_clust_size']
    work_dir = dict['out_dir']
    result_directory_id = dict['dir_name']

    # input_file = '/mnt/raid0/dgrigor/projects/TSS_CNN/Clustering/test_H9/H9.bam'
    # hg38_fa_fai = '/mnt/raid0/dgrigor/.GitHub/TSS_CNN/DeepTSS_gitLab/files/genomes/hg38.chrom.sizes'
    # quality = 10
    # map_distance = 50
    # threads = 16
    # tpm_threshold = 0.5
    # min_clust_size = 1
    # max_clust_size = 1000
    # work_dir = '/mnt/raid0/dgrigor/projects/TSS_CNN/Clustering/test_H9'
    # result_directory_id = 'test_tpm'

    # Set path for work/out dir to new folder
    work_dir = f'{work_dir}/{result_directory_id}'

    # Create work/results dir
    create_dir(work_dir)

    # Save filename with and without extension
    [filename, name] = _args_(input_file)

    # If ".bed" file was provided
    if os.path.splitext(input_file)[1] == ".bed":
        return(work_dir, input_file)

    # Get all chr1-22 & chrX, chrY
    chr = chromosomes()

    # Create tmp dir for leftovers
    tmp_dir = f'{work_dir}/{name}_tmpdir'
    create_dir(tmp_dir)

    # Bam to bed & Quality filter
    print(" (1/6) - Applying quality threashold...")
    # out = f'{tmp_dir}/{name}_quality.bed'
    # os.system(f"bedtools bamtobed -i {input_file} | awk -F '\t' -v OFS='\t' '{{ if ($5 > {quality}) print $1,$2,$3,$4,$5,$6 }}' > {out}")

    # Forward strand
    out_f = f'{tmp_dir}/{name}_quality_f.bed'
    os.system(f"samtools view -F 16 -@ {threads} -q {quality} {input_file}  | awk -F '\t' -v OFS='\t' '{{print $3, $4, $4+length($10)-1, $1, $5, \"+\"}}' > {out_f}")
    # Reverse strand
    out_r = f'{tmp_dir}/{name}_quality_r.bed'
    os.system(f"samtools view -f 16 -@ {threads} -q {quality} {input_file} | awk -F '\t' -v OFS='\t' '{{print $3, $4, $4+length($10)-1, $1, $5, \"-\"}}' > {out_r}")

    current = f'{tmp_dir}/{name}_quality.bed'
    os.system(f'cat {out_f} {out_r} > {current}')

    # Calculate ctss
    print(" (2/6) - CTSS...")
    tmp = []
    with open(current, 'r') as infile:
        with mp.Pool(threads) as p:
            tmp.extend( p.map(findAllRepresentatives, [line for line in infile]) )

    # Write ctss to file
    out = f'{tmp_dir}/{name}_quality_ctss.bed'
    with open(out, 'w') as f:
        for item in tmp:
            if item.split("\t")[0] in chr:  # Filter chr 0-22 & X & Y
                f.write("%s" % item)
    del tmp
    current = out

    # Sort bed file
    print(" (3/6) - Sorting bed...")
    out = f'{tmp_dir}/{name}_quality_ctss_sorted.bed'
    os.system(f'sort --parallel={threads} -k1,1 -k2,2n -V -s {current} > {out}')
    current = out

    ################################### Coverage ##############################################
    # Coverage for all candidate tss representatives
    print(" (4/6) - Computing coverage...")
    out_5 = f'{tmp_dir}/{name}_quality_ctss_sorted_coverage_5_tmp.bed'
    out_3 = f'{tmp_dir}/{name}_quality_ctss_sorted_coverage_3_tmp.bed'

    cov_input = [['+', out_5], ['-',out_3]]
    pool = mp.Pool(2)
    pool.starmap(genome_cov, [(row, current, hg38_fa_fai) for row in cov_input])
    pool.close()

    # Merge 5' and 3' coverage
    out = f'{tmp_dir}/{name}_quality_ctss_sorted_coverage_5.bed'
    os.system(f'awk -F \"\t\" -v OFS=\"\t\" \'{{print $1, $2, $3, $4, "+"}}\' {out_5} > {out}')
    os.system(f'rm {out_5}')
    # Sort
    out_5 = f'{tmp_dir}/{name}_quality_ctss_sorted_coverage_all_sorted_5.bed'
    os.system(f"sort --parallel={threads} -k1,1 -k2,2n -V -s {out} > {out_5}")

    out = f'{tmp_dir}/{name}_quality_ctss_sorted_coverage_3.bed'
    os.system(f'awk -F \"\t\" -v OFS=\"\t\" \'{{print $1, $2, $3, $4, "-"}}\' {out_3} > {out}')
    os.system(f'rm {out_3}')
    # Sort
    out_3 = f'{tmp_dir}/{name}_quality_ctss_sorted_coverage_all_sorted_3.bed'
    os.system(f"sort --parallel={threads} -k1,1 -k2,2n -V -s {out} > {out_3}")

    # Clustering
    print(" (5/6) - Clustering...")
    clusters_5, representatives_5 = clustering(out_5, map_distance, max_clust_size, min_clust_size)
    clusters_3, representatives_3 = clustering(out_3, map_distance, max_clust_size, min_clust_size)

    # Combine two strands for clusters and representatives
    clust = clusters_5 + clusters_3
    represent = representatives_5 + representatives_3

    if not clust:
        sys.exit(f'No clusters found. Try to reduce minimum cluster size')

    print(" (6/6) - Calculating tpm...")
    total_reads =  sum(1 for line in open(f'{tmp_dir}/{name}_quality.bed'))
    clusters, representatives = fix_clusters_represent(clust, represent, tpm_threshold, total_reads)

    # Save
    out_clusters = f'{work_dir}/{name}_clusters.tpm.bed'
    clusters.to_csv(out_clusters, header=None, index=False, sep="\t")

    out_representatives = f'{work_dir}/{name}_representative.tpm.bed'
    representatives.to_csv(out_representatives, header=None, index=False, sep="\t")

    # os.system(f'rm -R {tmp_dir}')

    print(f'Clustering time: {round(((time.time() - start_time) / 60), 2)} mins ---')

    return(work_dir, out_representatives)