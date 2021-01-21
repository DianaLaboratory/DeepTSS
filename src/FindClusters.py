#
# Utils for finding clusters and representatives.
#
def props(line):
    chr = line[0]
    start = int(line[1])
    depth = int(line[3])
    strand = line[4]
    cluster_start = int(line[1])
    cluster_end = int(line[2])
    reads = int(line[3])
    return(chr, start, depth, strand, cluster_start, cluster_end, reads)

def clustering(infile, map_dist, max_cluster_length, min_cluster_length):

    representatives = []
    clusters = []
    reads = 0
    i=0
    with open(infile) as f:
        for idx, line in enumerate(f):
            line = line.rstrip().split('\t')
            if idx == 0:
                # Load new line's values to variables
                chr, start, depth, strand, cluster_start, cluster_end, reads = props(line)
                cur_represent = [line[0], line[1], line[2], f'CTSS_{i}', line[3], line[4]]
            else:
                if (line[0] == chr) & (line[4] == strand) & ((int(line[2]) - start) < map_dist) & ((cluster_end - cluster_start) <= max_cluster_length):
                    cluster_end = int(line[2])
                    reads += int(line[3])
                    start = int(line[1])
                    # Find representative with max depth
                    if depth < int(line[3]):
                        depth = int(line[3])
                        cur_represent = [line[0], line[1], line[2], f'CTSS_{i}', line[3], line[4]]
                else:
                    # Append final cluster and representative
                    if ((int(cluster_end)-int(cluster_start)) >= min_cluster_length):
                        clusters.append([line[0], cluster_start, cluster_end, f'CTSS_{i}', reads, strand])
                        representatives.append(cur_represent)
                    i+=1
                    # Load new line's values to variables
                    chr, start, depth, strand, cluster_start, cluster_end, reads = props(line)
                    cur_represent = [line[0], line[1], line[2], f'CTSS_{i}', line[3], line[4]]
    return(clusters, representatives)
