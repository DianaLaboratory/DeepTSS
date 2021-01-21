from sklearn.preprocessing import StandardScaler
import math
import pandas as pd
import numpy as np
import re
import pyBigWig
import pickle
import multiprocessing as mp
import time
import os
import sys

from keras.models import load_model


class DNArepresent:
    def __init__(self, sequence, chrom, start, stop, strand, conservation_path):
        self.sequence = sequence.upper()
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.strand = strand
        self.conservation_path = conservation_path

    def list_of_zeros(self):
        x = [0]*len(self.sequence)
        y = [0]*len(self.sequence)
        z = [0]*len(self.sequence)
        return(x, y, z)

    def DNA_walk(self):
        x = []
        for i, f in enumerate(self.sequence):
            if i == 0:
                if f == 'C' or f == 'T':
                    x.append(1)
                else:
                    x.append(-1)
            else:
                if f == 'C' or f == 'T':
                    x.append(x[i-1]+1)
                else:
                    x.append(x[i-1]-1)
        return(x)

    def Z_curve(self):
        x,y,z = self.list_of_zeros()

        for i, f in enumerate(self.sequence):
            if f == 'T' or f == 'G':
                x[i] = x[i-1] + 1
            else:
                x[i] = x[i-1] - 1
            if f == 'A' or f == 'C':
                y[i] = y[i-1] + 1
            else:
                y[i] = y[i-1] - 1
            if f == 'A' or f == 'T':
                z[i] = z[i-1] + 1
            else:
                z[i] = z[i-1] - 1
        return(x, y, z)

    def paired_numeric(self):
        x = []
        for f in self.sequence:
            if f == 'A' or f == 'T':
                x.append(1)
            else:
                x.append(-1)
        return(x)

    def tetrahedron(self):
        x,y,z = self.list_of_zeros()

        for i, f in enumerate(self.sequence):
            if f == 'T':
                x[i] = 2*math.sqrt(2)/3
                y[i] = 0
                z[i] = -1/3
            if f == 'C':
                x[i] = -math.sqrt(2)/3
                y[i] = math.sqrt(6)/3
                z[i] = -1/3
            if f == 'G':
                x[i] = -math.sqrt(2)/3
                y[i] = -math.sqrt(6)/3
                z[i] = -1/3
            if f == 'A':
                x[i] = 0
                y[i] = 0
                z[i] = 1
        return(x, y, z)

    @classmethod
    def onehot_conversion_sequence(cls, letter):
        one_hot_map = {
            "A": np.asarray([1, 0, 0, 0],dtype=np.float32), "a": np.asarray([1, 0, 0, 0],dtype=np.float32),
            "C": np.asarray([0, 1, 0, 0],dtype=np.float32), "c": np.asarray([0, 1, 0, 0],dtype=np.float32),
            "G": np.asarray([0, 0, 1, 0],dtype=np.float32), "g": np.asarray([0, 0, 1, 0],dtype=np.float32),
            "T": np.asarray([0, 0, 0, 1],dtype=np.float32), "t": np.asarray([0, 0, 0, 1],dtype=np.float32),
            "N": np.asarray([0, 0, 0, 0],dtype=np.float32), "n": np.asarray([0, 0, 0, 0],dtype=np.float32)}
        return one_hot_map[letter]

    def one_hot_encoder(self):
        tmp = []
        for letter in self.sequence:
            tmp.append(self.onehot_conversion_sequence(letter))
        out = np.vstack(tmp)
        return (out)

    def bendability(self, size=3):
        out = [0, 0, 0]
        for x in range(0, len(self.sequence) - size):
            kmer = self.sequence[x:x + size]
            if 'N' in kmer:
                out.append('0')
            else:
                out.append(VALS_kmers[COLS_kmers.index(kmer)])
        if self.strand == '-':
            out = out[::-1]
        out =  np.vstack(out).astype(float)
        return(out)

    def propellerTwist(self, size=2):
        out = ['-12.6', '-12.6']
        for x in range(0, len(self.sequence) - size):
            kmer = self.sequence[x:x + size]
            if 'N' in kmer:
                out.append('-12.6')
            else:
                out.append(VALS_kmers_proptwst[COLS_kmers_proptwst.index(kmer)])
        if self.strand == '-':
            out = out[::-1]
        out =  np.vstack(out).astype(float)
        return(out)

    def conservation_calc(self):
        if self.conservation_path is not None:
            bw = pyBigWig.open(self.conservation_path)
            out = bw.values(self.chrom, self.start, self.stop)
            out = np.vstack(out)
            if self.strand == '-':
                out = np.flip(out)
            bw.close()
            return(out)
        else:
            return(False)

def Read_file(infile):
    with open(infile) as f:
        lines = f.readlines()
        cols = []
        vals = []
        for line in lines:
            cols.append(line.strip().split("\t")[0])
            vals.append(line.strip().split("\t")[1])
    return(cols, vals)


def read_fasta(path_fasta):
    fasta = pd.read_csv(path_fasta, header=None, sep="\t")
    fasta[['chr', 'strand']] = fasta[0].str.split("(", expand=True)
    fasta['strand'] = fasta['strand'].str[:-1]

    fasta[0] = fasta[0].str[:-3]
    fasta[['chr', 'start']] = fasta[0].str.split(":", expand=True)
    fasta[['start', 'stop']] = fasta['start'].str.split("-", expand=True)
    fasta['sequence'] = fasta[1]
    fasta = fasta.drop([0, 1], axis=1)

    # Cast to int
    fasta['start'] = fasta['start'].astype(int)
    fasta['stop'] = fasta['stop'].astype(int)
    return(fasta)

def fix_coords(bed, sequence_len, path_bed):
    bed['tmp'] = (((bed[1] + bed[2]) / 2) - (sequence_len / 2)).astype(int)
    bed[2] = (((bed[1] + bed[2]) / 2) + (sequence_len / 2)).astype(int)
    bed[1] = bed['tmp']
    bed = bed.drop(['tmp'], axis=1)
    bed = bed.rename(columns={0: 'chr', 1: 'start', 2: 'stop', 3: 'id', 4: 'tpm', 5: 'strand'})

    # Save pos/neg set (to utilize bedtools getfasta)
    bed.to_csv(path_bed, sep="\t", header=None, index=False)
    return(bed)

def features(index, row):
    # Create tss object
    tss = DNArepresent(row['sequence'], row['chr'], row['start'], row['stop'], row['strand'], conservation_path)

    # One-hot encoder sequence
    enc_seq = tss.one_hot_encoder()

    # Calculate DNA representations
    dnawalk = tss.DNA_walk()
    x, y, z = tss.Z_curve()
    prd_Num = tss.paired_numeric()
    r, g, b = tss.tetrahedron()
    # Stack vertically all representations
    dna_represent = np.vstack([dnawalk, x, y, z, prd_Num, r, g, b]).T

    # Bendability
    bend = tss.bendability()

    # Propeller_Twist
    propTwist = tss.propellerTwist()

    if conservation_path is not None:
        # Conservation
        conservation = tss.conservation_calc()
        cur = np.hstack([enc_seq, dna_represent, bend, propTwist, conservation])
    else:
        cur = np.hstack([enc_seq, dna_represent, bend, propTwist])

    return( cur )

def split_features_branch(X, conservation_path):
    # Split features for 3 branches
    X1 = X[:, :, [0, 1, 2, 3]]                # seq x4
    X2 = X[:, :, [4, 5, 6, 7, 8, 9, 10, 11]]  # features x8
    X3 = X[:, :, [12,13]]                     # bendabillity, PropellerTwist x2
    if conservation_path is not None:
        X4 = X[:, :, [14]]                        # conservation x1
        return(X1, X2, X3, X4)
    else:
        return (X1, X2, X3)

def standardize_feature(X, scaler=StandardScaler()):
    X = scaler.fit_transform(X.reshape(-1, X.shape[-1])).reshape(X.shape)
    return(X)

def export_file(result, prediction):
    cur = result.copy()
    cur['pred'] = prediction
    cur['tpm'] = cur['tpm'].round(3)
    cur['id'] = cur[['id', 'tpm']].astype(str).agg('||'.join, axis=1)
    cur = cur.drop(['sequence', 'tpm'], axis=1)
    cols = ['chr', 'start', 'stop', 'id', 'pred', 'strand']
    cur = cur[cols]
    cur['start'] = ((cur['start'].astype(int) + cur['stop'].astype(int))/2).astype(int)
    cur['stop'] = cur['start'].astype(int) + 1
    return(cur)

#------------------------------------------##------------------------------------------#
#------------------------------------------##------------------------------------------#
#------------------------------------------##------------------------------------------#


def main_predict(dict, representatives, work_dir):
    start_time = time.time()

    global ROOT_DIR
    global conservation_path
    global COLS_kmers
    global VALS_kmers
    global COLS_kmers_proptwst
    global VALS_kmers_proptwst

    out_dir = work_dir
    cores = dict['threads']
    conservation_path = dict['cons']
    inputBed_path = representatives
    hg38_fa = dict['hg38_fa']

    sequence_len = 600

###############################################
    # Current script dir
    ROOT_DIR = os.path.abspath(os.path.dirname(sys.argv[0]))

    scalers_path = os.path.join(ROOT_DIR, 'models', 'scalers')
    scalers_path_cons = os.path.join(ROOT_DIR, 'models', 'scalers_cons')

    model_path = os.path.join(ROOT_DIR, 'models', 'model17-10.hdf5')
    model_path_cons = os.path.join(ROOT_DIR, 'models', 'model28-09.hdf5')

    COLS_kmers, VALS_kmers = Read_file(f'{ROOT_DIR}/files/bendability.tsv')
    COLS_kmers_proptwst, VALS_kmers_proptwst = Read_file(f'{ROOT_DIR}/files/Propeller_twist.tsv')


    # Read representatives file (*.bed)
    inputBed = pd.read_csv(inputBed_path, header=None, sep="\t")

    # Fix coordinates start-sequence_len | stop+sequence_len from TSS (in case of coordinate > 1)
    path_bed = f'{out_dir}/{os.path.basename(inputBed_path).split(".")[0]}_CNN.scored_tmp.bed'
    inputBed = fix_coords(inputBed, sequence_len, path_bed)

    # Bedtools getfasta
    path_fasta = f'{os.path.dirname(path_bed)}/labeled.fasta'
    cmd = f'bedtools getfasta -fi {hg38_fa} -bed {path_bed} -s -tab -fo {path_fasta}'
    os.system(cmd)

    # Read the fasta
    fasta = read_fasta(path_fasta)

    # Merge fasta & bed
    result = pd.merge(left=inputBed, right=fasta, how='left', left_on=['chr', 'start', 'stop', 'strand'], right_on=['chr', 'start', 'stop', 'strand'])

    # Drop empty rows
    result = result.dropna()

    print(" (1/3) Calculating features...")
    # Calculate features
    with mp.Pool(cores) as p:
        tmp = np.stack(p.starmap(features, [row for row in result.iterrows()]))

    print(" (2/3) Scale...")
    # Split features/inputs for each branch
    if conservation_path is not None:
        X1, X2, X3, X4 = split_features_branch(tmp, conservation_path)
        X4[np.isnan(X4)] = 0
    else:
        X1, X2, X3 = split_features_branch(tmp, conservation_path)

    # scale
    if conservation_path is not None:
        scaler_X2 = pickle.load(open(f'{scalers_path_cons}/X2.pkl', 'rb'))
        scaler_X3 = pickle.load(open(f'{scalers_path_cons}/X3.pkl', 'rb'))
    else:
        scaler_X2 = pickle.load(open(f'{scalers_path}/X2.pkl', 'rb'))
        scaler_X3 = pickle.load(open(f'{scalers_path}/X3.pkl', 'rb'))

    X2 = standardize_feature(X2, scaler_X2)
    X3 = standardize_feature(X3, scaler_X3)

    print(" (3/3) Predict...")
    # Predict
    if conservation_path is not None:
        model = load_model(model_path_cons)
        prediction = model.predict([X1, X2, X3, X4])
    else:
        model = load_model(model_path)
        prediction = model.predict([X1, X2, X3])

    # Round output
    prediction = prediction.round(2)

    # Include predictions to final file | Create final file
    out = export_file(result, prediction)

    out_name = f'{out_dir}/{os.path.basename(inputBed_path).split(".")[0]}.Scored.bed'
    out.to_csv(out_name, header=None, index=False, sep="\t")

    # Delete leftovers
    cmd = f'rm -r {path_fasta} {path_bed}'
    os.system(cmd)
    print(f'\nPredict time: {round( ((time.time() - start_time)/60),2)} mins ---' )

    return(out_name)