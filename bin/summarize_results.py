#!/usr/bin/env python3

###############################################################################
# Ligase Fidelity Profiling
# Copyright (C) 2022 New England Biolabs, Inc.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
###############################################################################

import argparse
import re
import itertools
import pysam
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('bam1')
parser.add_argument('bam2')
parser.add_argument('-l','--left-bc', default='TTG([ACGT]{6})CGT')
parser.add_argument('-o','--overhang', default='TCC([ACGT]{4})GGA')
parser.add_argument('-r','--right-bc', default='ACG([ACGT]{6})CAA')
parser.add_argument('-n','--num-passes', type=int, default=3)
parser.add_argument('--output-fragments',  default='01_fragments.csv')
parser.add_argument('--output-overhangs',  default='02_overhangs.csv')
parser.add_argument('--output-barcodes',   default='03_barcodes.csv')
parser.add_argument('--output-barcodes-c', default='04_barcodes-c.csv')
parser.add_argument('--output-barcodes-p', default='05_barcodes-p.csv')
parser.add_argument('--output-matrix',     default='06_matrix.csv')
parser.add_argument('--output-fidelity',   default='07_fidelity.csv')
parser.add_argument('--output-mismatch-e', default='08_mismatch-e.csv')
parser.add_argument('--output-mismatch-m', default='09_mismatch-m.csv')

args = parser.parse_args()


def revcomp(seq):
    bases = {
        'A':'T',
        'C':'G',
        'G':'C',
        'T':'A',
        'a':'t',
        'c':'g',
        'g':'c',
        't':'a',
    }

    newseq = ''.join([bases[x] for x in seq])

    return(newseq[::-1])


def count_mismatch(s1,s2):
    count = 0

    for b1,b2 in zip(s1,s2):
        if b1 != b2:
            count += 1

    return(count)


def read_bam(bam):
    # compile regular expressions
    p1 = re.compile(args.left_bc)
    p2 = re.compile(args.overhang)
    p3 = re.compile(args.right_bc)

    bamfile = pysam.AlignmentFile(bam,check_sq=False)

    reads = []

    for read in bamfile:
        ### PacBio reads are reverse complement of an actual substrate
        query_sequence = revcomp(read.query_sequence)

        ### compile search patterns for three parts
        m1 = p1.search(query_sequence)
        m2 = p2.search(query_sequence)
        m3 = p3.search(query_sequence)

        ### all three parts must be present
        if (not m1) or (not m2) or (not m3):
            continue

        np = read.get_tag('np')

        # extract matched patterns
        left_bc  = m1.group(1)
        overhang = m2.group(1)
        right_bc = m3.group(1)

        ### accumulate extracted info
        reads.append([read.qname,np,left_bc,overhang,right_bc])
    
    ### convert to DataFrame
    df = pd.DataFrame(reads, columns = ['qname','np','left_bc','overhang','right_bc'])

    return(df)


def create_fragments_table(bam1, bam2):
    ### strand1 data
    df1 = read_bam(bam1)

    ## strand2 data
    df2 = read_bam(bam2)

    ### merge data from two strands
    df = df1.merge(df2, on='qname', how='inner', suffixes=('1', '2'))

    ### make sure that barcodes match between strands
    df['mm1'] = df.apply(lambda row: count_mismatch(row['left_bc1'],revcomp(row['right_bc2'])), axis=1)
    df['mm2'] = df.apply(lambda row: count_mismatch(row['right_bc1'],revcomp(row['left_bc2'])), axis=1)

    ### discard any reads with mismatch in barcodes
    df = df[(df['mm1'] == 0) & (df['mm2'] == 0)]
    df = df.drop(columns=['mm1','mm2'])

    ### filter based on numer of passes
    df = df[(df['np1'] >= args.num_passes) & (df['np2'] >= args.num_passes)]

    ### count number of mismatches for each overhang pair
    df['overhang_mismatch'] = df.apply(lambda row: count_mismatch(row['overhang1'],revcomp(row['overhang2'])), axis=1)

    return(df)


def summarize_overhangs(df):
    data = {}

    # go through all detected overhangs
    for index, row in df.iterrows():
        o1 = row['overhang1']
        o2 = row['overhang2']

        (o1,o2) = sorted([o1,o2])

        # count unique overhang pairs
        if (o1,o2) not in data:
            data[(o1,o2)] = 0
        
        data[(o1,o2)] += 1

    matrix = []

    # sort by frequency in descending order
    for (o1,o2) in sorted(data, key = lambda x: data[x], reverse = True):
        matrix.append([o1,o2,data[(o1,o2)]])

    # convert to pandas DataFrame
    df = pd.DataFrame(matrix, columns = ['O1','O2','Count'])

    return(df)


def summarize_barcodes(df):
    data = {}

    # go through all barcode sequences
    for index, row in df.iterrows():
        bc1 = row['right_bc1']
        bc2 = row['right_bc2']

        # count unique barcodes
        if bc1 not in data: data[bc1] = 0
        if bc2 not in data: data[bc2] = 0

        data[bc1] += 1
        data[bc2] += 1

    matrix = []

    # sort by frequency in descending order
    for bc in sorted(data, key = lambda x: data[x], reverse = True):
        matrix.append([bc,data[bc]])

    # convert to pandas DataFrame
    df = pd.DataFrame(matrix, columns = ['Barcode','Count'])

    return(df)


def create_overhang_matrix(df):
    ### determine overhang length
    olen = len(df['overhang1'].values[:1][0])

    ### generate all overhangs
    list1 = [''.join(x) for x in itertools.product('ACGT', repeat=olen)]
    list2 = [revcomp(x) for x in list1]

    ### init matrix
    data = {}

    for o1 in list1:
        data[o1] = {}
        for o2 in list2:
            data[o1][o2] = 0
    
    ### populate matrix with actual values
    for index,row in df.iterrows():
        o1 = row['overhang1']
        o2 = row['overhang2']

        data[o1][o2] += 1
        data[o2][o1] += 1

    matrix = []

    # create ovrhang matrix
    for o1 in list2:
        matrix.append([o1] + [data[o1][o2] for o2 in list1])

    # convert to pandas DataFrame
    df = pd.DataFrame(matrix, columns = ['Overhang'] + list1)

    return(df)


def create_barcode_table(df):
    # determine barcode length
    bclen = len(df['right_bc1'].values[:1][0])

    # init table
    data = {}

    cols = ['N%i' % i for i in range(1,bclen+1)]

    for b in 'ACGTN':
        data[b] = [0] * bclen

    # populate table with actual counts
    for index,row in df.iterrows():
        for i,b in enumerate(row['right_bc1']):
            data[b][i] += 1
        for i,b in enumerate(row['right_bc2']):
            data[b][i] += 1

    # convert to pandas DataFrame
    df = pd.DataFrame([
        ['A']+data['A'],
        ['C']+data['C'],
        ['G']+data['G'],
        ['T']+data['T']],
        columns=['Base']+cols)
    
    df['NN'] = df[cols].sum(axis=1)

    # calculate percentages
    df2 = df.copy()
    df2[cols+['NN']] = df2[cols+['NN']] / df2[cols+['NN']].sum()

    return(df,df2)


def create_fidelity_table(matrix):
    # list of overhangs
    olist = matrix.columns.values[1:]

    data = []

    for o1 in olist:
        o2wc = revcomp(o1)

        ### total number of ligation event for a given overhang
        total = np.sum([matrix[matrix['Overhang']==o1][olist]])

        ### number of correct ligation events (WC)
        correct = matrix[matrix['Overhang']==o1][o2wc].values[0]

        ### number of mismatch ligation events
        mismatch = np.sum([matrix[matrix['Overhang']==o1][o2] for o2 in olist if o2 != o2wc])

        ### list of mismatch overhangs
        mm_list = { o2:matrix[matrix['Overhang']==o1][o2].values[0] for o2 in olist if (o2 != o2wc) and (matrix[matrix['Overhang']==o1][o2].values[0] > 0) }
        mm_list_sorted = sorted(mm_list.keys(), key=lambda x: mm_list[x], reverse=True)
        mm_list_len = len(mm_list_sorted)

        mm_list_formatted = None

        if mm_list_len > 0:
            mm_list_formatted = '; '.join(['%s (%.0f%%)' % (o, mm_list[o]/mismatch*100) for o in mm_list_sorted[:5]])

        data.append([o1,total,correct,mismatch,correct/total,mm_list_len,mm_list_formatted])

    ### convert to DataFrame
    df = pd.DataFrame(data, columns=['Overhang','Total','Correct','Mismatch','Fidelity','Number mismatch overhangs','Top 5 mismatch overhangs'])

    return(df)


def create_mismatch_table(df):
    ### determine overhang length
    olen = len(df['O1'].values[:1][0])

    e = {}
    m = {}

    # init
    for b1 in 'ACGT':
        for b2 in 'ACGT':
            e[(b1,b2)] = 0
            m[(b1,b2)] = 0

    # examine each overhang pair
    for index, row in df.iterrows():
        o1 = row['O1']
        o2 = row['O2']
        count = row['Count']

        # overhang pair orientation 1
        top = o1
        bot = o2[::-1]

        # the edge position
        e[(bot[0],top[0])] += count

        # the middle position
        if olen > 2:
            m[(bot[1],top[1])] += count

        # overhang pair orientation 2
        top = o2
        bot = o1[::-1]

        # the edge position
        e[(bot[0],top[0])] += count

        # the middle position
        if olen > 2:
            m[(bot[1],top[1])] += count

    # convert to pandas DataFrame
    arr_e = []
    arr_m = []

    for b1 in 'ACGT':
        for b2 in 'ACGT':
            arr_e.append([b1,b2,e[(b1,b2)]])
            arr_m.append([b1,b2,m[(b1,b2)]])

    df_e = pd.DataFrame(arr_e, columns=['Bottom','Top','Count'])
    df_m = pd.DataFrame(arr_m, columns=['Bottom','Top','Count'])

    return(df_e,df_m)

fragments = create_fragments_table(args.bam1, args.bam2)
fragments.to_csv(args.output_fragments, sep=',', index=False)

overhangs = summarize_overhangs(fragments)
overhangs.to_csv(args.output_overhangs, index=False, sep=',')

barcodes = summarize_barcodes(fragments)
barcodes.to_csv(args.output_barcodes, index=False, sep=',')

matrix = create_overhang_matrix(fragments)
matrix.to_csv(args.output_matrix, index=False, sep=',')

barcode_table_counts, barcode_table_percentages = create_barcode_table(fragments)
barcode_table_counts.to_csv(args.output_barcodes_c, index=False, sep=',')
barcode_table_percentages.to_csv(args.output_barcodes_p, index=False, sep=',')

fidelity_table = create_fidelity_table(matrix)
fidelity_table.to_csv(args.output_fidelity, index=False, sep=',')

mismatch_table_e, mismatch_table_m = create_mismatch_table(overhangs)
mismatch_table_e.to_csv(args.output_mismatch_e, index=False, sep=',')
mismatch_table_m.to_csv(args.output_mismatch_m, index=False, sep=',')
