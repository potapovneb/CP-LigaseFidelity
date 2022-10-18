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
import pysam

parser = argparse.ArgumentParser()

parser.add_argument(
    'bam',
    help='Unaligned PacBio subreads BAM file'
)
parser.add_argument(
    '--adapter-len',
    type=int,
    default=45,
    help='Expected SMRTbell adapter length'
)
parser.add_argument(
    '--subread-len',
    type=int,
    default=98,
    help='Expected insert length'
)
parser.add_argument(
    '--smin',
    type=float,
    default=0.75,
    help='Maximum allowed subread length'
)
parser.add_argument(
    '--smax',
    type=float,
    default=1.25,
    help='Minimum allowed subread length'
)
parser.add_argument(
    '--amin',
    type=float,
    default=0.75,
    help='Maximum allowed adapter length'
)
parser.add_argument(
    '--amax',
    type=float,
    default=1.25,
    help='Minimum allowed adapter length'
)
parser.add_argument(
    '--outfile0',
    default='subreads.0.txt'
)
parser.add_argument(
    '--outfile1',
    default='subreads.1.txt'
)

args = parser.parse_args()

# Min/max subread length
smin = args.subread_len * args.smin
smax = args.subread_len * args.smax

# Min/max adapter length
amin = args.adapter_len * args.amin
amax = args.adapter_len * args.amax


def extract_continuous_stretch(subreads, outfile0, outfile1, np=3):
    # at least NP passes are required per strand
    if len(subreads) < np * 2:
        return

    # make sure subreads sorted by subread_start
    subreads = sorted(subreads, key = lambda x: x['subread_start'])

    # init contigs
    contig_id = 0
    contig_size = {}

    subreads[0]['contig'] = contig_id
    contig_size[contig_id] = 1

    for i in range(1,len(subreads)):
        # calculate adapter and subread length
        adapter_len = subreads[i]['subread_start'] - subreads[i-1]['subread_end']
        subread_len = subreads[i]['subread_end'] - subreads[i]['subread_start']

        # if either adapter or subread do not have expected length, increment contig ID
        if adapter_len < amin or adapter_len > amax or subread_len < smin or subread_len > smax:
            contig_id += 1

        # assign contig ID to a given subread
        subreads[i]['contig'] = contig_id

        # count number of subreads per contig
        if contig_id not in contig_size:
            contig_size[contig_id] = 0

        contig_size[contig_id] += 1

    # find the longest contig
    contig_size_sorted = sorted(contig_size, key = lambda x: contig_size[x], reverse=True)

    # at least NP passes are required per strand
    if contig_size[contig_size_sorted[0]] < np * 2:
        return

    # split strands in the longest contig
    count = 0

    for subread in subreads:
        if subread['contig'] == contig_size_sorted[0]:
            if count % 2 == 0:
                print(subread['qname'], file=outfile0)
            else:
                print(subread['qname'], file=outfile1)
            count += 1
    
    return


subreads = []
prev_zmw = None

outfile0 = open(args.outfile0, 'wt')
outfile1 = open(args.outfile1, 'wt')

bamfile = pysam.AlignmentFile(args.bam, check_sq=False)

for read in bamfile:
    movie, zmw, subread = read.qname.split('/')

    subread_start, subread_end = subread.split('_')
    subread_start = int(subread_start)
    subread_end = int(subread_end)

    if prev_zmw == None:
        prev_zmw = zmw

    if zmw != prev_zmw:
        # process subread for a given ZMW
        extract_continuous_stretch(subreads, outfile0, outfile1)

        # reset subreads
        subreads = [{'qname':read.qname, 'subread_start':subread_start, 'subread_end':subread_end, 'contig':None}]

        # reset ZMW
        prev_zmw = zmw
    else:
        # accumulate subreads for a given ZMW
        subreads.append({'qname':read.qname, 'subread_start':subread_start, 'subread_end':subread_end, 'contig':None})

# process the last chunk
extract_continuous_stretch(subreads, outfile0, outfile1)

outfile0.close()
outfile1.close()
