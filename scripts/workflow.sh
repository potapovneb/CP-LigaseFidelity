#!/bin/bash -l

conda activate cp-ligfi

### path to RADAR-seq repo
# LIGFI_DIR="/PATH/TO/LIGFI_INSTALLATION_DIRECTORY"
LIGFI_DIR="/mnt/home/potapov/projects/160202.LigaseFidelity/manuscript/github/ligase-fidelity/CP-LIGFI"

export PATH=$LIGFI_DIR/bin:$PATH

# ###############################################################################
# ### Download example PacBio BAM subreads file                               ###
# ###############################################################################

# wget https://sra-download.ncbi.nlm.nih.gov/traces/sra70/SRZ/019043/SRR19043774/lib26.bam
# ln -s /mnt/pacbio/smrtlink/data/runs/64163/r64163_20220707_183853/4_D01/m64163_220709_012719.subreads.bam movie.subreads.bam

# ###############################################################################
# ### Define input data files and output directory                            ###
# ###############################################################################

# input PacBio BAM subreads file
subreads=$PWD/movie.subreads.bam

# reference FASTA file
reference=$LIGFI_DIR/references/b4.fasta

# output directory
rundir=$PWD/ligfi_output

echo ""
echo "subreads   : $subreads"
echo "reference  : $reference"
echo "output dir : $rundir"

# ###############################################################################
# ### Demultiplexing (optional)                                               ###
# ###############################################################################

# echo ""
# echo "Demultiplexing reads"

# ### create output directory for demultiplexed reads
# jobdir=$rundir/00-demux
# mkdir -p $jobdir
# cd $jobdir

# ### demultiplex sequencing data
# lima $subreads $barcodes movie.bam \
#     --same \
#     --split-bam \
#     --num-threads 16 \
#     --log-level TRACE \
#     --log-file movie.log

# ###############################################################################
# ### Cluster top/bottom strands                                              ###
# ###############################################################################

# echo ""
# echo "Cluster reads"

### create output directory for clustered reads
jobdir=$rundir/01-cluster
mkdir -p $jobdir
cd $jobdir

# time cluster.pl --insert_length 98 ${subreads} | bzip2 - > clusters.csv.bz2

# ###############################################################################
# ### Split top/bottom strands                                                ###
# ###############################################################################

# echo ""
# echo "Split reads"

# time split.pl ${subreads} clusters.csv.bz2

# ###############################################################################
# ### Build CCS sequences                                                     ###
# ###############################################################################

# time samtools view -Sb subreads.fwd.sam > subreads.fwd.bam

# time ccs \
#     --report-file=subreads_ccs.fwd.csv \
#     --log-file=subreads_ccs.fwd.log \
#     --num-threads=8 \
#     --min-passes=1 \
#     subreads.fwd.bam subreads_ccs.fwd.bam

# time samtools view -Sb subreads.rev.sam > subreads.rev.bam

# time ccs \
#     --report-file=subreads_ccs.rev.csv \
#     --log-file=subreads_ccs.rev.log \
#     --num-threads=8 \
#     --min-passes=1 \
#     subreads.rev.bam subreads_ccs.rev.bam

# time TMPDIR=$PWD pbmm2 align $reference subreads_ccs.fwd.bam aligned_reads.fwd.bam \
#     --sort \
#     --min-concordance-perc 75.0 \
#     --sample \"\" \
#     --num-threads 16 \
#     --log-level TRACE \
#     --log-file aligned_reads.fwd.log

# time TMPDIR=$PWD pbmm2 align $reference subreads_ccs.rev.bam aligned_reads.rev.bam \
#     --sort \
#     --min-concordance-perc 75.0 \
#     --sample \"\" \
#     --num-threads 16 \
#     --log-level TRACE \
#     --log-file aligned_reads.rev.log

# time extract.pl \
#     --region 1-3,3-10/1,10-12 \
#     --region 45-47,47-52/1,52-54 \
#     --region 87-89,89-96/1,96-98 \
#     aligned_reads.fwd.bam $reference >fragments.fwd.csv 2>fragments.fwd.log

# time extract.pl \
#     --region 1-3,3-10/1,10-12 \
#     --region 45-47,47-52/1,52-54 \
#     --region 87-89,89-96/1,96-98 \
#     aligned_reads.rev.bam $reference >fragments.rev.csv 2>fragments.rev.log

# time bam2csv.pl subreads_ccs.fwd.bam zmws.fwd.csv
# time bam2csv.pl subreads_ccs.rev.bam zmws.rev.csv

# time reporter.pl \
#     --etype b4 \
#     --match 5 \
#     --np 5 \
#     --bcout barcodes.csv \
#     --ohout overhangs.csv \
#     .

# mktable_barcode.pl barcodes.csv > table-barcode_base_frequencies.csv
# mktable.pl --size 4 --prefix table overhangs.csv
