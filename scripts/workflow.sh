#!/bin/bash -l

### path to Ligase Fidelity repo
LIGASE_FIDELITY_DIR="/PATH/TO/LIGFI_INSTALLATION_DIRECTORY"

export PATH=$LIGASE_FIDELITY_DIR/bin:$PATH

###############################################################################
### Download example PacBio BAM subreads file                               ###
###############################################################################

# Download PacBio BAM file (m64163_220709_012719.subreads.bam) from NCBI SRA
# (PRJNA894239) to the working directory and rename to movie.subreads.bam

###############################################################################
### Define input data files and output directory                            ###
###############################################################################

# input PacBio BAM subreads file
subreads=$PWD/movie.subreads.bam

# output directory
rundir=$PWD/example

echo ""
echo "subreads   : $subreads"
echo "output dir : $rundir"

###############################################################################
### Cluster and split top/bottom strands                                    ###
###############################################################################

echo ""
echo "Cluster and split reads"

### create output directory for clustered reads
jobdir=$rundir/01-cluster
mkdir -p $jobdir
cd $jobdir

split.py \
    --subread-len 98 \
    --adapter-len 45 \
    --outfile0 subreads.0.txt \
    --outfile1 subreads.1.txt \
    ${subreads}

samtools view --threads 24 -N subreads.0.txt -o subreads.0.bam ${subreads}
samtools view --threads 24 -N subreads.1.txt -o subreads.1.bam ${subreads}

###############################################################################
### Build CCS sequences                                                     ###
###############################################################################

echo ""
echo "Build CCS sequences"

### create output directory for consensus reads
jobdir=$rundir/02-ccs
mkdir -p $jobdir
cd $jobdir

ccs \
    --num-threads=24 \
    --min-passes=3 \
    $rundir/01-cluster/subreads.0.bam subreads_ccs.0.bam

ccs \
    --num-threads=24 \
    --min-passes=3 \
    $rundir/01-cluster/subreads.1.bam subreads_ccs.1.bam

samtools index $rundir/02-ccs/subreads_ccs.0.bam
samtools index $rundir/02-ccs/subreads_ccs.1.bam

###############################################################################
### Summary tables                                                          ###
###############################################################################

echo ""
echo "Tabulate results"

### create output directory for summary tables
jobdir=$rundir/03-summary
mkdir -p $jobdir
cd $jobdir

summarize_results.py \
    --left-bc 'TTG([ACGT]{6})CGT' \
    --overhang 'TCC([ACGT]{4})GGA' \
    --right-bc 'ACG([ACGT]{6})CAA' \
    --num-passes 3 \
    $rundir/02-ccs/subreads_ccs.{0,1}.bam

plot_data.py $rundir/03-summary
