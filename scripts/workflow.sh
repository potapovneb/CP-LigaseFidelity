#!/bin/bash -l

conda activate cp-ligfi

### path to Ligase Fidelity repo
# LIGASE_FIDELITY_DIR="/PATH/TO/LIGFI_INSTALLATION_DIRECTORY"
LIGASE_FIDELITY_DIR="/mnt/home/potapov/projects/160202.LigaseFidelity/manuscript/github/ligase-fidelity/CP-LIGFI"

export PATH=$LIGASE_FIDELITY_DIR/bin:$PATH

###############################################################################
### Download example PacBio BAM subreads file                               ###
###############################################################################

# wget https://sra-download.ncbi.nlm.nih.gov/traces/sra70/SRZ/019043/SRR19043774/lib26.bam
ln -s /mnt/pacbio/smrtlink/data/runs/64163/r64163_20220707_183853/4_D01/m64163_220709_012719.subreads.bam movie.subreads.bam

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

echo "split"
time split.py \
    --subread-len 98 \
    --adapter-len 45 \
    --outfile0 subreads.0.txt \
    --outfile1 subreads.1.txt \
    ${subreads}

echo "samtools (0)"
time samtools view --threads 24 -N subreads.0.txt -o subreads.0.bam ${subreads}

echo "samtools (1)"
time samtools view --threads 24 -N subreads.1.txt -o subreads.1.bam ${subreads}

###############################################################################
### Build CCS sequences                                                     ###
###############################################################################

echo ""
echo "Build CCS sequences"

### create output directory for consensus reads
jobdir=$rundir/02-ccs
mkdir -p $jobdir
cd $jobdir

echo "ccs (1)"
time ccs \
    --num-threads=24 \
    --min-passes=3 \
    $rundir/01-cluster/subreads.0.bam subreads_ccs.0.bam

samtools index $rundir/02-ccs/subreads_ccs.0.bam

echo "ccs (2)"
time ccs \
    --num-threads=24 \
    --min-passes=3 \
    $rundir/01-cluster/subreads.1.bam subreads_ccs.1.bam

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

echo "summarize_results.py"
time summarize_results.py \
    --left-bc 'TTG([ACGT]{6})CGT' \
    --overhang 'TCC([ACGT]{4})GGA' \
    --right-bc 'ACG([ACGT]{6})CAA' \
    --num-passes 3 \
    $rundir/02-ccs/subreads_ccs.{0,1}.bam

###############################################################################
### Plot data (optional)                                                    ###
###############################################################################

echo ""
echo "Plot ligation fidelity data"

mkdir -p $rundir/03-summary/figures
cd $rundir/03-summary/figures
plot_data.py $rundir/03-summary
