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

# reference FASTA file
reference=$LIGASE_FIDELITY_DIR/references/b4.fasta

# output directory
rundir=$PWD/example/ligase_fidelity_output

echo ""
echo "subreads   : $subreads"
echo "reference  : $reference"
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

echo "cluster.pl"
time cluster.pl --insert_length 98 ${subreads} | bzip2 - > clusters.csv.bz2

echo "split.pl"
time split.pl ${subreads} clusters.csv.bz2

echo "samtools (1)"
time samtools view -Sb subreads.fwd.sam > subreads.fwd.bam

echo "samtools (1)"
time samtools view -Sb subreads.rev.sam > subreads.rev.bam

rm -f subreads.fwd.sam
rm -f subreads.rev.sam

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
    --report-file=subreads_ccs.fwd.csv \
    --log-file=subreads_ccs.fwd.log \
    --num-threads=8 \
    --min-passes=1 \
    $rundir/01-cluster/subreads.fwd.bam subreads_ccs.fwd.bam

samtools index $rundir/01-cluster/subreads.fwd.bam subreads_ccs.fwd.bam

echo "ccs (2)"
time ccs \
    --report-file=subreads_ccs.rev.csv \
    --log-file=subreads_ccs.rev.log \
    --num-threads=8 \
    --min-passes=1 \
    $rundir/01-cluster/subreads.rev.bam subreads_ccs.rev.bam

samtools index $rundir/01-cluster/subreads.rev.bam subreads_ccs.rev.bam

###############################################################################
### Summary tables                                                          ###
###############################################################################

echo ""
echo "Tabulate results"

### create output directory for summary tables
jobdir=$rundir/02-summary
mkdir -p $jobdir
cd $jobdir

echo "summarize_results.py"
time summarize_results.py $rundir/02-ccs/subreads_ccs.{fwd,rev}.bam
