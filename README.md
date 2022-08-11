# Ligase Fidelity Profiling

This is an official repository for computational workflow on profiling ligase fidelity using PacBio SMRT sequencing. Please see citing information below.

# Prerequisites

The computational workflow requires a number of tools to be installed and available from the command line in your system.

## PacBio tools

* ```lima``` (optional) - Demultiplexing PacBio data. Available as part of [pbbioconda](https://github.com/PacificBiosciences/pbbioconda).
* ```pbmm2``` - minimap2 SMRT wrapper for PacBio data. Available as part of [pbbioconda](https://github.com/PacificBiosciences/pbbioconda).
* ```ccs``` - Generating circular consensus sequences (CCS) based on PacBio sequencing data. Available as part of [pbbioconda](https://github.com/PacificBiosciences/pbccs).

## Third-party tools

* ```samtools``` - Handling high-throughput sequencing data. Available through [conda](https://anaconda.org/bioconda/samtools) or [GitHub](https://github.com/samtools/samtools).

## Custom PERL scripts

These custom scripts are provided as part of Ligase Fidelity GitHub repository. The scripts must be made available from the command line in your system (for example, by adding the scripts directory to ```$PATH``` environment variable).

* ```cluster.pl``` - Identify top/bottom subreads in PacBio sequencing reads.
* ```split.pl``` - Split top/bottom subreads to separate files.
* ```extract.pl``` - Extract fragments of sequenced amplicons.
* ```bam2csv.pl``` - Extract stats for PacBio sequencing reads.
* ```reporter.pl``` - Output overhang and barcode sequences.
* ```mktable_barcode.pl``` - Compute nucleotide frequencies in barcodes.
* ```mktable.pl``` - Tabulate ligase fidelity results.

# Computational workflow

Processing PacBio sequencing data in the Ligase Fidelity computational workflow proceeds through a series of steps, where output of one step serves as input for the next step. As a result, this workflow produces a set of tables summarizing ligase fidelity data.

The ```scripts/``` directory contains an example workflow (```workflow.sh```) that can be executed from the command line. Before executing the example workflow, make sure that the ```$LIGASE_FIDELITY_DIR``` environment variable in ```workflow.sh``` points to the Ligase Fidelity GitHub repository on your computer.

The computational steps are described below.

## Example PacBio BAM sequencing data

For the purpose of this example, you can download PacBio sequencing data for the b4 substrate sequencing.

```
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra70/SRZ/019043/SRR19043774/lib26.bam
```

For convenience, we define a set of environment variables to define location of the input data files and the output directory. Please update these variables according to your computational environment.

Input PacBio BAM subreads file:
```
subreads=$PWD/movie.subreads.bam
```

Reference FASTA file:
```
reference=$LIGASE_FIDELITY_DIR/references/b4.fasta
```

Output directory:
```
rundir=$PWD/ligase_fidelity_output
```

This repository comes with the reference FASTA files for the 4-base overhang substrate (```b4.fasta```).

## Cluster and split top/bottom strands

In this step, ```cluster.pl``` script detects a continuous stretch of insert sequences separated by adapter sequences (based on expected lengths of insert and adapter sequences) in PacBio long polymerase reads. Then ```split.pl``` script splits detected insert sequences to two separate groups.

```
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
```

## Build CCS sequences

PacBio ```ccs``` tool builds consensus sequences for "top" and "bottom" strands for each sequenced product.

```
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

echo "ccs (2)"
time ccs \
    --report-file=subreads_ccs.rev.csv \
    --log-file=subreads_ccs.rev.log \
    --num-threads=8 \
    --min-passes=1 \
    $rundir/01-cluster/subreads.rev.bam subreads_ccs.rev.bam
```

## Map CCS sequences

Consensus sequences are mapped to the expected reference using PacBio ```pbmm2``` tool.

```
### create output directory for mapped consensus reads
jobdir=$rundir/03-mapping
mkdir -p $jobdir
cd $jobdir

echo "pbmm2 (1)"
time TMPDIR=$PWD pbmm2 align $reference $rundir/02-ccs/subreads_ccs.fwd.bam aligned_reads.fwd.bam \
    --sort \
    --min-concordance-perc 75.0 \
    --sample \"\" \
    --num-threads 16 \
    --log-level TRACE \
    --log-file aligned_reads.fwd.log

echo "pbmm2 (2)"
time TMPDIR=$PWD pbmm2 align $reference $rundir/02-ccs/subreads_ccs.rev.bam aligned_reads.rev.bam \
    --sort \
    --min-concordance-perc 75.0 \
    --sample \"\" \
    --num-threads 16 \
    --log-level TRACE \
    --log-file aligned_reads.rev.log
```

## Fragment CCS sequences

Parts of the top and bottom strands corresponding to overhangs and internal barcode regions are detected using ```extract.pl``` script. See insert structure below for details.

```
### create output directory for fragmented consensus reads
jobdir=$rundir/04-fragments
mkdir -p $jobdir
cd $jobdir

echo "extract.pl (1)"
time extract.pl \
    --region 1-3,3-10/1,10-12 \
    --region 45-47,47-52/1,52-54 \
    --region 87-89,89-96/1,96-98 \
    $rundir/03-mapping/aligned_reads.fwd.bam $reference >fragments.fwd.csv 2>fragments.fwd.log

echo "extract.pl (2)"
time extract.pl \
    --region 1-3,3-10/1,10-12 \
    --region 45-47,47-52/1,52-54 \
    --region 87-89,89-96/1,96-98 \
    $rundir/03-mapping/aligned_reads.rev.bam $reference >fragments.rev.csv 2>fragments.rev.log

echo "bam2csv.pl (1)"
time bam2csv.pl $rundir/02-ccs/subreads_ccs.fwd.bam zmws.fwd.csv

echo "bam2csv.pl (2)"
time bam2csv.pl $rundir/02-ccs/subreads_ccs.rev.bam zmws.rev.csv
```

## Raw counts

The resulting fragment files produced in the previous step are filtered and all unique overhang (and barcode) sequences are tabulated.

```
### create output directory for overhang counts
jobdir=$rundir/05-results
mkdir -p $jobdir
cd $jobdir

echo "reporter.pl"
time reporter.pl \
    --etype b4 \
    --match 5 \
    --np 5 \
    --bcout barcodes.csv \
    --ohout overhangs.csv \
    $rundir/04-fragments
```

## Summary tables

The tabulated overhang (and barcode) sequences are processed to generate a set of output tables. See interpreation of results section below.

```
### create output directory for summary tables
jobdir=$rundir/06-summary
mkdir -p $jobdir
cd $jobdir

echo "mktable_barcode.pl"
mktable_barcode.pl $rundir/05-results/barcodes.csv > table-barcode_base_frequencies.csv

echo "mktable.pl"
mktable.pl --size 4 --prefix table $rundir/05-results/overhangs.csv
```

# Insert structure of the b4 substrate

The double stranded insert for the b4 substrate is schematically presented below. The Ligase Fidelity scripts are aimed at extracting the overhang and barcode sequences for each sequenced product. This is achieved by aligning consensus sequences of both top and bottom strands to the expected reference sequence and extracting sequence fragments corresponding to overhangs and barcodes. Additionally, the scripts extract the constant regions immediately adjacent to the overhang and barcode regions to ensure correcteness of the mapped consensus sequences.

```
         1         2         3         4         5         6         7         8         9
12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678
TTGNNNNNNCGTTGATCAATGGACGGCGCACTGGATCGCAGGTCTCCNNNNGGAGACCTGCGATCCAGTGCGCCGTCCATTGATCAACGNNNNNNCAA
---      ---                                ---    ---                                ---      ---

1..3    left constant region (1)
4..9    left internal barcode region
10..12  left constant region (2)

45..47  left constant region before overhang
48..51  overhang
52..54  right constant region after overhang

87..89  right constant region (1)
90..95  right internal barcode region
96..98  right constant region (2)
```

# Interpretation of results

# Citing

* Potapov V, Ong JL, Kucera RB, Langhorst BW, Bilotti K, Pryor JM, Cantor EJ, Canton B, Knight TF, Evans Jr TC, and Lohman GJS (2018). Comprehensive Profiling of Four Base Overhang Ligation Fidelity by T4 DNA Ligase and Application to DNA Assembly. *ACS Synth. Biol.* 7(11):2665–2674. https://doi.org/10.1021/acssynbio.8b00333
