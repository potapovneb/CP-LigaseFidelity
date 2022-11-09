# Ligase Fidelity Profiling

This is an official repository for computational workflow on profiling ligase fidelity using PacBio SMRT sequencing. Please see citing information below.

# Prerequisites

The computational workflow requires a number of tools to be installed and available from the command line in your system.

## Custom Python scripts

These custom scripts are provided as part of Ligase Fidelity GitHub repository. The scripts must be made available from the command line in your system (for example, by adding the scripts directory to ```$PATH``` environment variable).

* ```split.py``` - Split top and bottom subreads to separate files.
* ```summarize_results.py``` - Extract overhang pairs and generate output tables summarizing results.
* ```plot_data.py``` - Generate example ligase fidelity plots.

## PacBio tools

* ```ccs``` - Generating circular consensus sequences (CCS) based on PacBio sequencing data. Available as part of [pbbioconda](https://github.com/PacificBiosciences/pbccs).
* ```lima``` (optional) - Demultiplexing PacBio data. Available as part of [pbbioconda](https://github.com/PacificBiosciences/pbbioconda).

## Third-party tools

* ```samtools``` - Handling high-throughput sequencing data. Available through [conda](https://anaconda.org/bioconda/samtools) or [GitHub](https://github.com/samtools/samtools).
* ```pysam``` - Python module for handling sequencing files. Available through [conda](https://anaconda.org/bioconda/pysam) or [GitHub](https://github.com/pysam-developers/pysam).
* ```pandas``` - Python data analysis library. Available through [conda](https://anaconda.org/anaconda/pandas).
* ```numpy``` - Package for scientific computing with Python. Available through [conda](https://anaconda.org/anaconda/numpy).
* ```seaborn``` - Python data visualization library. Available through [conda](https://anaconda.org/anaconda/seaborn).
* ```matplotlib``` - Plotting library for the Python programming language. Available through [conda](https://anaconda.org/conda-forge/matplotlib).

# Computational workflow

Processing PacBio sequencing data in the Ligase Fidelity computational workflow proceeds through a series of steps, where output of one step serves as input for the next step. As a result, this workflow produces a set of tables summarizing ligase fidelity data.

The ```scripts/``` directory contains an example workflow (```workflow.sh```) that can be executed from the command line. Before executing the example workflow, make sure that the ```$LIGASE_FIDELITY_DIR``` environment variable in ```workflow.sh``` points to the Ligase Fidelity GitHub repository on your computer.

The computational steps are described below.

## Example PacBio BAM sequencing data

Download PacBio sequencing data (m64163_220709_012719.subreads.bam) for the b4 substrate from NCBI SRA (PRJNA894239) to the working directory and rename to movie.subreads.bam.

For convenience, we define a set of environment variables to define location of the input data files and the output directory. Please update these variables according to your computational environment.

Input PacBio BAM subreads file:
```
subreads=$PWD/movie.subreads.bam
```

Output directory:
```
rundir=$PWD/example
```

## Split top/bottom strands

In this step, ```split.py``` script splits subreads to two groups corresponding to the opposite strands of the double stranded ligation product. This is achived by examining the lengths of consecutive subreads/adapters in the PacBio polymerase reads. The default SMRTbell adapter length is 45 nt, and the subread length for the b4 substrate is 98 nt. The script looks for the longest sequence ...-[subread]-[adapter]-[subread]-[adapter]-... such that both subread and adapter lengths are within expected ranges. By default, 25% variation in expected subread and adapter lengths is allowed. This can be controlled through ```--smin```, ```--smax```, ```--amin```, ```--amax``` command line options. The subread names are saved to two separate files. Then, ```samtools``` tool is used to extract actual subread sequences and store them in two separate BAM files.

```
### create output directory for clustered reads
jobdir=$rundir/01-cluster
mkdir -p $jobdir
cd $jobdir

echo "split"
split.py \
    --subread-len 98 \
    --adapter-len 45 \
    --outfile0 subreads.0.txt \
    --outfile1 subreads.1.txt \
    ${subreads}

samtools view -N subreads.0.txt -o subreads.0.bam ${subreads}
samtools view -N subreads.1.txt -o subreads.1.bam ${subreads}
```

Note that ```samtools``` can use multiple CPU cores with the ```--threads``` option to speed up processing times.

## Build CCS sequences

PacBio ```ccs``` tool builds consensus sequences for opposite strands.

```
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
```

## Summary tables

The resulting consensus sequences for both strands are processed, all unique overhang and barcode sequences are tabulated to generate a set of output tables. For each consensus sequence the script locates the left barcode region, overhang, and right barcode region using standard pattern matching provided on the command line. For the b4 substrate, the left barcode region is six randomized bases flanked by TTG and CGT bases (```TTGNNNNNNCGT```). The corresponding pattern is ```TTG([ACGT]{6})CGT```. The overhang region is four randomized bases flanked by TCC and GGA (```TCCNNNNGGA```) and the corresponding pattern is ```TCC([ACGT]{4})GGA```. The right barcode region is six randomized bases flanked by ACG and CAA (```ACGNNNNNNCAA```) and the corresponding pattern is ```ACG([ACGT]{6})CAA```.

The full sequence of the b4 substrate is:
<u>TTG</u>NNNNNN<u>CGT</u>TGATCAATGGACGGCGCACTGGATCGCAGGTC<u>TCCNNNNGGA</u>GACCTGCGATCCAGTGCGCCGTCCATTGATCA<u>ACGNNNNNNCAA</u>.

When the consensus sequences of both strands are generated, the script applies a number of filters:
* overhang and barcode regions must strictly follow the expected patterns
* flanking bases in the opposite strands must match exactly
* at least 3 passes are required for each strand

```
TTGNNNNNNCGTTGATCAATGGACGGCGCACTGGATCGCAGGTCTCCNNNNGGAGACCTGCGATCCAGTGCGCCGTCCATTGATCAACGNNNNNNCAA (strand 0)
|||      |||                                |||    |||                                |||      |||
AACNNNNNNGCAACTAGTTACCTGCCGCGTGACCTAGCGTCCAGAGGNNNNCCTCTGGACGCTAGGTCACGCGGCAGGTAACTAGTTGCNNNNNNGTT (strand 1)
------------                                ----------                                ------------
  left bc                                    overhang                                   right bc
```

```
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
```

The last step of the workflow generates nine output files.

## Plot ligation fidelity data (optional)

The additional ```plot_data.py``` script generates four example plots based on the output tables produced in the previous step:
* ```06_matrix.png``` - Frequency heat map of all ligation events (log-scaled)
* ```07_fidelity.png``` - Stacked bar plot showing the frequency of ligation products containing each overhang
* ```08_mismatch-e.png``` - Frequency of specific base pair mismatches by position (the edge position)
* ```09_mismatch-m.png``` - Frequency of specific base pair mismatches by position (the middle position)

Please see "Interpretation of results" section below

```
cd $rundir/03-summary
plot_data.py .
```

## Interpretation of results

Note that all overhang and barcode sequences are always written in 5'-3' direction.

### [01_fragments.csv](example/03-summary/01_fragments.csv)
This is the raw output of the ligation fidelity data. Each line in this file gives:
* PacBio read name (```qname```)
* number of passes for the first strand (```np1```)
* sequence of the left barcode region (```left_bc1```), overhang (```overhang1```), and the right barcode region (```right_bc1```) in the first strand
* number of passes for the second strand (```np2```)
* sequence of the left barcode region (```left_bc2```), overhang (```overhang2```), and the right barcode region (```right_bc2```) in the second strand
* the number of mismatching bases for each overhang pair (```overhang_mismatch```).

The other output tables are built based on this raw ligation fidelity data.

### [02_overhangs.csv](example/03-summary/02_overhangs.csv)

This file tabulates the frequency of each detected overhang pair. The identity of the overhangs is provided in columns ```O1``` and ```O2```, the number of times the overhang pair was observed is provided in column ```Count```. For example, a line in this file like ```ACCG,CGGT,3909``` indciates this Watson-Crick pair was detected 3909 times in the sequencing run. Note, that overhang pair ```ACCG,CGGT``` can be viewed in two equivalent ways:
```
5' ACCG 3'
3' TGGC 5'
```

or

```
5' CGGT 3'
3' GCCA 5'
```

For the purpose of output, the overhangs in the pair are alphabetically sorted.

### [03_barcodes.csv](example/03-summary/03_barcodes.csv)

This files provides all unique barcode sequences and their frequency. This information is used to analyze nucleotide bias/composition of the substrates.

### [04_barcodes-c.csv](example/03-summary/04_barcodes-c.csv)

This file provides frequency of four bases (A, C, G, T) in every barcode position (N1, N3, N3, N4, N5, N6). The column NN provides combined frequency of four bases (irrespective of barcode position).

### [05_barcodes-p.csv](example/03-summary/04_barcodes-p.csv)

Same as above, but the frequency of four bases is provided as percentages.

### [06_matrix.csv](example/03-summary/06_matrix.csv)

This is a matrix represenatation of all ligation events presented in the ```02_overhangs.csv``` file. The top row and leftmost column provide identities of the overhang pairs, while the numbers in the matrix presents the ligation frequencies (the number of ligation events). The overhangs in the top row and leftmost column are ordered such that the diagonal corresponds to Watson-Crick (fully complementary) pairs. As mentioned above, each overhang pair can be represented in two equivalent ways. Therefore, in the matrix form each overhang pair is present twice, hence the total number of ligation events is doubled.

### [07_fidelity.csv](example/03-summary/07_fidelity.csv)

This table summarizes the number of total, correct, and mismatch ligation events for each overhang. The ratio of correct events to the total number of ligation events defines fidelity for a given overhang. Additionally, the table provides the total number of mismatch overhangs, and the five most frequent mismatch overhangs.

### [08_mismatch-e.csv](example/03-summary/08_mismatch-e.csv)

This table summarizes the frequency of mismatch bases in the "edge" position on the overhang pairs. Please check Figure 3A in the original publication below for details.

### [09_mismatch-m.csv](example/03-summary/09_mismatch-m.csv)

This table summarizes the frequency of mismatch bases in the "middle" position on the overhang pairs. Please check Figure 3B in the original publication below for details.

# Citing

* Potapov V, Ong JL, Kucera RB, Langhorst BW, Bilotti K, Pryor JM, Cantor EJ, Canton B, Knight TF, Evans Jr TC, and Lohman GJS (2018). Comprehensive Profiling of Four Base Overhang Ligation Fidelity by T4 DNA Ligase and Application to DNA Assembly. *ACS Synth. Biol.* 7(11):2665–2674. https://doi.org/10.1021/acssynbio.8b00333
