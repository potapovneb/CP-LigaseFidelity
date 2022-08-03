# Dependencies

## Global environments
```
conda activate pbbioconda
module load samtools-1.3.1
```

```
qsub \
        -v root="$PROJDIR",rundir="$rundir",rname="$reference",collectionPathUri="$collectionPathUri",instrument="SEQUEL" \
        -N "demux$sampleId" \
        -o "$rundir"/workflow.log \
        -j yes \
        "$PROJDIR"/scripts/workflow-demux.sh
```

## Perl scripts
```
### cluster reads
"$root"/bin/cluster.pl

### split forward and reverse reads
"$root"/bin/split.pl

### extract barcodes and overhangs
"$root"/bin/extract.pl

### extract results
"$root"/bin/reporter.pl

### barcode table
"$root"/bin/mktable_barcode.pl

### result tables
"$root"/bin/mktable.pl

### optional
"$root"/bin/split_overhangs_by_base.pl
/mnt/home/potapov/projects/160202.LigaseFidelity/blunt_matrix_2.pl
```

## PacBio tools
```
samtools
ccs
blasr
```



## Timings
```
Cluster reads

real    65m3.296s
user    106m29.048s
sys     5m41.955s
```

```
Split reads

real    54m20.970s
user    78m1.393s
sys     11m27.359s
```

```
SAM-to-BAM

real    15m20.012s
user    13m40.762s
sys     1m2.839s
```
# CP-LigaseFidelity
