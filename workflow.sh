#!/bin/bash

samples=$1

PROJDIR=$PWD

cd $PROJDIR

while read line
do
    if [[ ! $line =~ "SampleID" ]] && [[ ! $line =~ "#" ]]
    then
	### sample details
	sampleId=`echo $line | cut -d, -f1`
	reference=`echo $line | cut -d, -f2`
	collectionPathUri=`echo $line | cut -d, -f3`
	# collectionPathUri=`echo $line | cut -d, -f7`

	### preview run info
	echo ""
	echo "sampleId=$sampleId"
	echo "reference=$reference"
	echo "collectionPathUri=$collectionPathUri"

	### output directory
	rundir=`printf "%s/samples-kb/%05i" $PROJDIR $sampleId`
	# rm -rf "$rundir"
	mkdir -p "$rundir"

	# ### run demultiplexing script
	# qsub \
	#     -v root="$PROJDIR",rundir="$rundir",rname="$reference",collectionPathUri="$collectionPathUri",instrument="SEQUEL" \
	#     -N "demux$sampleId" \
	#     -o "$rundir"/workflow.log \
	#     -j yes \
	#     "$PROJDIR"/scripts/workflow-demux.sh

	if [[ ! "$reference" =~ ";" ]]
	then
	    ### run demultiplexing script
	    qsub \
		    -v root="$PROJDIR",rundir="$rundir",rname="$reference",collectionPathUri="$collectionPathUri",instrument="SEQUEL" \
		    -N "demux$sampleId" \
		    -o "$rundir"/workflow.log \
		    -j yes \
		    "$PROJDIR"/scripts/workflow-demux.sh
	else
	    IFS=';' read -r -a refArray <<< "$reference"
	    
	    for ref_i in "${refArray[@]}"
	    do
		rundir_i="$rundir/$ref_i"
		mkdir -p "$rundir_i"

		### run demultiplexing script (no log)
		qsub \
		    -v root="$PROJDIR",rundir="$rundir_i",rname="$ref_i",collectionPathUri="$collectionPathUri",instrument="SEQUEL" \
		    -N "demux$sampleId-$ref_i" \
		    -o /dev/null \
		    -j yes \
		    "$PROJDIR"/scripts/workflow-demux.sh
	    done
	fi
    fi
done < $samples
