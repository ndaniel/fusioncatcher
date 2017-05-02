#!/usr/bin/env bash
# This is a wrapper around 'seqtk subseq' for cases when there is not enough free memory
# it requires GNU parallel

# cat list-names-reads-filtered_genome.txt | gnu_parallel --part -k -j1 --block 10G seqtk subseq reads-filtered.fq - > reads_filtered_unique-mapped-genome.fq

SEQTKPATH="$1"
PARALLELPATH="$2"
NOSPLITS="$3"
FQIN="$4"
FIDS="$5"
FQOU="$6"


NEWSIZE=$((`stat --printf="%s" $FQIN`/$NOSPLITS + 1))

#NEWLINES=$((`wc -l FQIN`/$NOSPLITS + 1))



#cat IDS | parallel --pipe -k -j1 --block 10G seqtk subseq FQIN - > FQOU

#parallel --pipepart -k -j1 --block NEWSIZE -a FIDS seqtk subseq FQIN - > FQOU

#cat $FIDS | parallel --pipe -k -j1 --block $NEWSIZE seqtk subseq $FQIN - > $FQOU

if [[ "$PARALLELPATH" == "-" ]]; then
    if [[ "$SEQTKPATH" == "-" ]]; then
        parallel -a "$FIDS" --pipepart -k -j1 --block $NEWSIZE seqtk subseq "$FQIN" - > "$FQOU"
    else
        parallel -a "$FIDS" --pipepart -k -j1 --block $NEWSIZE "$SEQTKPATH/seqtk" subseq "$FQIN" - > "$FQOU"
    fi
else
    if [[ "$SEQTKPATH" == "-" ]]; then
        "$PARALLELPATH/parallel" -a "$FIDS" --pipepart -k -j1 --block $NEWSIZE seqtk subseq "$FQIN" - > "$FQOU"
    else
        "$PARALLELPATH/parallel" -a "$FIDS" --pipepart -k -j1 --block $NEWSIZE "$SEQTKPATH/seqtk" subseq "$FQIN" - > "$FQOU"
    fi
fi

# ORIGINAL
# parallel -a "$FIDS" --pipepart -k -j1 --block $NEWSIZE seqtk subseq "$FQIN" - > "$FQOU"




