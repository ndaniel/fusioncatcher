#!/usr/bin/env bash

infile1="$1"
infile2="$2"
oufile="/dev/stdout"

if [ $# -eq 3 ] ; then
  oufile="$3"
fi


LC_ALL=C \
paste \
-d '\n' \
<(LC_ALL=C cat $infile1 | paste - - - -) \
<(LC_ALL=C cat $infile2 | paste - - - -) \
| \
LC_ALL=C \
tr \
'\t' \
'\n' \
> $oufile
