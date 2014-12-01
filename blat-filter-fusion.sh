#!/usr/bin/env bash
SCRIPT=$(readlink -f $0)
SCRIPTPATH=$(dirname $SCRIPT)
FILEDB="$1"
FILEIN="$2"
PIPE="$3"
FILEOU="$4"
shift 4
EXTRA="$@"
mkfifo $PIPE
"blat" $EXTRA $FILEDB $FILEIN $PIPE &
"$SCRIPTPATH/blat-filter-fusion.py" $PIPE $FILEOU
exit 0
