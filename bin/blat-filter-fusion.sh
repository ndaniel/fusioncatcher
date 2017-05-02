#!/usr/bin/env bash
SCRIPT=$(readlink -f $0)
SCRIPTPATH=$(dirname $SCRIPT)
BLATPATH="$1"
FILEDB="$2"
FILEIN="$3"
PIPE="$4"
FILEOU="$5"
shift 5
EXTRA="$@"
mkfifo $PIPE
if [[ "$BLATPATH" == "-" ]]; then
  "blat" $EXTRA $FILEDB $FILEIN $PIPE &
else
  "$BLATPATH/blat" $EXTRA $FILEDB $FILEIN $PIPE &
fi
"$SCRIPTPATH/blat-filter-fusion.py" $PIPE $FILEOU
exit 0
