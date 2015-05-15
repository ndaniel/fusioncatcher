#!/usr/bin/env bash

fbin=""
fc=$(readlink $0)
if [[ -z $fc ]]; then
 fbin=$(dirname $0)
else
 fbin=$(dirname $fc)
fi

rm -rf "test_fusioncatcher"

"$fbin/fusioncatcher.py" --input "$fbin/../test" --output "test_fusioncatcher"

LC_ALL=C sort "$fbin/../test/summary_candidate_fusions.txt" > "test_fusioncatcher/summary_candidate_fusions_a.txt"

LC_ALL=C sort "test_fusioncatcher/summary_candidate_fusions.txt" > "test_fusioncatcher/summary_candidate_fusions_b.txt"

if diff "test_fusioncatcher/summary_candidate_fusions_a.txt" "test_fusioncatcher/summary_candidate_fusions_b.txt" >/dev/null ; then
  echo -e "\n\n\n\033[33;7m   Installation test went fine! Installation is ok!   \033[0m\n"
else
  echo -e "\n\n\n\033[33;7m   WARNING: Test is NOT ok! There is something wrong with FusionCatcher installation!   \033[0m\n"
fi

rm -f "test_fusioncatcher/summary_candidate_fusions_a.txt"
rm -f "test_fusioncatcher/summary_candidate_fusions_b.txt"
