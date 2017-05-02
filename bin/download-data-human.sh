#!/usr/bin/env bash
#
# This script downloads the built data files for human such that it is not 
# necessary to be built using fusioncatcher-build.py
#
# Assume that I am in fusioncatcher/bin
org="human_v87"
echo -e "\n\n\n\033[33;7m   Installing data files needed by FusionCatcher for human!   \033[0m\n"
# make sure that the working directory is the place where is this script which is fusioncatcher/bin/
hbin=""
fc=$(readlink $0)
if [[ -z $fc ]]; then
  hbin=$(dirname $0)
else
  hbin=$(dirname $fc)
fi
cd "$fbin"
# create the directory structure & cleaning
mkdir -p ../data
cd ../data
rm -rf current
rm -f "$org.tar.gz.*"
rm -f checksums.md5
rm -rf "$org"
ln -s "$org" current
wget --no-check-certificate "http://sourceforge.net/projects/fusioncatcher/files/data/$org.tar.gz.aa" -O "$org.tar.gz.aa"
wget --no-check-certificate "http://sourceforge.net/projects/fusioncatcher/files/data/$org.tar.gz.ab" -O "$org.tar.gz.ab"
wget --no-check-certificate "http://sourceforge.net/projects/fusioncatcher/files/data/$org.tar.gz.ac" -O "$org.tar.gz.ac"
wget --no-check-certificate "http://sourceforge.net/projects/fusioncatcher/files/data/$org.tar.gz.ad" -O "$org.tar.gz.ad"
wget --no-check-certificate "http://sourceforge.net/projects/fusioncatcher/files/data/checksums.md5" -O checksums.md5
md5sum -c checksums.md5
if [ "$?" -ne "0" ]; then
  echo -e "\n\n\n\033[33;7m   ERROR: The downloaded files from above have errors! MD5 checksums do not match! Please, download them again or re-run this script again!   \033[0m\n"
  exit 1
fi
cat "$org.tar.gz.*" > "$org.tar.gz"
rm -f "$org.tar.gz.*"
if ! tar -xzf "$org.tar.gz" -C .; then
    echo -e "\n\n\n\033[33;7m   ERROR: The downloaded files are corrupted!   \033[0m\n"
    exit 1
fi
rm -f "$org.tar.gz"
rm -f checksums.md5
echo -e "\n\n\n\033[33;7m   Installation went fine!   \033[0m\n"
exit 0
