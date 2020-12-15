#!/usr/bin/env bash

#
# FusionCatcher v1.30
#
cd "$(dirname "$(realpath -s "$0")")"

set -e
mkdir human_v102
ln -s human_v102 current
wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.aa -O human_v102.tar.gz.aa
wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ab -O human_v102.tar.gz.ab
wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ac -O human_v102.tar.gz.ac
wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ad -O human_v102.tar.gz.ad
wget --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.md5 -O human_v102.md5
md5sum -c human_v102.md5
if [ "$?" -ne "0" ]; then
  echo -e "\n\n\n\033[33;7m   ERROR: The downloaded files from above have errors! MD5 checksums do not match! Please, download them again or re-run this script again!   \033[0m\n"
  exit 1
fi
cat human_v102.tar.gz.* > human_v102.tar.gz
rm -f human_v102.tar.gz.*
if ! tar -xzf human_v102.tar.gz -C data; then
    echo -e "\n\n\n\033[33;7m   ERROR: The downloaded files are corrupted! Please, download them again or re-run this script again!   \033[0m\n"
    exit 1
fi
rm -f human_v102.tar.gz
rm -f human_v102.md5

