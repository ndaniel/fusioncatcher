#!/usr/bin/env bash


LC_ALL=C cat "$1" | echo $((`wc -l`/4))
