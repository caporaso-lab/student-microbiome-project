#!/bin/bash
# QIIME 1.7.0-dev, master@72514df
# R 2.12.0
# vegan 2.0-7
# SMP: bb888861d87c9d10b33121d829f672df01aace45
# map.txt MD5: 8de601815ea5f4f9b0900794407a42be
# Author: Jai Ram Rideout

# Run several methods in compare_categories.py over per-body-site
# time-series-only distance matrices using the PersonalID category.
DMDIR=../../beta-diversity
MAP=map.txt
CATEGORY=PersonalID
NUMPERMS=999
for ZIPFILE in $DMDIR/*weighted_unifrac_dm.*_ts_only.txt.gz
do
  DMFILE=$(echo $ZIPFILE | sed 's/.gz//')
  OUTDIR=$(basename $DMFILE | sed 's/.txt//')_stats
  mkdir $OUTDIR
  gunzip -c $ZIPFILE > $DMFILE
  compare_categories.py -i $DMFILE -m $MAP -c $CATEGORY -o $OUTDIR -n $NUMPERMS --method adonis
  compare_categories.py -i $DMFILE -m $MAP -c $CATEGORY -o $OUTDIR -n $NUMPERMS --method anosim
  compare_categories.py -i $DMFILE -m $MAP -c $CATEGORY -o $OUTDIR -n $NUMPERMS --method permdisp
done
