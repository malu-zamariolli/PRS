#!/bin/bash
gunzip -c $1 | \
awk '{seen[$3]++; if(seen[$3]==1){ print}}' | \
gzip - > $2

#Count number of time SNP ID was observed (seen[$3]++). If this is the first time seeing this SNP ID, print it.