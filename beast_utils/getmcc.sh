#!/bin/bash

# Check if a directory is provided as an argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

directory=$1

# Get arguments for burnin/ resample parameters
while getopts "burnin:resample:" var
do
   case "$var" in
       i) burnin=${OPTARG};;
       r) resample=${OPTARG};;
   esac
done

for combinedfile in *combined.trees; do 
  /Applications/BEAST\ v1.10.4/bin/treeannotator -burnin "$burnin" -heights mean $combinedfile "$(echo "$combinedfile" | awk '{sub(/combined.log$/, "mcc.tree")}1')";
done