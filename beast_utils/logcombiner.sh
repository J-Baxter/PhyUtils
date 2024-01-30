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

# Get a list of unique prefixes
unique_prefixes=($(ls "$directory" | sed 's/_.*//' | sort -u))

# Loop over unique prefixes
for prefix in "${unique_prefixes[@]}"; do
    echo "Prefix: $prefix"
    
    # Get files
    treefiles=($(find "$directory" -type f -name "${prefix}*trees" -print))
    logfiles=($(find "$directory" -type f -name "${prefix}*log" -print))
    
    # Create output filename
    outputtree="$(echo "$treefiles" | awk '{sub(/[0-9]$/, "combined")}1')" # to loop
    outputlog="$(echo "$treefiles" | awk '{sub(/[0-9]$/, "combined")}1')"  # to loop
    
    # concatenate input strings
    
    
    # You can perform additional actions for each unique prefix here
    /Applications/BEAST\ v1.10.4/bin/logcombiner -burnin "$burnin" -trees -resample "$resample" "$treefiles" "$outputtree"
    
    /Applications/BEAST\ v1.10.4/bin/logcombiner -burnin "$burnin" -resample "$resample" "$logfiles" "$outputlog"
    
done


