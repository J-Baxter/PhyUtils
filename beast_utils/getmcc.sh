#!/bin/bash

# Default values
burnin=25

# Function to display usage
function usage {
    echo "Usage: $0 -i input_file [-b burnin]"
    exit 1
}

# Parse command-line options
while getopts ":i:b:" opt; do
  case ${opt} in
    i )
      input_file=$OPTARG
      ;;
    b )
      burnin=$OPTARG
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      usage
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      usage
      ;;
  esac
done
shift $((OPTIND -1))

# Check if input file is provided
if [ -z "$input_file" ]; then
    echo "Input file is required"
    usage
fi

# Input filename without extension
input_file_basename=$(basename "$input_file")
input_file_noext="${input_file_basename%.*}"

# Output filename
output_file="${input_file_noext}_mcc.tree"

# Run treeannotator
treeannotator -burnin "$burnin" -heights median "$input_file" "$output_file"

# Optionally, check if the output file exists and inform the user
if [ -f "$output_file" ]; then
    echo "Output file '$output_file' generated successfully."
else
    echo "Error: Output file '$output_file' was not generated."
fi
