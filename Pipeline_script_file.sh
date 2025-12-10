#!/bin/bash

# Exit immediately if any command fails
set -e


echo " Initialising conda and activating environment "


# Load conda into this shell
source /opt/conda/etc/profile.d/conda.sh

# Activate your pipeline environment
conda activate pipeline

echo "Environment activated: $CONDA_DEFAULT_ENV"
echo ""



echo "=== Running Part 2: Variant Calling Pipeline ==="


python3 Part2.py

echo ""

echo "=== Part 2 Complete!                         ==="


