#!/bin/bash

# =========================================================
# PANDA Network Generation Script for LUAD and GTEx Samples
# This script uses the netzoopy panda command line tool
# =========================================================

# Define the common prior files for clarity
MOTIF_PRIOR="~/159009/data/motif_997TFs_ensembl.txt"
PPI_PRIOR="~/159009/data/ppi_997TFs.txt"

# Common PANDA arguments
PANDA_ARGS="--mode_process intersection --with_header"

## --- 1. Run LUAD Network Generation ---
echo "Starting PANDA run for LUAD (Cancer Samples)..."
netzoopy panda \
  -e luad_corrected_expression_final.txt \
  -m ${MOTIF_PRIOR} \
  -p ${PPI_PRIOR} \
  ${PANDA_ARGS} \
  -o ~/output_panda_LUAD_final.txt

if [ $? -eq 0 ]; then
  echo "LUAD PANDA network successfully generated."
else
  echo "Error: LUAD PANDA network generation failed."
fi

echo "------------------------------------------------------"

## --- 2. Run GTEx Network Generation ---
echo "Starting PANDA run for GTEx (Normal Samples)..."
netzoopy panda \
  -e gtex_corrected_expression_final.txt \
  -m ${MOTIF_PRIOR} \
  -p ${PPI_PRIOR} \
  ${PANDA_ARGS} \
  -o ~/output_panda_GTEx_final.txt

if [ $? -eq 0 ]; then
  echo "GTEx PANDA network successfully generated."
else
  echo "Error: GTEx PANDA network generation failed."
fi

echo "PANDA script finished."