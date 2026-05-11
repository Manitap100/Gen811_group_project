#!/bin/bash

#####################################################################
# Script: 02_run_leafcutter.sh
# Purpose: Run Leafcutter clustering and differential analysis
#####################################################################

set -e
set -u

echo "======================================================================"
echo "STEP 2: LEAFCUTTER DIFFERENTIAL SPLICING ANALYSIS"
echo "======================================================================"

# Define paths
PROJECT_DIR="/home/users/yx1040/Gen811_group_project"
JUNC_DIR="${PROJECT_DIR}/data/leafcutter/juncfiles"
OUTPUT_DIR="${PROJECT_DIR}/data/leafcutter"

# Leafcutter paths (adjust if needed)
LEAFCUTTER_DIR="${HOME}/leafcutter"
CLUSTER_SCRIPT="${LEAFCUTTER_DIR}/clustering/leafcutter_cluster_regtools.py"
DIFF_SCRIPT="${LEAFCUTTER_DIR}/scripts/leafcutter_ds.R"

# Check Leafcutter installation
if [ ! -f "${CLUSTER_SCRIPT}" ]; then
    echo "ERROR: Leafcutter not found at ${LEAFCUTTER_DIR}"
    echo "Please install Leafcutter or update LEAFCUTTER_DIR"
    exit 1
fi

cd ${OUTPUT_DIR}

echo ""
echo "[2.1] Clustering introns..."
echo "----------------------------------------------------------------------"

python ${CLUSTER_SCRIPT} \
    -j ${JUNC_DIR}/juncfile_list.txt \
    -m 50 \
    -o leafcutter \
    -l 500000

echo "  ✓ Clustering complete"

if [ ! -f "leafcutter_perind.counts.gz" ]; then
    echo "ERROR: Clustering failed"
    exit 1
fi

echo ""
echo "[2.2] Differential splicing analysis..."
echo "----------------------------------------------------------------------"

Rscript ${DIFF_SCRIPT} \
    --num_threads 4 \
    --FDR 0.05 \
    leafcutter_perind.counts.gz \
    ${JUNC_DIR}/groups_file.txt \
    --output_prefix diff_splicing

echo "  ✓ Analysis complete"

# Rename for clarity
if [ -f "diff_splicing_cluster_significance.txt" ]; then
    cp diff_splicing_cluster_significance.txt diff_splicing_results.txt
fi

echo ""
echo "[2.3] Summary..."
echo "----------------------------------------------------------------------"

TOTAL=$(wc -l < diff_splicing_results.txt)
TOTAL=$((TOTAL - 1))
SIG=$(awk 'NR>1 && $5<0.05' diff_splicing_results.txt | wc -l)

echo "Total testable: ${TOTAL}"
echo "Significant (FDR<0.05): ${SIG}"

echo ""
echo "======================================================================"
echo "STEP 2 COMPLETE"
echo "======================================================================"
echo "Output:"
echo "  - leafcutter_perind.counts.gz"
echo "  - diff_splicing_results.txt"
echo "======================================================================"
