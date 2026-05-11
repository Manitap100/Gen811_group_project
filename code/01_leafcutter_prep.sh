#!/bin/bash

#####################################################################
# Script: 01_leafcutter_prep.sh
# Purpose: Prepare data for Leafcutter analysis
# - Extract junctions from BAM files using regtools
# - Create sample list and group files
#####################################################################

set -e
set -u

echo "======================================================================"
echo "STEP 1: LEAFCUTTER PREPARATION"
echo "======================================================================"

# Define paths
PROJECT_DIR="/home/users/yx1040/Gen811_group_project"
BAM_DIR="${PROJECT_DIR}/data/RNASEQ/alignment/star"
OUTPUT_DIR="${PROJECT_DIR}/data/leafcutter"
JUNC_DIR="${OUTPUT_DIR}/juncfiles"

# Create output directories
mkdir -p ${JUNC_DIR}

# Define samples (excluding KD3 and KD4)
WT_SAMPLES="WT1 WT2 WT3 WT4"
KD_SAMPLES="KD1 KD2"

echo ""
echo "[1.1] Extracting junctions from BAM files..."
echo "----------------------------------------------------------------------"

# Extract junctions for each sample
for sample in ${WT_SAMPLES} ${KD_SAMPLES}; do
    echo "Processing: ${sample}"
    
    BAM_FILE="${BAM_DIR}/${sample}_Aligned.sortedByCoord.out.bam"
    JUNC_FILE="${JUNC_DIR}/${sample}.junc"
    
    if [ ! -f "${BAM_FILE}" ]; then
        echo "ERROR: BAM file not found: ${BAM_FILE}"
        exit 1
    fi
    
    # Extract junctions using regtools
    # -a 8: minimum anchor length
    # -m 50: minimum intron size
    # -M 500000: maximum intron size
    # -s XS: use XS tags from STAR aligner
    regtools junctions extract \
        -a 8 \
        -m 50 \
        -M 500000 \
        -s XS \
        ${BAM_FILE} \
        -o ${JUNC_FILE}
    
    echo "  ✓ ${JUNC_FILE}"
    echo "  Junctions: $(wc -l < ${JUNC_FILE})"
done

echo ""
echo "[1.2] Creating sample list..."
echo "----------------------------------------------------------------------"

JUNCFILE_LIST="${JUNC_DIR}/juncfile_list.txt"
> ${JUNCFILE_LIST}

for sample in ${WT_SAMPLES} ${KD_SAMPLES}; do
    echo "${JUNC_DIR}/${sample}.junc" >> ${JUNCFILE_LIST}
done

echo "  ✓ ${JUNCFILE_LIST}"
echo "  Samples: $(wc -l < ${JUNCFILE_LIST})"

echo ""
echo "[1.3] Creating groups file..."
echo "----------------------------------------------------------------------"

GROUPS_FILE="${JUNC_DIR}/groups_file.txt"
> ${GROUPS_FILE}

for sample in ${WT_SAMPLES}; do
    echo "${sample} WT" >> ${GROUPS_FILE}
done

for sample in ${KD_SAMPLES}; do
    echo "${sample} KD" >> ${GROUPS_FILE}
done

echo "  ✓ ${GROUPS_FILE}"
cat ${GROUPS_FILE}

echo ""
echo "======================================================================"
echo "STEP 1 COMPLETE"
echo "======================================================================"
echo "Output: ${JUNC_DIR}"
echo "Ready for Leafcutter clustering"
echo "======================================================================"
