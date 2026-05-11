#!/bin/bash

#####################################################################
# Master script to run complete Leafcutter analysis pipeline
#####################################################################

set -e
set -u

echo "======================================================================"
echo "COMPLETE LEAFCUTTER ANALYSIS PIPELINE"
echo "======================================================================"
echo "Start time: $(date)"
echo ""

# Check if in correct directory
if [ ! -f "01_leafcutter_prep.sh" ]; then
    echo "ERROR: Please run this script from the code/ directory"
    exit 1
fi

# Stage 1: Leafcutter preparation
echo ""
echo ">>> STAGE 1: LEAFCUTTER PREPARATION"
bash 01_leafcutter_prep.sh
if [ $? -ne 0 ]; then
    echo "ERROR: Stage 1 failed"
    exit 1
fi

# Stage 2: Leafcutter differential analysis
echo ""
echo ">>> STAGE 2: LEAFCUTTER ANALYSIS"
bash 02_run_leafcutter.sh
if [ $? -ne 0 ]; then
    echo "ERROR: Stage 2 failed"
    exit 1
fi

# Stage 3: Effect size calculation
echo ""
echo ">>> STAGE 3: EFFECT SIZE CALCULATION"
python 03_calculate_effectsize.py
if [ $? -ne 0 ]; then
    echo "ERROR: Stage 3 failed"
    exit 1
fi

# Stage 4: Quality filtering
echo ""
echo ">>> STAGE 4: QUALITY FILTERING"
python 04_quality_filter.py
if [ $? -ne 0 ]; then
    echo "ERROR: Stage 4 failed"
    exit 1
fi

# Stage 5: Splice site validation
echo ""
echo ">>> STAGE 5: SPLICE SITE VALIDATION"
python 05_validate_splice_sites.py
if [ $? -ne 0 ]; then
    echo "ERROR: Stage 5 failed"
    exit 1
fi

# Stage 6: Ranking candidates
echo ""
echo ">>> STAGE 6: RANKING CANDIDATES"
python 06_rank_candidates.py
if [ $? -ne 0 ]; then
    echo "ERROR: Stage 6 failed"
    exit 1
fi

# Stage 7: Generate figures
echo ""
echo ">>> STAGE 7: GENERATING FIGURES"
python 07_generate_figures.py
if [ $? -ne 0 ]; then
    echo "ERROR: Stage 7 failed"
    exit 1
fi

echo ""
echo "======================================================================"
echo "PIPELINE COMPLETE"
echo "======================================================================"
echo "End time: $(date)"
echo ""
echo "Results location:"
echo "  Tables:  ../results/tables/"
echo "  Figures: ../results/figures/"
echo ""
echo "Key output files:"
echo "  - effect_sizes.txt"
echo "  - quality_filtered.txt"
echo "  - AS_candidates.txt"
echo "  - TSS_TES_candidates.txt"
echo "  - ranked_AS_candidates.txt"
echo "  - AS_vs_TSSTES_ratio.png"
echo "  - effect_size_distribution.png"
echo "  - top_candidates_heatmap.png"
echo "======================================================================"
