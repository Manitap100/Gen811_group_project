# Alternative Splicing Analysis in *Toxoplasma gondii* BDP1 Knockdown

## Authors
Yangfanxu  Manita Pariyar
GEN811
May 2026

## Background

BDP1 is a bromodomain-containing protein in *Toxoplasma gondii* that was identified to interact with chromatin remodeling factor CHD1 and splicosome subunit PRP8 through Co-IP-MS experiments (Fleck et al. 2023). Bromodomain proteins recognize acetylated histones and regulate chromatin structure. We were interested in whether BDP1 depletion affects RNA processing, particularly alternative splicing.

Alternative splicing generates multiple mRNA isoforms from a single gene and can be influenced by chromatin structure through RNA polymerase II elongation rate (Bentley 2014). The BDP1-CHD1 interaction suggests a potential link between chromatin state and splicing regulation in *Toxoplasma*.

Our RNA-seq data compared wild-type parasites (4 replicates) to BDP1 knockdown parasites (2 replicates, using an AID-inducible degron system). We used Leafcutter (Li et al. 2018), an annotation-free method, to detect differential splicing events directly from junction reads without requiring prior knowledge of transcript isoforms.

## Methods

### Data
RNA-seq BAM files were generated using STAR aligner on the RON computing cluster at UNH. Our experiment used the *Toxoplasma gondii* ME49 strain with the reference genome from ToxoDB release 68. Files are located in `data/RNASEQ/alignment/star/`.

We analyzed:
- 4 wild-type samples (WT1-WT4)
- 2 BDP1 knockdown samples (KD1-KD2)
- KD3 and KD4 were excluded due to library construction failure

### Software
We used a conda environment (`leafcutter_analysis`) on the RON cluster containing:
- Leafcutter v0.2.9 for differential splicing analysis
- regtools for junction extraction from BAM files
- SAMtools v1.17 for BAM file processing
- Python 3.8 with pandas, biopython for data analysis
- R 4.2 for statistical testing

### Analysis Pipeline

**Step 1: Junction Extraction** (`01_leafcutter_prep.sh`)
- Extracted splice junctions from BAM files using regtools
- Generated junction files for each sample
- Created sample list and group files for differential analysis

**Step 2: Leafcutter Clustering** (`02_run_leafcutter.sh`)
- Clustered junctions that share splice sites
- Performed differential splicing analysis using Fisher's exact test
- Identified significant events at FDR < 0.05

**Step 3: Effect Size Calculation** (`03_calculate_effectsize.py`)
- Calculated PSI (Percent Spliced In) for each junction:
PSI = junction_reads / cluster_total_reads × 100
ΔPSI = PSI_KD - PSI_WT
- Selected representative junction per cluster (maximum |ΔPSI|)

**Step 4: Quality Filtering** (`04_quality_filter.py`)
- Applied filtering criteria:
  - |ΔPSI| ≥ 10% (biological significance)
  - WT_reads ≥ 20, KD_reads ≥ 20 (sufficient coverage)
  - total_reads ≥ 100 (adequate depth)

**Step 5: Splice Site Validation** (`05_validate_splice_sites.py`)
- Checked for canonical splice sites (AG...GT or AG...GC) using the reference genome
- Classified events as:
  - AS candidates: with canonical splice sites (true splicing events)
  - TSS/TES candidates: without splice sites (transcription start/end changes)

**Step 6: Ranking and Visualization** (`06_rank_candidates.py`, `07_generate_figures.py`)
- Ranked AS events by effect size
- Generated summary figures


## Findings

We identified **545 significant differential splicing clusters** (FDR < 0.05) after BDP1 depletion. After applying quality filters, **162 high-confidence events** remained.
