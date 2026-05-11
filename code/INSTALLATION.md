# Installation Instructions

## Step 1: Create Conda Environment

```bash
# From the Gen811_group_project directory
bash setup_environment.sh
```

This will:
- Create a new conda environment called `leafcutter_gen811`
- Install all required Python packages
- Install R and required R packages
- Install SAMtools, bedtools, regtools
- Verify all installations

## Step 2: Install Leafcutter

```bash
# Activate the environment
conda activate leafcutter_gen811

# Install Leafcutter from source
bash install_leafcutter.sh
```

This will:
- Clone Leafcutter from GitHub to `~/leafcutter`
- Install the R package
- Verify key scripts are present

## Step 3: Verify Installation

```bash
# Activate environment
conda activate leafcutter_gen811

# Check Python
python --version
python -c "import pandas, numpy, matplotlib, seaborn, Bio; print('Python packages OK')"

# Check command-line tools
samtools --version
bedtools --version
regtools --version

# Check R
R --version
Rscript -e "library(optparse); library(edgeR); print('R packages OK')"

# Check Leafcutter
ls ~/leafcutter/clustering/leafcutter_cluster_regtools.py
ls ~/leafcutter/scripts/leafcutter_ds.R
```

## Step 4: Update Scripts with Leafcutter Path

If Leafcutter is installed in a different location, update `code/02_run_leafcutter.sh`:

```bash
# Line 18-19 in 02_run_leafcutter.sh
LEAFCUTTER_DIR="${HOME}/leafcutter"  # Update this path if needed
```

## Troubleshooting

### Issue: regtools not found
```bash
conda install -c bioconda regtools
```

### Issue: R packages missing
```bash
conda activate leafcutter_gen811
Rscript -e "install.packages('optparse', repos='http://cran.r-project.org')"
Rscript -e "BiocManager::install('edgeR')"
```

### Issue: Leafcutter scripts not executable
```bash
cd ~/leafcutter
chmod +x clustering/*.py
chmod +x scripts/*.R
```

## Environment Management

### Activate environment
```bash
conda activate leafcutter_gen811
```

### Deactivate environment
```bash
conda deactivate
```

### Remove environment (if needed)
```bash
conda env remove -n leafcutter_gen811
```

### Export environment (for sharing)
```bash
conda env export > environment_full.yml
```

## RON Cluster Specific Notes

If on RON cluster:
1. Load conda module first: `module load anaconda`
2. Create environment in your home directory
3. Check available disk space: `quota -s`
4. Junction files can be large (~500MB each)

## Testing the Installation

Run a quick test:
```bash
conda activate leafcutter_gen811
cd ~/Gen811_group_project/code
bash 01_leafcutter_prep.sh
# Should complete without errors
```
