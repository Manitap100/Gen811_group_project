#!/bin/bash

echo "======================================================================"
echo "SETTING UP LEAFCUTTER ANALYSIS ENVIRONMENT"
echo "======================================================================"

# Create conda environment
echo ""
echo "[1] Creating conda environment from environment.yml..."
conda env create -f environment.yml

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create conda environment"
    exit 1
fi

echo ""
echo "[2] Activating environment..."
conda activate leafcutter_gen811

echo ""
echo "[3] Verifying installation..."

# Check Python packages
echo ""
echo "Python packages:"
python -c "import pandas; print('  ✓ pandas', pandas.__version__)"
python -c "import numpy; print('  ✓ numpy', numpy.__version__)"
python -c "import matplotlib; print('  ✓ matplotlib', matplotlib.__version__)"
python -c "import seaborn; print('  ✓ seaborn', seaborn.__version__)"
python -c "import Bio; print('  ✓ biopython', Bio.__version__)"

# Check command-line tools
echo ""
echo "Command-line tools:"
which samtools && echo "  ✓ samtools $(samtools --version | head -1)"
which bedtools && echo "  ✓ bedtools $(bedtools --version)"
which regtools && echo "  ✓ regtools $(regtools --version 2>&1 | head -1)"

# Check R
echo ""
echo "R version:"
R --version | head -1

# Check R packages
echo ""
echo "R packages:"
Rscript -e "if(require('optparse')) cat('  ✓ optparse\n')"
Rscript -e "if(require('ggplot2')) cat('  ✓ ggplot2\n')"
Rscript -e "if(require('edgeR')) cat('  ✓ edgeR\n')"

echo ""
echo "======================================================================"
echo "INSTALLATION COMPLETE"
echo "======================================================================"
echo ""
echo "To activate this environment in the future:"
echo "  conda activate leafcutter_gen811"
echo ""
echo "To deactivate:"
echo "  conda deactivate"
echo "======================================================================"
