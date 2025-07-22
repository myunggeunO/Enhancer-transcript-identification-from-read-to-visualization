#!/bin/bash

# 1. Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda

# 2. Set PATH
export PATH="$HOME/miniconda/bin:$PATH"
echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc

# 3. Initialize conda
source "$HOME/miniconda/bin/activate"
conda init bash
source ~/.bashrc

# 4. Update conda
conda update -n base -c defaults conda -y

# 5. Install mamba in the base environment
conda install -n base -c conda-forge mamba -y

# 6. Create a virtual environment named "enhancer-env"
mamba create -n enhancer-env python=3.9 numpy pandas r-base -y

# 7. Launch a new shell with the "enhancer-env" environment activated
echo "Environment 'enhancer-env' created. Opening new shell with it activated."
exec bash --rcfile <(echo 'source ~/.bashrc && conda activate enhancer-env')
