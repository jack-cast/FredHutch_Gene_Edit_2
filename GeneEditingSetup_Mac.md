# Gene Editing Analysis Setup Instructions

## For Mac Users

### One-Time Setup

1. **Install Homebrew** (if not already installed):
   ```bash
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```

2. **Install Python and basic tools**:
   ```bash
   brew install python3 git wget
   ```

3. **Install Miniconda**:
   ```bash
   cd /tmp
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
   bash Miniconda3-latest-MacOSX-x86_64.sh
   source ~/.bash_profile  # or ~/.zshrc if using zsh
   ```

### For Each New Project

1. **Create project structure**:
   ```bash
   cd ~/Desktop  # or wherever you want your project
   mkdir GeneEditing
   cd GeneEditing
   mkdir NewProject
   cd NewProject
   mkdir fastq_files
   ```

2. **Copy input files**:
   - Copy `sequence.txt` to main directory
   - Copy `geneEditStart.py` to main directory
   - Copy your FASTQ files to `fastq_files/`

3. **Run initial setup**:
   ```bash
   python3 geneEditStart.py sequence.txt
   ```

4. **Download analysis scripts**:
   ```bash
   cd scripts
   wget https://raw.githubusercontent.com/jack-cast/FredHutch_Gene_Edit_2/main/geneEditFilter.py
   wget https://raw.githubusercontent.com/jack-cast/FredHutch_Gene_Edit_2/main/geneEditOutput.py
   wget https://raw.githubusercontent.com/jack-cast/FredHutch_Gene_Edit_2/main/geneEditCommon.py
   # Copy the updated genFilter2_local.py to this directory
   ```

5. **Set up environment**:
   ```bash
   cd ..
   python3 -m venv gene_editing_env
   source gene_editing_env/bin/activate
   pip install biopython pandas numpy
   conda activate base
   conda install -c bioconda pear emboss
   ```

6. **Run analysis**:
   ```bash
   cd scripts
   python3 genFilter2_local.py mk_run.txt [your_username] --run
   ```

