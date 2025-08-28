# Gene Editing Analysis Setup Instructions

## For Windows Users (WSL2)

### One-Time WSL2 Setup

1. **Install WSL2** (if not already done):
   ```powershell
   # Run in PowerShell as Administrator
   dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
   dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
   # Restart computer
   wsl --set-default-version 2
   ```

2. **Install Ubuntu** from Microsoft Store

3. **Install Miniconda** (first time only):
   ```bash
   cd /tmp
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   source ~/.bashrc
   ```

### For Each New Project

1. **Create project structure**:
   ```bash
   # Navigate to your Windows folder via WSL
   cd /mnt/c/Users/[YourUsername]/Desktop/GeneEditing/
   mkdir NewProject
   cd NewProject
   mkdir fastq_files scripts stitched_reads aligned_reads results
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
   python3 genFilter2_local.py mk_run.txt [your_username] --run --nopear
   ```

7. **Check results**:
   ```bash
   ls -la ../results/
   ```