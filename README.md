# FredHutch_Gene_Edit_2
Fred Hutch CGT Custom gene editing code version 2

## Quick Start

1. **Setup Instructions:**
   - [Windows (WSL2) Setup](setup/GeneEditingSetup_Windows.md)
   - [Mac Setup](setup/GeneEditingSetup_Mac.md)

2. **Required Files:**
   - Your FASTQ files (paired-end, R1/R2)
   - `sequence.txt` with your gene editing configuration ([example](setup/example_sequence.txt))

3. **Run Analysis:**
   ```bash
   python3 geneEditStart.py sequence.txt
   cd scripts
   python3 genFilter2_local.py mk_run.txt your_username --run
