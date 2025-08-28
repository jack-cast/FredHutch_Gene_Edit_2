#!/usr/bin/env python3
"""
Local version of genFilter2.py - runs gene editing analysis locally on Windows/Mac/Linux
"""

import csv
import os
import subprocess
import sys
import argparse
import shutil
from pathlib import Path

def find_pear_executable():
    """Find PEAR executable across different platforms"""
    # First try to find it in PATH
    pear_path = shutil.which('pear')
    if pear_path:
        # Make sure it's the bioinformatics PEAR, not PHP PEAR
        try:
            result = subprocess.run([pear_path, '--help'], capture_output=True, text=True)
            if 'Paired-End reAd mergeR' in result.stdout:
                return pear_path
        except:
            pass
    
    # If not found in PATH, try common conda locations
    possible_paths = [
        os.path.expanduser('~/miniconda3/bin/pear'),
        os.path.expanduser('~/anaconda3/bin/pear'),
        os.path.expanduser('~/opt/anaconda3/bin/pear'),  # Common Mac location
        '/opt/conda/bin/pear',
        '/usr/local/bin/pear',
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            try:
                result = subprocess.run([path, '--help'], capture_output=True, text=True)
                if 'Paired-End reAd mergeR' in result.stdout:
                    return path
            except:
                continue
    
    return 'pear'  # fallback to system PATH

# Find PEAR executable at startup
PEAR_EXECUTABLE = find_pear_executable()

def run_command(cmd, description):
    """Run a command and handle errors"""
    print(f"Running: {description}")
    
    try:
        result = subprocess.run(cmd, shell=True, check=True, 
                              capture_output=True, text=True)
        print(f"  ✓ Success")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ✗ Failed with return code {e.returncode}")
        if e.stderr:
            print(f"  Error: {e.stderr.strip()}")
        return False

def process_sample(sample_info, source_dir, result_dir, run_pear=True, run_filter=True, run_output=True):
    """Process a single sample"""
    
    condition, config_file, file_R1, file_R2 = sample_info
    
    print(f"\n{'='*50}")
    print(f"Processing: {condition}")
    print(f"{'='*50}")
    
    # Set up paths - we're running from scripts/ directory
    if source_dir == "LocalTest" or source_dir == "NewProject":
        base_path = ".."
    else:
        base_path = f"../{source_dir}"
    
    # File paths
    test = condition
    file_r1 = f"{base_path}/fastq_files/{file_R1}"
    file_r2 = f"{base_path}/fastq_files/{file_R2}"
    file_p_out = f"{base_path}/stitched_reads/{test}"
    
    # Check input files exist
    if not os.path.exists(file_r1):
        print(f"  ✗ R1 file not found: {file_R1}")
        return False
    if not os.path.exists(file_r2):
        print(f"  ✗ R2 file not found: {file_R2}")
        return False
    
    success = True
    
    # Step 1: PEAR - merge paired reads
    if run_pear:
        pear_cmd = f"{PEAR_EXECUTABLE} -f {file_r1} -r {file_r2} -o {file_p_out} -v 10"
        if not run_command(pear_cmd, f"PEAR merging reads"):
            success = False
    
    # Step 2: Gene editing filter (if the script exists)
    if run_filter and success:
        filter_script = f"{base_path}/scripts/geneEditFilter.py"
        if os.path.exists(filter_script):
            filter_cmd = f"python3 {filter_script} {base_path}/ {test}.assembled {result_dir} {test} {config_file}"
            if not run_command(filter_cmd, f"Gene editing filter"):
                success = False
        else:
            print(f"  ! Filter script not found - skipping filter step")
    
    # Step 3: Gene editing output (if the script exists)
    if run_output and success:
        output_script = f"{base_path}/scripts/geneEditOutput.py"
        if os.path.exists(output_script):
            output_cmd = f"python3 {output_script} {base_path}/ {test}.assembled {result_dir} {test} {config_file}"
            if not run_command(output_cmd, f"Gene editing output"):
                success = False
        else:
            print(f"  ! Output script not found - skipping output step")
    
    return success

def main():
    parser = argparse.ArgumentParser(description='Run gene editing pipeline locally')
    parser.add_argument('inputFile', type=str, help='control file specifying input sequencing files')
    parser.add_argument('userName', type=str, help='username (for compatibility)')
    parser.add_argument('--run', dest='runScript', action='store_true', 
                       help='actually run the analysis')
    parser.add_argument('--nopear', dest='runPear', action='store_false', default=True,
                       help='skip PEAR step')
    parser.add_argument('--nofilter', dest='runFilter', action='store_false', default=True,
                       help='skip gene editing filter step')
    parser.add_argument('--nooutput', dest='runOutput', action='store_false', default=True,
                       help='skip gene editing output step')
    
    args = parser.parse_args()
    
    print(f"Gene Editing Pipeline - Local Execution")
    print(f"PEAR executable: {PEAR_EXECUTABLE}")
    
    # Read the configuration file
    files = {}
    source_dir = ""
    result_dir = ""
    
    with open(args.inputFile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            if not row or row[0].startswith('#'):
                continue
            elif row[0] == "SOURCE_DIR":
                source_dir = row[1]
            elif row[0] == "RESULT_DIR":
                result_dir = row[1]
            else:
                key = row[0]
                files[key] = row
    
    if not source_dir:
        print("ERROR: Could not find SOURCE_DIR in input file")
        return 1
    
    if not result_dir:
        result_dir = "results"
    
    print(f"Source directory: {source_dir}")
    print(f"Result directory: {result_dir}")
    print(f"Samples to process: {len(files)}")
    
    if not args.runScript:
        print("\n=== DRY RUN MODE ===")
        print("Use --run to actually execute")
        for sample_name, sample_info in files.items():
            print(f"  {sample_name}")
        return 0
    
    # Process each sample
    print(f"\n=== PROCESSING {len(files)} SAMPLES ===")
    successful = 0
    failed = 0
    
    for sample_name, sample_info in files.items():
        if process_sample(sample_info, source_dir, result_dir, 
                         args.runPear, args.runFilter, args.runOutput):
            successful += 1
            print(f"  ✓ {sample_name} completed successfully")
        else:
            failed += 1
            print(f"  ✗ {sample_name} failed")
    
    print(f"\n{'='*50}")
    print(f"FINAL SUMMARY")
    print(f"{'='*50}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"{'='*50}")
    
    return 0 if failed == 0 else 1

if __name__ == "__main__":
    sys.exit(main())