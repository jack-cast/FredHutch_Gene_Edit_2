#
# geneEdidStart.py [inputFile]
#
#    run from the base of a geneEditing experiment
#       expects one dir fastq_files with R1 and R2 reads
#       input file is:
#    
# start_primer
# reference
# end_primer
# guide
#
from sys import argv
import csv
import os
import time
import argparse
from pathlib import Path
from os import listdir
from os.path import isfile, join

ge_base   = "/fh/fast/adair_j/grp/jadair_shared/PROJECTS/GeneEditing/"
ge_script = "/fh/fast/adair_j/grp/jadair_shared/PROJECTS/GeneEditing/scripts"
 
#-----------------------------------------------------------------------------------------------
# main start
#-----------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Build script files to run gene editing pipeline')
parser.add_argument('inputFile',type=str,help='control file specifying input sequencing files and editing detials')
args = parser.parse_args()
inputFile = args.inputFile
#
# read input file
#
locusList = []
fh = (line.strip() for line in open(inputFile,'r'))
try:
  while True:
    locus = next(fh)
    start_primer = next(fh)
    refseq = next(fh)
    end_primer = next(fh)
    guideseq = next(fh)
    guideseq = guideseq.lower()
    print(f'locus        : {locus}')
    print(f'start_primer : {start_primer}')
    print(f'refseq       : {refseq}')
    print(f'end_primer   : {end_primer}')
    print(f'guideseq     : {guideseq}')
    locusList.append((locus,start_primer,refseq,end_primer,guideseq))
except StopIteration:
  pass
#
# make sure directory structure exists
#
my_file = Path("./fastq_files")
if not my_file.is_dir():
  print('fastq_files directory required')
  quit()

my_file = Path("./scripts")
if not my_file.is_dir():
  os.mkdir('./scripts')


my_file = Path("./stitched_reads")
if not my_file.is_dir():
  os.mkdir('./stitched_reads')


my_file = Path("./aligned_reads")
if not my_file.is_dir():
  os.mkdir('./aligned_reads')


my_file = Path("./results")
if not my_file.is_dir():
  os.mkdir('./results')

#
# get all fastq files
#
#
#  look for all .....ge.csv files in dir
#
fid = '_R1_'
mypath = './fastq_files/'
onlyFiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
r1Files = [f for f in onlyFiles if fid in f]
#
# make mk_run.txt file
#
path, srcDir = os.path.split(os.getcwd())
if len(locusList) > 1:
  config = 'geConfig_' + '[FIX]' + '.txt'
  print("!!! user must fix mk_run.txt to set correct geConfig locus for each file")
else:
  locus,t1,t2,t3,t4 = locusList[0]
  config = 'geConfig_' + locus + '.txt'

with open('./scripts/mk_run.txt','w') as fh:
  fh.write('#\n')
  fh.write('# text file to generate scripts for gene edit processing\n')
  fh.write('#\n')
  fh.write(f'SOURCE_DIR,{srcDir}\n')
  fh.write(f'RESULT_DIR,results\n')
  fh.write('#\n')
  fh.write('#[condition],[vector],[srcfile_R1],[srcfile_R2]...\n')
  fh.write('#\n')
  for fr1 in r1Files:
    loc = fr1.find('_R1_')
    fr2 = fr1[:loc] + '_R2_' + fr1[(loc+4):]
    # replace _R1_ with _R2_
    fh.write(f'{fr1[:loc]},{config},{fr1},{fr2}\n')
#
# make geConfig.txt file
#

for (locus,start_primer,refseq,end_primer,guideseq) in locusList:

  with open('./scripts/geConfig_' + locus + '.txt','w') as fh:
    fh.write('#\n')
    fh.write(f'name:{locus}\n')
    fh.write(f'start_primer:{start_primer}\n')
    fh.write(f'refseq:{refseq}\n')
    fh.write(f'guide_seq:{guideseq}\n')
    fh.write(f'end_primer:{end_primer}\n')
#
# make sh file to run
#
with open('./scripts/run.sh','w') as fh:
  fh.write('python3 genFilter2.py mk_run.txt jcastell --run\n')
