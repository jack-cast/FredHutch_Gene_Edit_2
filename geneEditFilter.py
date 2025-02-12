#
# easyGeneEdit.py
#
# Lauren and Mark   3-19-18
# Dnyanada          11-13-18, 11-16-18, 11-27-18
#
# useage: easyGeneEdit.py [dir] [fileNameBase] [primer set name]
#
# must run pear previous to running this script to assemble double-ended reads
# <single-ended reads may be placed in stitched_reads sub-directory
# expects fastq file in sub-directory /stitched_reads/fileBaseName.fastq
# 
# This script will filter gene-edit reads for Quality and exact match of primers
# Good Q + primer reads are trimmed <by primer or by optional length>
# Identical reads are combined and marked with an ID
# Combined reads are written to a fasta file
# Needle is run on the fasta file to compute alignments (ml emboss)
# Needle results are loaded and the alignment string is matched to the sequence ID
# Aligned sequences are compared to the reference sequence to identify substitution errors 
# 
import sys
import math
import re
import csv
import subprocess
import os
import uuid
from subprocess import call
from collections import OrderedDict
from collections import defaultdict
from collections import namedtuple
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Emboss.Applications import NeedleCommandline
import multiprocessing
import argparse
from geneEditCommon import *
#from Bio import SeqIO
#--------------------------------------------------------------------------------
#
# checkQuality: how many bases are below Q 
# 
#--------------------------------------------------------------------------------
def checkQuality(qualSeq):
  #
  # count how many bases are below low quality
  #
  fail = 0
  i = 0
  for char in qualSeq:
    q = ord(char) - ord('!')
    if (q < 15):
      fail += 1
  return fail
#--------------------------------------------------------------------------------
# checkPrimers
#
#   check for start and end primer with an allowed error if 2 bases, no indel
#     this is meant in gene editing to be the entire sequence from primer to
#     start of target, and from end of target to end of last primer
#
#   return trimmed sequence or ""
#
#--------------------------------------------------------------------------------
def checkPrimers(apt,start_primer,end_primer,refseq):
    #
    # exact match
    #
    startPrimerLoc = apt.find(start_primer)
    
    if startPrimerLoc != 0:
      #
      # slow check...allow error
      # globalxs = global alignment
      # x = No parameters. Identical characters have a score of 1, otherwise 0. (Parameter for matches)
      # s = Same open and extend gap penalties for both sequences. (Parameter for gap penalties = -4, -1)
      #
      l = len(start_primer) - 2
      alignments = pairwise2.align.globalxs(apt[0:len(start_primer)],start_primer,-4,-1)

      score = alignments[0][2]
      if (score < l):
        #print("Fail start primer",score,len(start_primer))
        #print(alignments[0][0])    
        #print(alignments[0][1])    
        #print(alignments[0][2])    
        #print(alignments[0][3])    
        #print(alignments[0][4])    
        return("")
    #
    # exact match
    #
    endPrimerLoc   = apt.find(end_primer)
    #print("end primer loc = {}".format(endPrimerLoc))
    if endPrimerLoc >= (len(start_primer) + len(refseq)/2):
         #print("good exact match")
         return(apt[len(start_primer) : endPrimerLoc])
    #
    # mismatch check, allow 1 error in 20      
    #
    errors = int(len(end_primer)/20)
    l = len(end_primer) - errors
    alignments = pairwise2.align.localxs(apt[len(start_primer):],end_primer,-4,-1)

    score = alignments[0][2]
    endPrimerLoc = alignments[0][3]
    #print("align score = {} len = {}".format(score,len(end_primer)))    
    if (score >= l):
       endPrimerLoc += len(start_primer) # removed for alignment
       #print("Good end primer alignment, loc = {}".format(endPrimerLoc))
       return(apt[len(start_primer) : endPrimerLoc])

    #print("Fail end primer",score,len(end_primer))
    #print(alignments[0][0])    
    #print(alignments[0][1])    
    return("")
#-------------------------------------------------------------------------------------
#
# logFailures - write primer fails to file if count > threshold
#
#
#-------------------------------------------------------------------------------------
def logFailures(logFile,stat,failList,start_primer,refseq,end_primer):
  sortFail = sorted(failList,key=lambda al: al[1],reverse=True)  # largest first
  with open(logFile,'w') as fh:
    fh.write("Total input sequences         = {}\n".format(stat.totalInputSequences))
    fh.write("Fail Q                        = {}\n".format(stat.failQuality))
    fh.write("Pass Primer chech & Alignment = {}\n".format(stat.goodAlignments))
    for k,v in sortFail:
      if (v >= 20):
        fh.write(">count = {}\n".format(v))
        glob = start_primer+refseq+end_primer
        minL = min(len(k),len(glob))
        kL = []
        for i in range(0,minL):
          if(glob[i] == k[i]):
             kL.append(k[i])
          else:
             kL.append(k[i].lower())
        k = "".join(kL)
        l1 = len(start_primer)
        l2 = l1 + len(refseq)
        fh.write("ref start primer {}\n".format(glob[0:l1]))
        fh.write("tst start primer {}\n".format(k[0:l1]))
        fh.write("ref target seq   {}\n".format(glob[l1:l2]))
        fh.write("tst target seq   {}\n".format(k[l1:l2]))
        fh.write("ref end primer   {}\n".format(glob[l2:]))
        fh.write("tst end primer   {}\n".format(k[l2:]))


#---------------------------------------------------------------------------
#
# execNeedle
#    manage and execute call to needle
#
#
#---------------------------------------------------------------------------
def execNeedle(seqDict,resultList):
  if len(seqDict.items()) == 0:
    print('execNeedle: no sequences found')
    return
  #
  # prepare to run needle
  #
  try:
    tmpFile = str(uuid.uuid4())
    refFile = str(uuid.uuid4())
    samFile = str(uuid.uuid4())
    with open(tmpFile,'w') as fh:
      for id,(seq,count) in seqDict.items():
        fh.write(">{} \n".format(id))
        fh.write("{}\n".format(seq))

    with open(refFile,'w') as fh:
      fh.write(">REFSEQ\n")
      fh.write("{}\n".format(g_refseq))
    p = subprocess.Popen(["needle",
          "-asequence",refFile,
          "-bsequence",tmpFile,
          "-gapopen","10.0",
          "-gapextend","0.5",
          "-aformat3","sam",
          "-outfile",samFile], stdout=subprocess.PIPE)
    r = p.communicate()
    #
    # read needle output
    #
    with open(samFile,'r') as fh:
        while True:
          l = fh.readline().strip()
          if l == "": break
          l = l.split('\t')
          if l[0][0] == "@": continue
          seqID = l[0]
          cigar = l[5]
          s,count = seqDict[seqID]
          resultList.append((True,s,count,cigar))
  except Exception as e:
    print("Error running needle, EMBOSS loaded?")
  finally:
    if os.path.exists(tmpFile):
      os.remove(tmpFile)
    if os.path.exists(refFile):
      os.remove(refFile)
    if os.path.exists(samFile):
      os.remove(samFile)
  return()

#---------------------------------------------------------------------------
#   processWorker - multiprocessing worker routine to check start/end primers
#   and run needle alignment
#
#
#
#
#---------------------------------------------------------------------------

def processWorker(i):
  trimDict = defaultdict(int)
  seqDict = {}
  rl = []
  print("Worker process id for {0}: {1}".format(i, os.getpid())) 
  dataList = g_dispatchList[i] 
  #
  # check for matching start and end sequence
  # sequences that pass are added to trimDict
  # to consolidate matching sequences
  #
  seqFail = defaultdict(int)
  for k,v in dataList:
    trim = checkPrimers(k,start_primer,end_primer,refseq)
    if trim == "":
      rl.append((False,k,v,""))
    else:
      trimDict[trim] += v
  #
  # assign each sequece an ID prior to alignment with needle
  #
  seqIndex = 0
  for k,v in trimDict.items():
      seqKey = "{}_{}_{}".format(i,seqIndex,v)
      seqIndex += 1
      seqDict[seqKey] = (k,v)
  #
  # subprocess needle
  #
  execNeedle(seqDict,rl)
  return(rl)

#---------------------------------------------------------------------------
#
# read source fastq and check quality
# 
# combine identical sequences
#
# return list of passing sequences sorted by number of times seen
#
#---------------------------------------------------------------------------
def readAndCheckQuality(inputFile,stat):
  qd  = defaultdict(int)
  #
  # read fastq file, check Q, combine indentical
  #
  with open(inputFile,'r') as fh:
    while True:
      id = fh.readline().strip()
      if id == "": break
      seq = fh.readline().strip()
      p   = fh.readline().strip()
      q   = fh.readline().strip()
      fail = checkQuality(q)
      if fail >= 4: 
        stat.failQuality += 1
      else:
        qd[seq] += 1
       
      stat.totalInputSequences += 1
  #
  # sortQ = list of dictionary k,v pairs sorted by count(v)
  #
  sortQ = sorted(qd.items(),key=lambda al: al[1],reverse=True)  # largest first
  return(sortQ)

#---------------------------------------------------------------------------
#
# prepareMultiprocessing
#
#
#---------------------------------------------------------------------------
def prepareMuliprocessing(sortQ,nProc):
  batchSize = int(len(sortQ)/nProc)
  mylist = []
  index = 0
  #
  # break list of sequences to process up into nProc groups
  #
  for i in range(0,nProc):
    mylist.append(i)
    if i == (nProc-1):
      l = sortQ[index:]
    else:
      l = sortQ[index:index+batchSize]
    g_dispatchList.append(l)
    index += batchSize
  
  p = multiprocessing.Pool(nProc) 
  # map list to target function 
  if True:
    result = p.map(processWorker, mylist) 
  else:  
    result=[]
    result.append(processWorker(0))
    result.append(processWorker(1))
    result.append(processWorker(2))
    result.append(processWorker(3))
    result.append(processWorker(4))
    result.append(processWorker(5))
    result.append(processWorker(6))
    result.append(processWorker(7))
  #
  # save failures
  # combine identical good sequences
  #
  goodDict = {}
  failList = []
  for i in range(0,nProc):
    d = result[i]
    for p,seq,count,cigar in d:
      if p:
        try:
          oldCount,oldCigar = goodDict[seq]
          oldCount += count
          goodDict[seq] = (oldCount,oldCigar)
        except KeyError:
          goodDict[seq] = (count,cigar)
      else:
        failList.append((seq,count))
  return((goodDict,failList))

#---------------------------------------------------------------------------
#
# writeAlignment
#
#---------------------------------------------------------------------------
def writeAlignment(goodDict,outFile):
  sortQ = sorted(goodDict.items(),key=lambda x:x[1],reverse=True)
  with open(outFile,'w') as fh:
    for seq,(count,cigar) in sortQ:
      fh.write("{}\t{}\t{}\n".format(seq,count,cigar))

#---------------------------------------------------------------------------
#   Main Entry
#
#
# @M03100:356:000000000-C4NGM:1:1101:17289:2061 1:N:0:2
# CCCACACTATCTCAATGCAAATATCTGTCTGAAACGGTCCCTGGCTAAACTCCACCCAT...
# +
# >>>1>11C1D@DD33B3B111B1FGHBGGHDBHHHCAFGGHHHHHHHHHHHHHHHGGGH...
#
#
#
#---------------------------------------------------------------------------
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Check quality adn primer sequences of gene editing fastq files, then align to reference')
  parser.add_argument('baseDir',type=str,help='Gene Editing base directory, sub-directories are stitched_reads, aligned_reads and results')
  parser.add_argument('inputFileBase',type=str,help='name of the input assembled fastq file')
  parser.add_argument('outputDir',type=str,help='sub-directory under base')
  parser.add_argument('outputFileBase',type=str,help='name of output file')
  parser.add_argument('paramFile',type=str,help='goCongif.txt file containing primer and reference sequences as well as hdr and barcodes')
  parser.add_argument("-n", "--nproc", type=int,help="number of processes to use",default=8)
  args = parser.parse_args()
  baseDir = args.baseDir
  inputFileBase = args.inputFileBase
  outputDir = args.outputDir
  outputFileBase = args.outputFileBase
  paramFile = args.paramFile  
  g_proc = args.nproc
  if g_proc < 1: g_proc = 1
  if g_proc > 8: g_proc = 8
  print("g_proc = {}".format(g_proc))
  g_dispatchList = []
  #
  # read params
  #
  set_name,start_primer,refseq,end_primer,cut_site,donor_template,hdr,barcode,sbarcode,pam = loadParameters(paramFile)
  g_refseq = refseq
  #
  # open fastaq
  #
  inFile      = baseDir + "stitched_reads/" + inputFileBase + ".fastq"
  logFile     = baseDir + "aligned_reads/" + outputFileBase + "_ge_error.txt"
  outFile     = baseDir + "aligned_reads/" + outputFileBase + "_ge.tsv"
  #
  # statistic variables
  #
  stat = CStat()
  #
  # read input fastq file, filter on read quality, combine identical
  #
  sortQ = readAndCheckQuality(inFile,stat)
  # 
  # (organize multi-processing)
  # filter on presence of start and end primers, trim off primers
  # combine identical 
  # run needle alignment
  #
  goodDict,failList = prepareMuliprocessing(sortQ,g_proc)
  #
  # write out sorted list of sequence-alignments
  #  
  writeAlignment(goodDict,outFile)
  #
  # write out list of failed sequences
  #
  for seq,(count,cigar) in goodDict.items():
    stat.goodAlignments += count 
  logFailures(logFile,stat,failList,start_primer,refseq,end_primer)
