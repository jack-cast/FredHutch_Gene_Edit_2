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
from geneEditCommon import *
##from Bio import SeqIO
##alignments = pairwise2.align.globalxs(seq1.seq, seq2.seq,-4,-1)

#------------------------------------------------------------------------------------------------------------
# summarize
#   go through all sequences and summarize number of WT, sub, indel, HDR, barcode, sbarcode
#
#------------------------------------------------------------------------------------------------------------
def summarize(finalGroup,refseq,statt,outputFile,bcFile):
  #
  #
  # indel rate: this is the rate of all indels INCLUDING NHEJ, HDR and BARCODE
  # sub rate: substitutions that occur on non-indel sequences  (have WildType CIGAR)
  #
  #
  ct = statt.goodAlignments   # total number of good alignments
  g_total    = 0
  g_wt       = 0
  g_sbarcode = 0
  g_barcode  = 0
  g_hdr      = 0
  g_indel    = 0
  g_sub      = 0
  ssH = defaultdict(int)
  ssB = defaultdict(int)
  pamInBarcode = defaultdict(int)
  #
  # barcodes stored based on pamID
  #
  #ssD = defaultdict(int)
  barcodeDict = {}
  #
  # classify each sequence(group)
  #
  for fg in finalGroup:
    g_total += fg.count
    #
    if fg.sbarcode != "":
      ssB[fg.sbarcode] += fg.count
      stat.SBARCODE    += fg.count
      g_sbarcode       += fg.count
    elif fg.barcode != "":
      try:
        bcDefDict = barcodeDict[fg.pamID]
        bcDefDict[fg.barcode] += fg.count
      except KeyError:
        bcDefDict = defaultdict(int)
        bcDefDict[fg.barcode] += fg.count
        barcodeDict[fg.pamID] = bcDefDict
      g_barcode       += fg.count
      #ssD[fg.barcode] += fg.count
      statt.BARCODE   += fg.count
      # for barcode only, check PAM
      if fg.pamID != "":
        pamInBarcode[fg.pamID] += fg.count
      else:
        pamInBarcode["PamNone"] += fg.count
    elif fg.hdr != "":
      g_hdr += fg.count
      ssH[fg.hdr] += fg.count
      stat.HDR += fg.count
    elif len(fg.cgl) > 1:
      g_indel += fg.count
    elif fg.subs != []:
      g_sub += fg.count
    else:
      g_wt += fg.count
  statt.sub = g_sub
  statt.indel = g_indel
  print("SUMMARIZE count = {}, {}".format(g_total,outputFile))
  if g_total == 0:
    return
  #
  # these are global
  #
  fRefseq   = float(g_wt) / float(g_total)
  fIndel    = float(g_indel) / float(g_total)
  fSub      = float(g_sub) / float(g_total)
  fBARCODE  = float(g_barcode) / float(g_total)
  fSBARCODE = float(g_sbarcode) / float(g_total)
  fHDR      = float(g_hdr) / float(g_total)
  fBaseErrorRate = (float(g_total-g_wt)/float(g_total)) / float(len(refseq))
  sortPam = sorted(pamInBarcode.items(),key=lambda x:x[1],reverse=True)
  with open(outputFile,'w') as fh:
    fh.write("GoodAlign,WildType,LengthRef,Indel,Substitution,HDR,Barcode,failAlign")
    for k,v in sortPam:
      fh.write(",{}".format(k))
    fh.write("\n")
    fh.write("{0:9d},".format(statt.goodAlignments))
    fh.write("{0:9d},".format(g_wt))
    fh.write("{0:9d},".format(len(refseq)))
    fh.write("{0:9d},".format(g_indel))
    fh.write("{0:9d},".format(g_sub))
    fh.write("{0:9d},".format(g_hdr))
    fh.write("{0:9d},".format(g_barcode))
    fh.write("{0:9d}".format(statt.failAlignments))
    for k,v in sortPam:
      fh.write(",{}".format(v))
    #fh.write("{0:9d}".format(g_sbarcode))
    fh.write("\n")

    fh.write("# good alignment             = {0:9d}\n".format(statt.goodAlignments))
    fh.write("# match refseq               = {0:9d}\n".format(g_wt))
    fh.write("# number of indel            = {0:9d}\n".format(g_indel))
    fh.write("# number of substitution     = {0:9d}\n".format(g_sub))
    fh.write("# number of HDR seq          = {0:9d}\n".format(g_hdr))
    fh.write("# number of BARCODE seq      = {0:9d}\n".format(g_barcode))
    for k,v in sortPam:
      fh.write("#    pam ID {0:s}             = {1:9d}\n".format(k,v))
    fh.write("#--- % of total sequences that pass Quality, primer and alignment checks\n")
    fh.write("# freguency of good sequence = {0:9.2f}\n".format(float(statt.goodAlignments)/float(g_total)))
    fh.write("#--- % of sequences with good alignent that:\n")
    fh.write("# match refseq               = {0:9.2f}\n".format(fRefseq))
    fh.write("# mock per base error rate   = {0:11.4f}\n".format(fBaseErrorRate))
    fh.write("# indel                      = {0:9.4f}\n".format(fIndel))
    fh.write("# substitution               = {0:9.4f}\n".format(fSub))
    #fh.write("#--- % of sequences with good alignment that have HDR\n")
    fh.write("# freguency of HDR           = {0:9.4f}\n".format(fHDR))
    #fh.write("#--- % of sequences with good alignment that have BARCODE\n")
    fh.write("# freguency of BARCODE       = {0:9.4f}\n".format(fBARCODE))
    #fh.write("#--- % of sequences with good alignment that have SBARCODE\n")
    fh.write("# freguency of SBARCODE      = {0:9.4f}\n".format(fSBARCODE))
  #
  # write barcode and sbarcode
  #
  with open(bcFile,'w') as fh:
    fh.write("#Fixed base delimited barcodes\n")
    for pamID,dd in barcodeDict.items():
      ssSort = sorted(dd.items(),key = lambda x: x[1],reverse=True)
      for k,v in ssSort:
        fh.write("{}     {},{}\n".format(pamID,k,v))
    fh.write("#Substitution-only barcodes\n")
    ssSort = sorted(ssB.items(),key = lambda x: x[1],reverse=True)
    for k,v in ssSort:
      fh.write("{},{}\n".format(k,v))
#------------------------------------------------------------------------------------------------------------
# alignment:
#    
#  format sequence data for output
#
#  To print the reference, alignment and indels
#  TCCCTTCCTCTTTTCTGCTCACACAGGAAGCCCTGGAAGCTGCTTCCTCAGACATGCCGCTGCTGCTACTGCTGCCCCTGCTTCCCTTCCGTGAGTGGCTGTGG
#  ||||||||||||*||||||||||||||||||||*||||||**||*|||||||||||||*|||||||||||||||||||||||        ||||*||||||*||
#  TCCCTTCCTCTTcTCTGCTCACACAGGAAGCCCgGGAAGCctCTgCCTCAGACATGCCaCTGCTGCTACTGCTGCCCCTGCT--------GTGAaTGGCTGcGG
#  
#------------------------------------------------------------------------------------------------------------

def alignment(csvF,txtfile,finalGroup,refseq,statt):
    ct = statt.goodAlignments   # total number of good alignments
    with open(csvF, 'w') as xl:
        xl = csv.writer(xl, delimiter=",", quotechar = '"', quoting=csv.QUOTE_MINIMAL)
        xl.writerow(["Sequence","Count","Frequency","Match/Deletions","Substitution","Insertion","RefSeqL","SSODN"])
        g_total = 0
        for fg in finalGroup:
          g_total += fg.count


        assert g_total == ct


        with open(txtfile,'w') as fh:
            for fg in finalGroup:
                #(k,seq,count,cigar,cgl,subs,insertions,hdr,barcode) 
                #
                # Append subs to cigar
                #
                fh.write(">{}\n".format(fg.k))
                fh.write("{}\n".format(fg.count))
                fh.write(str(float(fg.count)/float(ct))+"\n")
                if fg.hdr != "":
                  fh.write("HDR: {}\n".format(fg.hdr))
                if fg.barcode != "":
                  fh.write("BARCODE: {}\n".format(fg.barcode))
                if fg.sbarcode != "":
                  fh.write("SBARCODE: {}\n".format(fg.sbarcode))
                if fg.pamID != "":
                  fh.write("PAM: {}\n".format(fg.pamID))
                freq = (str(float(fg.count)/float(ct)))
                #totalFreq = (str(float(fg.count)/float(g_total)))
                #print(count,numAllSeq,totalFreq)
                # use space to separate individual subs, use : to separate
                # mutations
                subGen = (" ".join(map(str,sg)) for sg in fg.subs)
                sub_list = ":".join(subGen) 
                #sub_list = ":".join(map(str, fg.subs))
                ins_list = "".join(fg.insertions)
                #
                # Trimming the seq
                # 
                # To print the reference according to the alignment.
                #
                fh.write("ref:")
                i = 0
                j = 0
                for m in fg.cgl:
                    cmdLength = int(m[1])
                    if m[0] == "M" or m[0] == "D":
                      for t in range(0,cmdLength):
                        fh.write(refseq[j])
                        j += 1
                    elif m[0] == 'I':
                        space = int(m[1])
                        for s in range(0,space):
                            fh.write(" ")

                fh.write("\n")  

                #
                # To print alignment: Match = |, D = -, Sub = *
                #
                fh.write("    ")
                i = 0
                j = 0
                for m in fg.cgl:
                    cmdLength = int(m[1])
                    if m[0] == "M":
                      for t in range(0,cmdLength):
                        if fg.seq[i] != refseq[j]:
                          fh.write("*")
                        else:
                            fh.write("|")
                        i += 1
                        j += 1
                    elif m[0] == "D":
                        #cmdLength = int(m[1])
                        j = j + cmdLength
                        for u in range(0, cmdLength):
                            fh.write(" ")
                    elif m[0] == 'I':
                        space = int(m[1])
                        for s in range(0,space):
                            fh.write(" ")
                        i = i + cmdLength

                fh.write("\n")  

                #
                # To print alignment: D = -, Substitution = lower case
                #
                fh.write("ts :")
                i = 0
                j = 0
                al_seq = ""
                for m in fg.cgl:
                    cmdLength = int(m[1])
                    if m[0] == "M":
                      for t in range(0,cmdLength):
                        if fg.seq[i] != refseq[j]:
                          fh.write(fg.seq[i].lower())
                          al_seq = al_seq + fg.seq[i].lower()
                        else:
                            fh.write(fg.seq[i])
                            al_seq = al_seq + fg.seq[i]
                        i += 1
                        j += 1
                    elif m[0] == "D":
                        #cmdLength = int(m[1])
                        j = j + cmdLength
                        for u in range(0, cmdLength):
                            fh.write('-')
                            al_seq = al_seq + "-"
                    elif m[0] == 'I':
                      insertSeq = fg.seq[i:i+cmdLength]
                      fh.write(insertSeq)
                      al_seq = al_seq + insertSeq
                      i = i + cmdLength

                fh.write("\n")

                xl.writerow([al_seq,fg.count,freq,fg.cigar,sub_list,ins_list,len(refseq),fg.hdr,fg.barcode])
                                
                #
                # indel cigar
                #
                fh.write("INDEL,")
                first = True
                for m in fg.cgl:
                  if first:
                    fh.write("{},{}".format(m[0],m[1]))
                  else:
                    fh.write(",{},{}".format(m[0],m[1]))
                  first = False
                fh.write("\n")
                #
                # substitutions
                #
                fh.write("SUB,")
                first = True
                for s in fg.subs:
                  if first:
                    fh.write(str(s[0]) + "," + s[1] + "," + s[3])
                  else:
                    fh.write("," + str(s[0]) + "," + s[1] + "," + s[3])
                  first = False
                fh.write("\n")

    print("Results written in the files!")
    return

#-------------------------------------------------------------------------------------------------
# findFixedHDR
# seq_group helpfer function to search for a fixed hdr pattern at a given location 
#
# hdr pattern specification is startLoction:preSequence:hdrSequence
#
#  ie: ACGT      TTCCCTT
#      (pre)     (hdr)
#
#   allow one mutation in hdr seq?
#
#-------------------------------------------------------------------------------------------------
def findFixedHDR(seq,hdrPatterns,allowMutation):
  hdr = ""
  for (nStart,pPre,pHDR) in hdrPatterns:
    p1 = nStart
    p2 = p1 + len(pPre)
    p3 = p2 + len(pHDR)
    template = seq[p2:p3]
    # does pre-sequence match?
    if pPre == seq[p1:p2]:
      # hdr exact match
      if (pHDR == template):
        return(template)
      if not allowMutation:return("")
      if len(pHDR) != len(template): return("")
      #
      # allow one mutation
      #
      for i in range(0,len(pHDR)):
        l1 = list(template)
        l1[i] = 'N'
        l2 = list(pHDR)
        l2[i] = 'N'
        if l1 == l2:
          return(template)
  return("")

#-------------------------------------------------------------------------------------------------
# findBarcode
# seq_group helpfer function to search hdr barcode with fixed bases before and after code 
#
# barcode pattern specification is startLoction:preSequence:barcodeLength:postSeq
#
#  ie: ...CC..AC    TT        ACACACAC   TT
#      (pre)        (fixed)   (barcode)  (post)
#
#-------------------------------------------------------------------------------------------------
def findBarcode(seq,barcodePatterns):
  #
  # does this seq contain BARCOEDE defined as:
  #   [start location, fixed matching template, variable bases, fixed ending template
  #   51,TCCTCCTGACTT,6,TT
  #
  #
  barcode = ""
  for (preSeq,preFix,nCode,postFix) in barcodePatterns:
    # check pre-sequence
    if len(seq) < (len(preSeq) + len(preFix) + nCode + len(postFix)):return "" 
    for loc in range(0,len(preSeq)):
      if preSeq[loc] != '.':
        if preSeq[loc] != seq[loc]:return("")
    # check prefix
    p1 = len(preSeq)
    p2 = p1 + len(preFix)
    p3 = p2 + nCode
    p4 = p3 + len(postFix)
    s1 = seq[p1:p2]
    s2 = seq[p2:p3]
    s3 = seq[p3:p4]

    if s1 == preFix and s3 == postFix:
      #print("found Barcode: preSeq  = {}".format(preSeq))
      #print("found Barcode: preFix  = {}".format(s1))
      #print("found Barcode: barcode = {}".format(s2))
      #print("found Barcode: postFix = {}".format(s3))
      #print("                         {}".format(seq))
      return(s2)
  
  return("")

#-------------------------------------------------------------------------------------------------
# findSubBarcode
# seq_group helpfer function to search hdr barcode with only substitution changes 
#
# barcode pattern specification is startLoction:preSeq:barcodeLength:postSeq
#
#  ie: TT      ACACACAC   TT
#      (pre)   (barcode)  (post)
#
#-------------------------------------------------------------------------------------------------
def findSubBarcode(seq,sbarcodePatterns):
  sbarcode = ""
  for (nStart,pPre,nLength,pPost) in sbarcodePatterns:
    p1 = nStart
    p2 = nStart + len(pPre)
    p3 = p2 + nLength
    p4 = p3 + len(pPost)

    if seq[p1:p2] != pPre: continue
    if seq[p3:p4] != pPost: continue
    sbRef = refseq[p2:p3]
    sbSeq = seq[p2:p3]
    # if sequences are the same, no barcode
    if sbSeq == sbRef: continue
    # sequence must differ from wild type at every base to be a barcode (allow 1 mut?)
    match = False
    # different for every base compared to refseq      
    for i in range(0,len(sbSeq)):
      if sbSeq[i] == sbRef[i]:
        match = True
        break
    if match == False:
       return(sbSeq)

  return("")

#-------------------------------------------------------------------------------------------------
# findPamID
#
#     pam = [[xxxxxxxxxx,name]]
#
#       where the first sequence is a pam ID and the second param is a name
#
#
#-------------------------------------------------------------------------------------------------
def findPamID(seq,pam):
  for (pamRef,pamID) in pam:
    i = seq.find(pamRef)
    if i != -1:
      return(pamID)
  return("")   
#-------------------------------------------------------------------------------------------------
# seq_group
#   Identify substitutions based on alignment
#   Create sequence groups by combining information for each seqeunce.
#
#  seqList  (sequnce,n,cigar)
#  refseq             reference DNA sequence
#  hdrPatterns        optional fixed HDR insert
#          - list of (startLocation, sequence)
#  barcodePatterns    optianal barcode HDR
#          - list of (startLocation, sequence, barcodeLen, stopSeq)
#
#-------------------------------------------------------------------------------------------------
def seq_group(seqList,refseq,hdrPatterns,barcodePatterns,sbarcodePatterns,pam):
    #
    #
    #
    sortD = []
    badAlignmentList = []
    for seq,count,cigar in seqList:
        cgList = parseInsertCigar(cigar)
        #if len(cgList) <= 20:
        if len(cgList) <= 4:
          sortD.append((seq,count,cigar,cgList))
        else:
          badAlignmentList.append((cigar,count))
          print("Bad alignment for sequence")
          print(seq)
          print(cigar)
    #
    # finalGroups will be the final list of sequence groups
    #
    finalGroups =[]
    #
    # for each sequence group do a base-by-base comparison to refseq 
    #
    total_good_count = 0
    for seq,count,cigar,cgl in sortD:
      #
      # substitution mutations and insertion sequences
      #
      subs = []
      insertions = []
      align_seq = []
      #
      # base by base comparison
      #    i = test sequence index
      #    j = reference sequence index
      #    cgl = list of (command,length) from CIGAR  M = match, I = insert, D = deletion
      #
      i = 0
      j = 0
 
      for m in cgl:
        cmdLength = int(m[1])
        if m[0] == "M":
          for t in range(0,cmdLength):
            #print("comp {} {}".format(i,j))
            if seq[i] != refseq[j]:
              #print("substitue at ",i,seq[i],j,refseq[j])
              mut = (j+1,refseq[j],i+1,seq[i])  # Increase the position no. by 1 (i and j) because the count in python begins at 0.
              subs.append(mut)
              #print(subs)
            i += 1
            j += 1
        elif m[0] == "D":
          j = j + cmdLength
        elif m[0] == 'I':
          insertSeq = seq[i:i+cmdLength]
          insertions.append(insertSeq)
          i = i + cmdLength

      total_good_count += count       
      #
      # does this sequence contain a constant HDR?
      #
      hdr = findFixedHDR(seq,hdrPatterns,True)
      #
      # does this sequence contain a barcode hdr with fixed ID bases
      #
      barcode = findBarcode(seq,barcodePatterns)
      #
      # another type is sbarcode for (substitution-only barcode) 
      # no fixed ID bases
      #
      sbarcode = findSubBarcode(seq,sbarcodePatterns)
      #
      # for barcode sequences, split into goups based in pam mutations
      #
      pamID = findPamID(seq,pam)
      #
      # combine for final output
      #
      fg = CResult('-',seq,count,cigar,cgl,subs,insertions,hdr,barcode,sbarcode,pamID)

      finalGroups.append(fg)

    return [finalGroups, total_good_count, badAlignmentList]



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
#
# useage: python filterFasta.py file.fastq out.fastq
#
#
if __name__ == "__main__": 
  if len(sys.argv) != 6:
    print("python geOutput.py baseDir inputFileBase outputDir outputFileBase paramFile")
    exit(0)
  script,baseDir,inputFileBase,outputDir,outputFileBase,paramFile = sys.argv
  #
  # read params
  #
  set_name,start_primer,refseq,end_primer,cut_site,donor_template,hdr,barcode,sbarcode,pam = loadParameters(paramFile)
  #
  # open fastaq
  #
  inFile      = baseDir + "aligned_reads/" + outputFileBase + "_ge.tsv"
  csvFile     = baseDir + outputDir + outputFileBase + "_ge.csv"
  bcFile      = baseDir + outputDir + outputFileBase + "_bc.txt"
  mutFile     = baseDir + outputDir + outputFileBase + "_ge.txt"
  statFile    = baseDir + outputDir + outputFileBase + "_stat.csv"
  #
  # statistic variables
  #
  stat = CStat()
  #
  # read input from aligned
  #
  dl = [] 
  with open(inFile,'r') as fh:
    while True:
      seq = fh.readline().strip()
      if seq == "": break
      seq,n,cigar = seq.split('\t')
      dl.append((seq,int(n),cigar))
  #
  # identify indel and substitutions
  # use needle alignments to identify indels, also compare to refseq to identify subsitutions
  # identify hdr and/or barcode hdr
  #
  finalgroup,stat.goodAlignments,failAlign = seq_group(dl,refseq,hdr,barcode,sbarcode,pam)
  #
  # add up alignment fails
  #
  stat.failAlignments = 0
  for cigar,count in failAlign:
    stat.failAlignments += count 
  #
  # how many seq match refseq
  #
  for fg in finalgroup:
    if fg.seq == refseq:
      stat.matchRefSeq += fg.count
  if stat.goodAlignments > 0:
    mockErrorRate = (float(stat.goodAlignments)-float(stat.matchRefSeq))/float(stat.goodAlignments)
  else:
    mockErrorRate = 0.0
  #
  # create output files
  #
  summarize(finalgroup,refseq,stat,statFile,bcFile)
  out = alignment(csvFile,mutFile,finalgroup,refseq,stat)
