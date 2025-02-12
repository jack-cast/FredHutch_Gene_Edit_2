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

#
# CResult : aligned result after needle and compare to refseq
#
#    one CResult class is computed for each group of identical sequences by the
#    seq_group routine
#
#
class CResult:
  def __init__(self,k,seq,count,cigar,cgl,subs,insertions,hdr,barcode,sbarcode,pamID):
    self.k = k
    self.seq = seq
    self.count = count
    self.cigar = cigar
    self.cgl = cgl
    self.subs = subs
    self.insertions = insertions
    self.hdr = hdr
    self.barcode = barcode
    self.sbarcode = sbarcode
    self.pamID = pamID
#
#
# CStat : save stats
#
#
class CStat:
  def __init__(self):
    self.totalInputSequences = 0
    self.failQuality         = 0
    self.failPrimerCheck     = 0
    self.passPrimerCheck     = 0
    self.goodAlignments      = 0
    self.failAlignments      = 0
    self.matchRefSeq         = 0
    self.BARCODE             = 0
    self.SBARCODE            = 0
    self.indel               = 0
    self.sub                 = 0
    self.HDR                 = 0
# 
# Named tuples to store CIGAR string token parse info
#
cigToken = namedtuple('cigToken',['cmd','length','rest'])
CigEntry = namedtuple('CigEntry', ['cmd', 'length'])
#CigEntry = namedtuple('CigEntry', ['cmd', 'length'], verbose=False)
#--------------------------------------------------------------------------------
#
# parseInsertCigar - given cigar string call cigFindNextToken to break it
# into commands of form ('M', length)
#
#--------------------------------------------------------------------------------
def cigFindNextToken(cig):
    #
    # look for integers then a single (M,I,D)
    #
    r = re.compile("([0-9]+)([a-zA-Z]+)(.+)")
    m = r.match(cig)
    if m:
        return(cigToken(m.group(2),int(m.group(1)),m.group(3)))
    else:
        r = re.compile("([0-9]+)([a-zA-Z]+)")
        m = r.match(cig)
        if m:
            return(cigToken(m.group(2),int(m.group(1)),""))
        else:
            return(cigToken(False,0,""))
#--------------------------------------------------------------------------------
#
# parseInsertCigar - this routine translates a CIGAR string into a list 
# of cigEntry tuples (cmd,length)  
#
#--------------------------------------------------------------------------------
def parseInsertCigar(cig):
    #
    # 186M8I67M 
    #
    # like to return  (M 186)(I 8) (M 67)
    #
    #(186,M) (8,I) (67M)
    retList = []
    #print("Start ",cig)
    while True:
        ct = cigFindNextToken(cig)
        if (ct.cmd == False):
            return retList
        else:
            cigEntry = CigEntry(ct.cmd,ct.length)
            retList.append(cigEntry)

            cig = ct.rest
            if cig == "":
                return retList

    return False


#---------------------------------------------------------------------------
#
#  read inital parameter file
#
#
#
#---------------------------------------------------------------------------
def loadParameters(paramFile):
  if not os.path.isfile(paramFile):
    raise Exception("parameter file not found {}".format(paramFile))
  name = ""
  start_primer = ""
  refseq = ""
  end_primer = ""
  guide_seq = ""
  donor_template = []
  hdr = []
  barcode = []
  sbarcode = []
  pam = []
  with open(paramFile,'r') as fh:
    while True:
      line = fh.readline().strip()
      if line == "": break
      line = line.split(':')
      for i,j in enumerate(line):
         line[i] = line[i].strip()      
      if line[0] == "#": continue
      if line[0] == "name":
        name = line[1]
      if line[0] == "start_primer":
        start_primer = line[1]
      if line[0] == "refseq":
        refseq = line[1]
      if line[0] == "end_primer":
        end_primer = line[1]
      if line[0] == "guide_seq":
        guide_seq = line[1].strip()
      if line[0] == "donor_template":
        donor_template.append((int(line[1]),line[2]))
      if line[0] == "hdr":
        preSeq = line[1].strip()
        preLoc = refseq.find(preSeq)
        if preLoc == -1:
          raise Exception("hdr preSequence not found in reference {}".format(preSeq))
        #print(preLoc,preSeq,line[2])
        hdr.append((preLoc,preSeq,line[2]))
      #
      # barcode line has specific sequence to search for.  "." are don't care
      # 
      # barcode:.......A...CCG:TT:6:TT
      # line = ['barcode','.....A...CCG','TT',6,'TT']
      #
      if line[0] == "barcode":
        barcode.append((line[1],line[2],int(line[3]),line[4]))
      if line[0] == "sbarcode":
        sbarcode.append((int(line[1]),line[2],int(line[3]),line[4]))
      if line[0] == "pam":
        pam.append((line[1],line[2]))
  #
  # verify
  #
  if (start_primer == "") or (refseq == "") or (end_primer == "") or (guide_seq == ""):
    print(start_primer)
    print(refseq)
    print(guide_seq)
    print(end_primer)
    raise Exception("primer codes not present in parameter file")
  #
  # determine cut site from guide_seq
  #
  # cut site is first capital letter in guide_seq
  #
  cut_offset = -1
  for i in range(0,len(guide_seq)):
    if guide_seq[i].isupper():
      cut_offset = i
      break
  #
  # no longer require cut site
  #
  #if cut_offset == -1:
  #  raise Exception("no cut location for guide {}".format(guide_seq))
  guide_seq_cap = guide_seq.upper()
  guide_loc = refseq.find(guide_seq_cap)  
  if guide_loc == -1:
    raise Exception("guide seq not found in ref seq {}".format(guide_seq_cap))  
  cut_site = guide_loc + cut_offset

  return name,start_primer,refseq,end_primer,cut_site,donor_template,hdr,barcode,sbarcode,pam
