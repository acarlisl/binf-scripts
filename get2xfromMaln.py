#!/usr/bin/env python
# Created by Allison Carlisle on 03-28-2013
# acarlisl@soe.ucsc.edu

# Usage: get2xfromMaln.py -m malnFile [-o outputFile]

# takes in a maln file and outputs a consensus with poor supported bases as N
# default 2
"""rmDup.py
     Usage:
          ./get2xfromMaln.py -m malnFile [-o output file]
     Outputs to standard out or specified output file
     Can also import and call as a function
"""
from optparse import OptionParser
import sys
import gzip
import subprocess
import tempfile

# the wee main class that just does argument parsing and calls the rmMiaDup function
def main(argv=None):
   parser = OptionParser()
   
   parser.add_option("-o","--output", nargs=1, metavar="OUTFILE", dest="out", help="output file (wil be maln format)")
   parser.add_option("-m","--maln", nargs=1, default=None, dest="malnFile", help="Maln mia output file")
   parser.add_option("-s","--support", nargs=1, default=2, dest="support", help="Integer specifying support")
   (args, opts) = parser.parse_args() 
   if not args.malnFile:
    print >> sys.stderr, "usage: rmMiaDup.py -m malnFile [-o outPutFile] [-s INT_DESIRED_SUPPORT]"
    exit(1)
   
   passOut=None
   support=0
   if args.out:
      passOut=args.out
   try:
       support=int(args.support)
   except:
       print >> sys.stderr, "bad support number! Must be positive integer."
       return -1
   getMinCov(args.malnFile, support, passOut)

   return 0

# reads the file, outputs consensus with at least the minimum number
# of supporting reads
# may be called as part of another python program
def getMinCov(malnFile, support=2, outPutFileName=None):
  outPutFile=sys.stdout
  if outPutFileName:
     outPutFile=open(outPutFileName, 'w')
  f41File =  tempfile.NamedTemporaryFile(delete=False)
  args = ["ma", "-M", malnFile, "-f", str(41)]
  maP = subprocess.Popen(args,bufsize=-1, stdout=f41File, stderr=subprocess.PIPE )
  #maP.wait()
  while maP.returncode == None:
    maP.poll()
    if "sucka" in subprocess.PIPE:
      print >> sys.stderr, "ma out of memory", malnFile
      exit(1)
  if (maP.returncode != 0):
    print >> sys.stderr, "ma exited badly", malnFile
    exit(1)
  sequence=""
  consSequence=""
  lineNum=0
  avgCovAgree=0
  avgCov=0
  with open(f41File.name) as f:
    for line in f:
      lineNum+=1
      line=line.strip()
      line=line.split()
      base=line[2]
      if base=="-":
        continue
      f41PosKey={'A':4,'C':5,'G':6,'T':7,'N':3}
      covAgree=int(line[f41PosKey[base]])
      avgCovAgree+=covAgree
      consSequence+=base
      avgCov+=int(line[3])
      if covAgree >= support:
         sequence+=base
      else:
         sequence+='N'
    avgCov=avgCov/len(sequence)
  tmpName=malnFile.split("/")[-1] # get rid of the path
  tmpName=tmpName.split(".")[0] # get rid of dots
  print >> outPutFile, ">"+tmpName+"_cov_"+str(avgCov)
  print >> outPutFile, quickWrap("".join(sequence))
  outPutFile.close()
  return

def quickWrap(text, width=60):
    return '\n'.join(text[i:i+width] for i in
                     range(0, len(text), width))

# End of Program   
if __name__ == "__main__":
   sys.exit(main())
    
