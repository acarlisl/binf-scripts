#!/usr/bin/env python
# Created by Allison Carlisle on 03-28-2013
# acarlisl@soe.ucsc.edu

# Usage: rmMiaDup.py -m malnFile [-o outputFile]

# takes in a maln file and removes all reads that have NUM_INPUTS < 2
"""rmDup.py
     Usage:
          ./rmMiaDup.py -m malnFile [-o output file]
     Outputs to standard out or specified output file
     Can also import and call as a function
"""
from optparse import OptionParser
import sys
import gzip
import time
import logging as log

# the wee main class that just does argument parsing and calls the rmMiaDup function
def main(argv=None):
   parser = OptionParser()
   
   parser.add_option("-o","--output", nargs=1, metavar="OUTFILE", dest="out", help="output file (wil be maln format)")
   parser.add_option("-m","--maln", nargs=1, default=None, dest="malnFile", help="Maln mia output file")
   (args, opts) = parser.parse_args() 
   if not args.malnFile:
    print >> sys.stderr, "usage: rmMiaDup.py -m malnFile [-o outPutFile]"
    exit(1)
   
   passOut=None
   if args.out:
      passOut=args.out

   startTime=time.time()
   rmMiaDup(args.malnFile, passOut)
   print >> sys.stderr, "total time seconds", (time.time()-startTime)

   return 0

# reads the file, removes all reads that do NOT have NUM_INPUTS > 1
# may be called as part of another python program
def rmMiaDup(malnFile, outPutFileName=None):
  outPutFile=sys.stdout
  if outPutFileName:
     outPutFile=open(outPutFileName, 'w')
  f=None
  f=open(malnFile)
  if f: #file is good and open
    reachSeqs=False
    read=""
    keepRead=False
    for line in f:
      if not reachSeqs:
        print >> outPutFile, line.strip('\n')
        if "__ALNSEQS__" in line:
          reachSeqs=True
      else:
        # we need to read in all the info about the read before we decide to keep or chuck it
        if line[:2]=="ID":
          if keepRead:
            print >> outPutFile, read.strip('\n')
          read=""
          keepRead=False
        if line[:10]=="NUM_INPUTS":
          if int(line[11:]) > 1:
            keepRead=True
        read+=line
    f.close()
  outPutFile.close()
  return

# End of Program   
if __name__ == "__main__":
   sys.exit(main())
    
