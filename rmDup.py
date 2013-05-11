#!/usr/bin/env python
# Created by Allison Carlisle on 03-26-2013
# acarlisl@soe.ucsc.edu

# Usage: rmDup.py -r reads [-o outputFile]

# takes in a fasta or fastq file and removes all identically duplicate reads
# keeps the read with the higher score
# technically this will work with any squences formatted like a fasta or fastq file
"""rmDup.py
     Usage:
          ./rmDup.py -r reads [-o output file]
     Outputs to standard out or specified output file
     Can also import and call as a function
"""
from optparse import OptionParser
import sys
import gzip
import time
import logging as log

# the wee main class that just does argument parsing and calls the rmDup function
def main(argv=None):
   parser = OptionParser()
   
   parser.add_option("-o","--output", nargs=1, metavar="OUTFILE", dest="out", help="output file of sequence")
   parser.add_option("-r","--reads", nargs=1, default=None, dest="readFiles", help="Reads (fasta or fastq)")
   (args, opts) = parser.parse_args() 
   if not args.readFiles:
    print >> sys.stderr, "usage: rmDup.py -r reads [-o outPutFile]"
    exit(1)
   
   passOut=None
   if args.out:
      passOut=args.out

   startTime=time.time()
   rmDup([args.readFiles], passOut)
   print >> sys.stderr, "total time seconds", (time.time()-startTime)

   return 0

# reads the file, removes all duplicate reads, prints the final set
# may be called as part of another python program
def rmDup(readFileList, outPutFileName=None):
   outPutFile=sys.stdout
   if outPutFileName:
       outPutFile=open(outPutFileName, 'w')

   avgReadLen=0
   totalReads=0
   tossedReads=0
   fileNumber=0
   readDict={}
   isFasta=True
   for fileIn in readFileList:
      fileNumber+=1
      f=None
      f=ParseFastQ(fileIn)
      if f: #file is good and open
         log.info("Reading "+fileIn+" "+str(time.time()))
         numberLines=0
         currentContig=""
         currentSeq=""
         currentQuality=""
         # go through lines in contig file
         for line in f:
            if line=="":
              continue
            numberLines+=1
            currentSeq=line[1]# read sequence
            currentContig=line[0][1:]# name
            currentQuality=line[3]
            avgReadLen+=len(currentSeq)
            totalReads+=1
            if numberLines == 1:
              if line[0][0] == "@":
                isFasta=False
            if currentSeq in readDict:
              if isFasta:
                tossedReads+=1
                continue
              if scoreRead(readDict[currentSeq][2]) > scoreRead(currentQuality):
                readDict[currentSeq]=[currentContig,currentSeq,currentQuality]
              else:
                tossedReads+=1
            else:
              readDict[currentSeq]=[currentContig,currentSeq,currentQuality]
   avgReadLen=float(avgReadLen)/totalReads
   print >> sys.stderr, "total reads before:", totalReads, " avg read length:", avgReadLen, " duplicates:", tossedReads, " kept:", (totalReads-tossedReads), " duplicate ratio:", float(tossedReads)/totalReads
   hChar="@"
   if isFasta:
      hChar=">"
   for key in readDict:
      print >> outPutFile, hChar+readDict[key][0]+'\n'+readDict[key][1]
      if not isFasta:
        print >> outPutFile, '+\n'+readDict[key][2]
   outPutFile.close()
   return

def scoreRead(sequence):
  score=0
  for c in list(sequence):
    score+=ord(c)
  return score

def quickWrap(text, width=60):
    return '\n'.join(text[i:i+width] for i in
                     range(0, len(text), width))

# helper fasta and fastq parsing stuff. Usually lives in a seperate file,
# but has been included for your convenience.
# from https://gist.github.com/1866279
# some changes made by me to support reading fasta, but filling in a bogus 
#    quality score value
class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
 
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
        self._fastq=True # this will be changed if when reading it's discovered to ve a fasta
        
    def __iter__(self):
        return self
    
    def randomScore(self, length):
        return "U"*length
    
    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        if self._fastq:
            for i in range(4):
                prevLine=self._file.tell()
                line = self._file.readline()
                if not line:
                    elemList.append(None)
                    continue
                if (i == 0) and (line[0]==">"):
                    self._fastq = False
                    self._hdSyms=headerSymbols=['>','+']
                    self._file.seek(prevLine)
                    break
                self._currentLineNumber += 1 ## increment file position
                if line:
                    elemList.append(line.strip('\n'))
                else: 
                    elemList.append(None)
        if not self._fastq:# this is NOT an else to pick up in case it discovers it's actually a fasta
            for i in range(2):
                line = self._file.readline()
                self._currentLineNumber += 1 ## increment file position
                if line:
                    elemList.append(line.strip('\n'))
                else: 
                    elemList.append(None)
            if elemList[1]:
                elemList.append("+")
                elemList.append(self.randomScore(len(elemList[1]))) # a fake random score
                # NOTE: this is done intentionally for another program... I know it seems weird.
            else:
                elemList.append(None)
                elemList.append(None)
            
        
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
        
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)

# End of Program   
if __name__ == "__main__":
   sys.exit(main())
    
