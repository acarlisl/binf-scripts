#!/usr/bin/env python
# Created by Allison Carlisle on 05-01-2013
# acarlisl@soe.ucsc.edu

from optparse import OptionParser
import argparse
from matplotlib import pyplot as plt
import math
import sys
import gzip

# the wee main class that just does argument parsing and calls the function
def main(argv=None):
    """plotCoverageFromMaln.py
       Usage:
            ./plotCoverageFromMaln.py -m malnFile -o outputFileRoot [-n number_plots]
       an extension of get2xfromMaln.py
       takes in maln file(s) and outputs coverage plot(s)
       outputs a csv (maln file individual, coverage) to stdout

       automatically groups the individuals into number_plots groups
       trying to keep those with similar coverage together
    """
    #parser = OptionParser()
    #parser.add_option("-o","--output", nargs=1, metavar="OUT_NAME", dest="out", 
    #                  help="output file base for creating the plot names (OUT_NAME.INT.png)")
    #parser.add_option("-m","--maln", nargs=2, default=None, dest="maln_file",
    #                  help="Maln file(s)")
    #parser.add_option("-n","--numplots", nargs=1, default=2, dest="num_plots",
    #                  help="Integer specifying number of plots")
    #(args, opts) = parser.parse_args() 
    parser = argparse.ArgumentParser(description="""an extension of get2xfromMaln.py
       takes in maln file(s) and outputs coverage plot(s)
       outputs a csv (maln file individual, coverage) to stdout

       automatically groups the individuals into number_plots groups
       trying to keep those with similar coverage together""")
    parser.add_argument("-o","--output", nargs=1, metavar="OUT_NAME", dest="out", 
                      help="output file base for creating the plot names (OUT_NAME.INT.png)")
    parser.add_argument("-m","--maln", nargs='+', default=None, dest="maln_file",
                      help="Maln file(s)")
    parser.add_argument("-n","--numplots", nargs=1, default=2, dest="num_plots",
                      help="Integer specifying number of plots")
    args = parser.parse_args() 
    if (not args.maln_file) or (not args.out):
        print >> sys.stderr, "usage: plotCoverageFromMaln.py -m malnFile -o outputFileRoot [-n int_number_plots]"
        exit(1)  
    num_plots = 0
    try:
        num_plots = int(args.num_plots[0])
    except:
        print >> sys.stderr, "bad number of plots! Must be positive integer."
        exit(1)
    if num_plots > len(args.maln_file):
        print >> sys.stderr, "Must have at least as many maln files as plots!"
        exit(1)
    plotCoverageFromMaln(args.maln_file, num_plots, args.out[0])
    
    return 0

def plotCoverageFromMaln(maln_list, num_plots, image_stub_name):
    """Reads the file(s), sorts them into groups by coverage, plots the coverage,
            and outputs a matrix of (maln file, coverage at base pos) to std out
    """
    num_files = 0
    maln_coverage = [] # [(median_coverage,maln_file_name,[coverage_list]), ...]
    for file in maln_list:
        num_files+=1
        if file.endswith('.gz'):
            f = gzip.open(file)
        else:
            f = open(file)
        sequence = ""
        cov_count = {}
        num_lines = 0
        if f:
          reach_seqs = False
          reach_assembly = False
          span = [None,None]
          for line in f:
            num_lines+=1
            if not reach_seqs:
              if not reach_assembly:
                if line[:4] == 'SIZE':
                    reach_assembly = True
                elif "__ALNSEQS__" == line.strip():
                    reach_seqs = True
              else:
                  if line[:4] == 'GAPS':
                      sequence = list(sequence[5:])
                      reach_assembly = False
                      continue
                  sequence+=line.strip()
              
            else:
              if line[:5] == 'START':
                  span[0] = int(line.split()[1])
              if line[:3] == 'END':
                  span[1] = int(line.split()[1])
                  assert None not in span, "Bad entry in maln file line "+str(num_lines)
                  for x in range(*span):
                      if (x not in cov_count):
                          cov_count[x] = 0
                      cov_count[x]+=1
                  span = [None,None]        
          f.close()
          cov_list = cov_count.items()
          med = median(cov_list) # NOTE THAT THIS SORTS cov_list!
          cov_list = [x[1] for x in cov_list] # gets just the coverage
          print >> sys.stderr, cov_list[0:10]
          # if you switch to average you need to sort cov_count to get the 
          #    print to sys.stdout in the correct order
          print >> sys.stdout, ",".join([f.name]+[str(x) for x in cov_list])
          maln_coverage+=[(med,f.name,cov_list)]
    maln_coverage.sort()
    print >> sys.stderr, maln_coverage[0][0], maln_coverage[1][0]
    maln_coverage = groupList(maln_coverage, num_plots)
    # resume here! We need to plot!
    for i, data in enumerate(maln_coverage):
      out_img = image_stub_name+str(i)+'.png'
      plotData(data, out_img)

    return

def plotData(data, out_img):
    colors = colorList(len(data))
    # [(median_coverage,maln_file_name,[coverage_list]), ...]
    plot_args = []
    legend_args = []
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    for i, col in enumerate(colors):
        plot_args+=[range(0,len(data[i][2])),data[i][2],col]
        legend_args+=[data[i][1]]
    ax.plot(*plot_args)
    #ax.xlabel("Base Position")
    #ax.ylabel("Coverage")
    lgd = ax.legend(legend_args, loc='upper center', bbox_to_anchor=(0.5,-0.1))
    ax.grid('on')
    fig.savefig(out_img, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()
    # the memory required for a figure is not completely released until the 
    # figure is explicitly closed with close()
    return

def quickWrap(text, width=60):
    """Wraps the input text to the specified number of characters.
    """
    return '\n'.join(text[i:i+width] for i in range(0, len(text), width))

def average(values):
    """Computes the arithmetic mean of a list of numbers. 
        This function 'borrowed' from the python docs.

    >>> print average([20, 30, 70])
    40.0
    """
    return sum(values, 0.0) / len(values)

def median(values):
    """Finds the median (middle number) of a list of numbers.
    
    >>> print median([10,15,19,32,44,8,2])
    15.0
    """
    values.sort()
    middle = len(values)/2
    if len(values)%2:
        return float(values[middle][1])
    else:
        return sum([float(x[1]) for x in values[middle-1:middle+1]], 0.0)/2

def groupList(input_list, num_groups):
    """Splits the list into num_groups groups. Assumes the input 
        list is sorted. Note that num_groups must be <= len(list).
        Adds to solution from 
        http://stackoverflow.com/questions/1624883/alternative-way-to-split-a-list-into-groups-of-n/1625023#1625023

    >>> print groupList([1,2,3,4,5,6,7], 3)
    [(1,2,3),(4,5,6),(7,)]
    """
    from itertools import izip_longest
    assert num_groups <= len(input_list), "Must have less than %d groups" % len(input_list)
    out_list = list(izip_longest(*(iter(input_list),) * (len(input_list)/num_groups) ))
    out_list[-1] = tuple(x for x in out_list[-1] if x)
    return out_list

def colorList(num_colors, depth=1):
    """Returns a list of HSV colors maximally spaced around the color wheel.
      Ends up looking like a rainbow for many num_colors.
      The depth parameter moves each hue through saturation and value:
          eg. [(pink,red,dark_red), (light_yellow, yellow, dark_yellow)...]
      If not specified, you just get the color (eg, [red, yellow, ...])

      This is adapted from a javaScript function I wrote for another program.

      >>> colorList(3,2)
      [('#660000', '#ffc6c6'), ('#006600', '#c6ffc6'), ('#000066', '#c6c6ff')]

      >>> colorList(2)
      ['#ff0000', '#00ffff']
    """
    min_val = 0.4 # value less than ~ 40% gets hard to see
    min_sat = 0.2 # saturation less than ~ 20% gets hard to see
    col_jump = 1.0/num_colors
    val_jump = 0
    hue = 0.0
    satFun = lambda x:(-1.3*x)+1.52 # plug in value to get saturation
    # this formula set by my personal preference of how colors should move
    # through saturation and value
    out=[]
    if depth == 1:
        for i in range(0,num_colors):
            out+=[hsvToHex(hue,1,1)]
            hue+=col_jump
    else:
        val_jump = (1.0-min_val)/(depth-1)
        
        for i in range(0,num_colors):
            value = min_val
            colors = []
            for j in range(0,depth):
                 colors+=[hsvToHex(hue,satFun(value),value)]
                 value+=val_jump
            out+=[tuple(colors)]
            hue+=col_jump
    return out

def hsvToHex(h, s, v):
    """ Takes in hue, saturation, value on a 0-1 scale and returns 
         hex color (eg #a4b26b).
         Used hsv to rgb conversion from:
         http://stackoverflow.com/questions/2353211/hsl-to-rgb-color-conversion
    """
    i = math.floor(h * 6)
    f = h * 6 - i
    p = v * (1 - s)
    q = v * (1 - s * f)
    t = v * (1 - s * (1 - f))
    options = {0: lambda : rgbToHex(int(v*255), int(t*255),int(p*255)),
               1: lambda : rgbToHex(int(q*255), int(v*255),int(p*255)),
               2: lambda : rgbToHex(int(p*255), int(v*255),int(t*255)),
               3: lambda : rgbToHex(int(p*255), int(q*255),int(v*255)),
               4: lambda : rgbToHex(int(t*255), int(p*255),int(v*255)),
               5: lambda : rgbToHex(int(v*255), int(p*255),int(q*255))}
    return options[i]()

def rgbToHex(r, g, b):
    """Takes in red, green, blue component on a 0-255 scale and Returns
        hex color (eg #a4b26b).
    """
    return '#'+''.join([str(hex(x)[2:]).rjust(2,'0') for x in [r,g,b]])





# End of Program   
if __name__ == "__main__":
   sys.exit(main())
    
