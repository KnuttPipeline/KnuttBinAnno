#!/usr/bin/env python

# Modified from the dbCAN Meta Server dbCAN2 Driver Script (Stand Alone Version)

import argparse
import os
from subprocess import check_output,STDOUT

parser = argparse.ArgumentParser(description='dbCAN2 HotPEP splitter')
parser.add_argument('inputFile', help='CDS input in FASTA format')
parser.add_argument('chunks',type=int,help='Number of chunks to produce (number of threads)')
parser.add_argument('dir', help='Target directory inside the  hotpep dir')

args = parser.parse_args()

count = int(check_output("tr -cd '>' < %s | wc -c" % args.inputFile, shell=True, stderr=STDOUT))
numThreads = min(args.chunks,count)
count_per_file = count / numThreads
directory = args.dir
num_files = 1
num_genes = 0
out = open(os.path.join(directory,"orfs%s.txt" % (str(num_files))), 'w')
with open(args.inputFile, 'r') as f:
   for line in f:
      if line.startswith(">"):
         num_genes += 1
         if num_genes > count_per_file and num_files != numThreads:
            out.close()
            num_files += 1
            num_genes = 0
            out.close()
            out = open(os.path.join(directory,"orfs%s.txt" % (str(num_files))), 'w')
      out.write(line)

