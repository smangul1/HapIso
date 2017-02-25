"""
Created on Mon Jan 04 18:45:18 2016

@author: Harry Yang
"""

import pysam
import sys
import numpy as np
import csv
from Bio import SeqIO
from scipy.cluster.vq import whiten
import Bio.pairwise2 as pair

from pylab import plot,show
from numpy import vstack,array
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq
import collections
from clustering import mean_shift_clustering
from collections import Counter

if len(sys.argv) < 6 :
    print "HapIso - Haplotype-specific Isoform reconstrunction"
    print "ERROR COMPARISON with LSC"
    print "Written by Harry Yang and Serghei Mangul"
    print "[1] - Gene name"
    print "[2] - String from Hercules"
    print "[3] - String from LCS"
    print "[4] - Output file"
    print "If any question: email harry2416@gmail.com"
    sys.exit(1)
name = sys.argv[1]
hercules = sys.argv[2]
lsc = sys.argv[3]
out = sys.argv[4]


alignment = pair.align.globalxx(hercules, lsc)
print "RESULT   ", name, alignment[1][2]


