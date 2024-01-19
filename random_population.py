#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 14:37:47 2024

@author: christopher
"""


import random
import numpy as np
from optparse import OptionParser

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


__author__ = "Christopher Klapproth"
__institution__= "University Leipzig"
__credits__ = []
__license__ = "GPLv2"
__version__="0.5.0"
__maintainer__ = "Christopher Klapproth"
__email__ = "christopher@bioinf.uni-leipzig.de"
__status__ = "Development"


nucleotides = ['A', 'C', 'G', 'U']

def_min_length = 50
def_max_length = 100
def_n = 1000


def write_fasta(sequences, filename):
    count = 0
    records = []
    
    for seq in sequences:
        count += 1
        record = SeqRecord(
                seq=Seq(seq),
                id='random_RNA_%s' % count,
                name="",
                description=""
                # description="synthetic nucleotide sequence %s" % synthetic_population,
                )
        records.append(record)
        
    SeqIO.write(records, open(filename, 'w'), format='fasta')
    

def print_fasta(sequences):
    count = 0
    
    for seq in sequences:
        count += 1
        print('>random_RNA_%s' % count)
        print(seq)


def random_RNA(min_length, max_length):
    length = np.random.randint(low=min_length, high=max_length)
    seq = str().join([random.SystemRandom().choice(nucleotides) for i in range(0, length)])
    return seq


def main():
    usage = "\nrandom_population.py [options]"
    parser = OptionParser(usage=usage, version="__version__")
    
    # Parse out cmd line arguments
    parser.add_option("-o", "--output", action="store", default="", type="string", dest="output", help="File to write to. If empty, defaults to stdout.")
    parser.add_option("-m", "--min-length", action="store", default=def_min_length, type="int", dest="min_length", help="Minimal length for sampled sequences.")
    parser.add_option("-l", "--max-length", action="store", default=def_max_length, type="int", dest="max_length", help="Maximum length for sampled sequences.")
    parser.add_option("-n", "--num-samples", action="store", default=def_n, type="int", dest="samples", help="Number of samples to generate.")
    options, args = parser.parse_args()
    
    # Generate N random RNA sequences
    seqs = []
    
    for n in range(0, options.samples):
        seqs.append(random_RNA(min_length=options.min_length, max_length=options.max_length))
    
    if len(options.output) == 0:
        print_fasta(seqs)
    else:
        write_fasta(sequences=seqs, filename=options.output)

if __name__ == '__main__':
    main()
