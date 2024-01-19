#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 11:12:23 2024

@author: christopher
"""


import sys
import subprocess
import random
import numpy as np
import pandas as pd
import hashlib
from optparse import OptionParser

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import RNA
from anytree import AnyNode, RenderTree
from anytree.search import find
from anytree.iterators.levelorderiter import LevelOrderIter


__author__ = "Christopher Klapproth"
__institution__= "University Leipzig"
__credits__ = []
__license__ = "GPLv2"
__version__="0.1.0"
__maintainer__ = "Christopher Klapproth"
__email__ = "christopher@bioinf.uni-leipzig.de"
__status__ = "Development"


################################
#
#   Default parameters
#
################################

def_generations = 5
def_mutation_rate = 4
def_children = 4
def_selection = float(1/def_children)

nucleotides = [
    'A', 'C', 'G', 'T',
    ]

mutations_dict = {
    'A': ['C', 'G', 'T', '-'],
    'C': ['A', 'G', 'T', '-'],
    'G': ['A', 'C', 'T', '-'],
    'T': ['A', 'G', 'C', '-'],
    }


################################
#
#   File I/O functions
#
################################

def load_seq_records(filename):
    aln_handle = SeqIO.parse(handle=open(filename, 'r'), format='fasta')
    records = [record for record in aln_handle]
    return records


def load_aligned_records(filename):
    aln_handle = SeqIO.parse(handle=open(filename, 'r'), format='clustal')
    records = [record for record in aln_handle]
    return records


def write_candidate_seqs(seqs, filename):
    records = []
    count = 0

    for seq in seqs:
        count += 1
        record = SeqRecord(
        Seq(seq.replace('-', '')),
        id="candidate_%s" % count,
        name="",
        description="synthetic nucleotide sequence %s" % count,
        )
        records.append(record)
        
    SeqIO.write(records, open(filename, 'w'), format='fasta')
    return 


def write_population_file(records, filename):
    SeqIO.write(records, open(filename, 'w'), format='fasta')
    return 


def print_current_candidates(seqs, scores):
    for i in range(0, len(seqs)):
        print('#{number} \t{sequence} \t{score}'.format(
            number = i+1,
            sequence = seqs[i],
            score = scores[i],
            ))
        
        
def map_seq_to_score(aligned_records, scores):
    score_dict = {}
    
    for i in range(0, len(aligned_records)):
        seq = str(aligned_records[i].seq.ungap('-').ungap('.')).upper().replace('U', 'T')
        encode = hashlib.md5(seq.encode('utf-8')).hexdigest()
        score_dict[encode] = scores[i]
    
    return score_dict

        
################################
#
#   External program calls
#
################################

def clustalw_alignment(fasta_file):
    cmdline = "clustalw2 {population}".format(
        population = fasta_file
        )
    returnvalue = subprocess.run(cmdline.split())
    return fasta_file.replace('.fa', '.aln')


def build_cm_model(fasta_file, cm_file):
    cmdline = "cmbuild --noss -F {model} {fasta}".format(
        model = cm_file,
        fasta = fasta_file,
        )
    returnvalue = subprocess.run(cmdline.split())
    return
  

def infernal_align(fasta_file, cm_file, cm_scores, cm_alignment):
    cmdline = "cmalign --sfile {scores} --cpu 8 --outformat {form} -o {outfile} {model} {seqfile}".format(
        scores = cm_scores,
        form = 'AFA',
        outfile = cm_alignment,
        model = cm_file,
        seqfile = fasta_file,
        )
    returnvalue = subprocess.run(cmdline.split())
    
    aligned_records = [r for r in SeqIO.parse(open(cm_alignment, 'r'), format='fasta')]
    return aligned_records
        

################################
#
#   Feature calculations
#
################################

def rna_mfe(node):
    fc = RNA.fold_compound(node.seq)
    mfe_struct, mfe = fc.mfe()
    node.mfe = mfe


def CG_content(node):
    cg = (node.seq.count('C') + node.seq.count('G')) / len(node.seq)
    node.CG_content = cg


def pairwise_sequence_identity(seq1, seq2):
    count = 0
    for i in range(0, len(seq1)):
        if seq1[i] == seq2[i]:
            count += 1
    return count / len(seq1)


def mean_pairwise_identity(records):
    length = len(records)
    ids = []
    
    for i in range(0, length-1):
        for j in range(i+1, length):
            ids.append(pairwise_sequence_identity(records[i].seq, records[j].seq))
    
    return np.mean(ids)


def select_best_scoring(score_dict, fraction):
    pairs = [(str(key), score_dict.get(key)) for key in score_dict.keys()]
    sorted_pairs = sorted(pairs, key=lambda x: x[1], reverse=True)
    length = len(pairs)
    
    if (fraction <= 1):
        cutoff_index = int(length*fraction)
    elif (fraction < length):
        cutoff_index = int(fraction)
    else:
        cutoff_index = length -1
    
    return sorted_pairs[:cutoff_index]
        

def parse_cmscores(cm_scores):
    scores = []
    
    with open(cm_scores, 'r') as f:
        for line in f.readlines():
            if line.startswith('#'):
                continue
            ls = line.split()
            scores.append(float(ls[7]))
            
    return scores


################################
#
#   Tree operations
#
################################


def initialize_tree():
    return AnyNode(id='root', parent=None, generation=-1, 
                   hash=hashlib.md5(str().encode('utf-8')).hexdigest())   


def add_children_nodes(seqs, parent, generation, child_count):
    children = []
    
    for i in range(0, len(seqs)):
        name = 's_{generation}_{child}'.format(
            generation=generation,
            child=child_count+i,
            )
        seq = seqs[i]
        children.append(AnyNode(id=name, seq=seq, 
                                generation=generation,
                                hash=hashlib.md5(seq.encode('utf-8')).hexdigest()))
    parent.children = children
    return child_count+1


def get_generation_nodes(node, generation):
    iterator = LevelOrderIter(node, filter_=lambda n: n.generation == generation)
    return [node for node in iterator]


def assign_features_to_nodes(nodes, cm_score_dict):
    for node in nodes:
        rna_mfe(node)
        CG_content(node)
        node.cm_score = cm_score_dict.get(node.hash)


################################
#
#   Genetic operators
#
################################


def remove_gaps(seq):
    return seq.replace('-', '')


def point_mutate(nt):
    return random.SystemRandom().choice(mutations_dict.get(nt, nucleotides))


def point_crossover(seq1, seq2):
    index = np.random.randint(low=0, high=len(seq2))
    local_nt = seq2[index]
    return local_nt, index


def mutate(seqs, tree, generation, mutation_rate, n_children=def_children):
    return_seqs = []
    N = 0   # Total number of children added to Tree
    
    for seq in seqs:
        hash_ = hashlib.md5(remove_gaps(seq).encode('utf-8')).hexdigest()
        node = find(tree, lambda n: n.hash == hash_)
        mutated_variants = []
        children = 0
        
        # Keep generating recombined children till the target number is satisfied
        while (children < n_children):
            seq_ls = list(seq)
        
            for i in range(0, mutation_rate):
                mutation_type = random.SystemRandom().choice(['point', 'crossover'])
                
                if (mutation_type == 'point'):
                    index = np.random.randint(low=0, high=len(seq_ls))
                    local_nt = point_mutate(seq_ls[index])
                
                if (mutation_type == 'crossover'):
                    # TODO: Exclude own ident sequence
                    crossover_seq = random.SystemRandom().choice(seqs)
                    local_nt, index = point_crossover(seq_ls, crossover_seq)
                    
                seq_ls[index] = local_nt
            
            # Add the mutated sequence to the pool
            mut_seq = remove_gaps(str().join(seq_ls))
            mutated_variants.append(mut_seq)
            children += 1
            
        # Update children in Tree
        N = add_children_nodes(seqs=mutated_variants, parent=node, generation=generation, child_count=N)
        
        # Append mutated sequences to return List
        return_seqs += mutated_variants
    
    return return_seqs


################################
#
#   Main
#
################################


def main():
    usage = "\nseqgen.py [options] [alignment fasta]"
    parser = OptionParser(usage=usage, version="__version__")
    args = sys.argv
    
    if len(args) < 2:
        print(usage)
        sys.exit(1)
    
    # Parse out cmd line arguments
    parser.add_option("-o", "--output", action="store", default="", type="string", dest="output", help="Results file to write to. If empty, defaults to stdout.")
    parser.add_option("-g", "--generations", action="store", default=def_generations, type="int", dest="iterations", help="Integer defining number of cycles.")
    parser.add_option("-M", "--model", action="store", default='', type="string", dest="cov_model", help="Path to covariance model used as reference. If none is provided, will be calculated from start population.")
    parser.add_option("-m", "--mutationrate", action="store", default=def_mutation_rate, type="int", dest="mutation_rate", help="Float defining the base mutation chance per sequence and generation.")
    parser.add_option("-s", "--size", action="store", default=def_selection, type="float", dest="selection", help="Generation size to be selected for the next cycle. Can also be a fraction.")
    parser.add_option("-c", "--children", action="store", default=def_children, type="int", dest="children", help="Initial number of children nodes to generate per generation and element (before scoring).")
    parser.add_option("-l", "--label", action="store", default='RNA', type="string", dest="label", help="Data label for result table.")
    options, args = parser.parse_args()
    
    # Prepare base parameters
    aln_file = args[-1]
    iterations = options.iterations
    mutation_rate = options.mutation_rate
    n_children = options.children
    fraction = options.selection
    label = options.label
    population_file = 'population.fa'
    cm_file = 'cov_model.cm'
    candidates_file = 'candidates.fa'
    score_file = 'cm_scores.txt'
    cm_alignment = 'cm_alignment.fa'


    # TODO: Deal with U/T in input file
    # TODO: Deal with input duplicates !!

    # Write the base RNA population file
    records = load_seq_records(aln_file)
    write_population_file(records, population_file) 
    
    # Alignment of current population
    population_aln = clustalw_alignment(fasta_file=population_file)
    records = load_aligned_records(filename=population_aln)
    ungapped_seqs = [str(r.seq.ungap('-')) for r in records]
    
    # Build an initial CM model iff none is provided
    if len(options.cov_model) == 0:
        build_cm_model(fasta_file=population_file, cm_file=cm_file)
    else:
        cm_file = options.cov_model
        
    # Initialize tree
    root = initialize_tree()
    add_children_nodes(seqs=ungapped_seqs, parent=root, generation=0, child_count=0)
    
    # Score at zero generation (t=0)
    cm_aligned_records = infernal_align(fasta_file=population_file, cm_file=cm_file,
                                     cm_scores=score_file, cm_alignment=cm_alignment)
    scores = parse_cmscores(cm_scores=score_file)
    score_dict = map_seq_to_score(cm_aligned_records, scores=scores)
    
    # Add values to tree
    children = get_generation_nodes(root, generation=0)
    assign_features_to_nodes(children, cm_score_dict=score_dict)
    
    
    # Start Generation:
    for generation in range(1, iterations+1):
        print('Generation #%s' % generation)
        synthetic_population = 0
        
        # Calculate alignment of current sequence pool
        population_aln = clustalw_alignment(fasta_file=population_file)
        
        # Load current population and mutate sequences
        records = load_aligned_records(filename=population_aln)
        sequences = [str(record.seq) for record in records]
        mutated_seqs = mutate(seqs=sequences, tree=root, generation=generation,
                              mutation_rate=mutation_rate,
                              n_children=n_children)
        
        # Write candidates to Fasta:
        write_candidate_seqs(seqs=mutated_seqs, filename=candidates_file)
        
        # Align mutated entities to CM model and score them
        cm_aligned_records = infernal_align(fasta_file=candidates_file, cm_file=cm_file,
                                         cm_scores=score_file, cm_alignment=cm_alignment)
        scores = parse_cmscores(cm_scores=score_file)
        score_dict = map_seq_to_score(cm_aligned_records, scores=scores)
        children = get_generation_nodes(root, generation=generation)
        assign_features_to_nodes(children, cm_score_dict=score_dict)
        
        # Select best candidate(s) and add to current population
        best_candidates = select_best_scoring(score_dict=score_dict, fraction=fraction)
        best_hashes = [x[0] for x in best_candidates]
        records = []
        
        for child in children:
            if (child.hash not in best_hashes):
                continue
            sequence=child.seq
            synthetic_population += 1
            synthetic_record = SeqRecord(
            seq=Seq(sequence),
            id='synthetic_RNA_%s' % synthetic_population,
            name="",
            description='covariance model score: %s' % child.cm_score,
            )
            records.append(synthetic_record)
        
        # Write population file and restart the cycle
        write_population_file(records, population_file)
        records = []
    
    # Test: render tree
    print(RenderTree(root))
    
    '''
    data = {
        'Generation': generations,
        'Minimum free energy': mean_mfe,
        'Mean covariance model score': mean_cm_score,
        'Mean pairwise identity': mean_identity,
        'Mean CG content': mean_cg,
        'Label': [label]*len(generations),
        }
    df = pd.DataFrame(data=data)
    
    if len(options.output) == 0:
        print('\t'.join([c for c in df.columns]))
        for i in range(0, len(df)):
            print('\t'.join([df[c].iloc[i] for c in df.columns]))
    else:
        df.to_csv(options.output, sep='\t')
    '''  

if __name__ == '__main__':
    main()
