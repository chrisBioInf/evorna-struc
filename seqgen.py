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
import hashlib
from optparse import OptionParser

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import RNA
from anytree import AnyNode, RenderTree
from anytree.search import find
from anytree.iterators.levelorderiter import LevelOrderIter
from anytree.exporter import JsonExporter


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

def_generations = 10
def_mutation_rate = 1
def_children = 2
def_selection = 0.25

protected_start = 33
protected_end = 38
protected_range = set([x for x in range(protected_start, protected_end)])


nucleotides = [
    'A', 'C', 'G', 'U',
    ]

mutations_dict = {
    'A': ['C', 'G', 'U'],
    'C': ['A', 'G', 'U'],
    'G': ['A', 'C', 'U'],
    'U': ['A', 'G', 'U'],
    }

mutation_types = ['point', 'deletion', 'insertion',]

default_weights = [0.5, 0.25, 0.25]  # Should sum up to 1


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


def join_to_records(seqs, ids):
    records = []
    
    for i in range(0, len(seqs)):
        seq = seqs[i]
        record = SeqRecord(
        Seq(seq.replace('-', '')),
        id=ids[i],
        name="",
        description="",
        )
        records.append(record)
        
    return records


def write_candidate_seqs(records, filename):       
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
        name = aligned_records[i].id
        score_dict[name] = scores[i]
    
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
    cmdline = "cmalign --nonbanded --sfile {scores} --cpu 8 --outformat {form} -o {outfile} {model} {seqfile}".format(
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
#   Feature calculation
#
################################

def rna_mfe(node):
    fc = RNA.fold_compound(node.seq)
    mfe_struct, mfe = fc.mfe()
    return round(mfe, 2)


def CG_content(node):
    cg = (node.seq.count('C') + node.seq.count('G')) / len(node.seq)
    return round(cg, 2)


def hamming_distance(s1, s2):
    length = len(s1)
    sum_of_distances = sum(c1 != c2 for c1, c2 in zip(s1, s2))
    return (sum_of_distances / length)


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


def initialize_tree(options, weights):
    return AnyNode(id='root', parent=None, generation=-1, 
                   mutations=options.mutation_rate,
                   n_children=options.children,
                   selection=options.selection,
                   crossover_length=options.crossover_length,
                   weights=';'.join([str(w) for w in weights]),
                   hash=hashlib.md5(str().encode('utf-8')).hexdigest())   


def add_children_nodes(seqs, parent, generation, child_count, mutations=[]):
    children = []
    node_ids = []
    
    for i in range(0, len(seqs)):
        child_count += 1
        name = 's_{generation}_{child}'.format(
            generation=generation,
            child=child_count,
            )
        node_ids.append(name)
        seq = seqs[i]
        new_node = AnyNode(id=name, seq=seq, 
                                generation=generation,
                                hash=hashlib.md5(seq.encode('utf-8')).hexdigest())
        if (len(mutations) > i):
            new_node.mutation = ' '.join(mutations[i])
        children.append(new_node)
        
    parent.children = children
    return node_ids, child_count


def get_generation_nodes(node, generation):
    iterator = LevelOrderIter(node, filter_=lambda n: n.generation == generation)
    return [node for node in iterator]


def update_tree(records, parent_nodes):
    node_dict = {}
    children = []
    count = 0
    print(records)
    
    for node in parent_nodes:
        node_dict[node.id] = node
    
    for record in records:
        count += 1
        parent_id = record.id.split('-')[0]
        parent_node = node_dict.get(parent_id)
        seq = str(record.seq)
        generation = parent_node.generation +1
        name = 's_%s_%s' % (generation, count)
        new_node = AnyNode(id=name, seq=seq, 
                           generation=generation,
                           parent=parent_node,
                           hash=hashlib.md5(seq.encode('utf-8')).hexdigest())
        children.append(new_node)
    
    return children


def assign_similarity_scores(nodes, aligned_records):
    hamming_dict = {}
    
    for i in range(0, len(aligned_records)):
        distances = []
        for j in range(0, len(aligned_records)):
            if (i == j):
                continue
            distances.append(hamming_distance(aligned_records[i].seq, aligned_records[j].seq))
        
        hamming_dict[aligned_records[i].id] = np.mean(distances)
    
    for node in nodes:
        hamming_dist = hamming_dict.get(node.id, 0)
        node.mean_hamming_distance = hamming_dist


def assign_features_to_nodes(nodes, cm_score_dict):
    for node in nodes:
        node.mfe = rna_mfe(node)
        node.CG_content = CG_content(node)
        node.cm_score = round(cm_score_dict.get(node.id), 2)


################################
#
#   Genetic operators
#
################################


def remove_gaps(seq):
    return seq.replace('-', '')


def point_mutate(nt):
    return random.SystemRandom().choice(mutations_dict.get(nt, nucleotides))


def insertion():
    return random.SystemRandom().choice(nucleotides)


def crossover(records, crossover_length=4):
    to_mutate = []
    
    for record in records:
        shuffle_record = SeqRecord(
        record.seq,
        id=record.id,
        name="",
        description="",
        )
        to_mutate.append(shuffle_record)
    
    np.random.shuffle(to_mutate)
    crossover_pairs = [(to_mutate[i-1], to_mutate[i]) for i in range(1, len(to_mutate), 2)]
    crossed_records = []
    
    for (record_a, record_b) in crossover_pairs:
        ls_seq_a = [nt for nt in record_a.seq]
        ls_seq_b = [nt for nt in record_b.seq]
        
        start_index = np.random.randint(low=0, high=min(len(ls_seq_a)-(crossover_length+1),
                                                        len(ls_seq_b)-(crossover_length+1)))
        attempts = 0
        
        while (start_index in protected_range) or (start_index+crossover_length in protected_range):
            attempts += 1
            start_index = np.random.randint(low=0, high=min(len(ls_seq_a)-(crossover_length+1),
                                                            len(ls_seq_b)-(crossover_length+1)))
            
            if  (attempts < 32):
                break
        
        stop_index = start_index +crossover_length
        local_nt_a = ls_seq_a[start_index:stop_index]
        local_nt_b = ls_seq_b[start_index:stop_index]
        
        for i in range(0, crossover_length):
            ls_seq_a[i+start_index] = local_nt_b[i]
            ls_seq_b[i+start_index] = local_nt_a[i]
            
        record_a.seq = Seq(str().join(ls_seq_a))
        record_b.seq = Seq(str().join(ls_seq_b))
        crossed_records.append(record_a)
        crossed_records.append(record_b)
    
    return crossed_records + records
    

def mutate(records, tree, generation, mutation_rate, weights, n_children=def_children):
    return_seqs = []
    return_ids = []
    N = 0   # Counter for total number of children added to Tree
    
    for record in records:
        node = find(tree, lambda n: n.id == record.id)
        seq = str(record.seq)
        mutated_variants = []
        mutation_paths = []
        children = 0
        
        # Keep generating recombined children till the target number is satisfied
        while (children < n_children):
            mutations = []
            seq_ls = list(seq)
        
            for i in range(0, mutation_rate):
                mutation_type = np.random.choice(mutation_types,
                                                 p=weights)
                index = np.random.randint(low=0, high=len(seq_ls))
                attempts = 0
                
                while (index in protected_range):
                    attempts += 1
                    index = np.random.randint(low=0, high=len(seq_ls))
                    
                    if  (attempts < 32):
                        break
                
                if (mutation_type == 'point'):
                    local_nt = point_mutate(seq_ls[index])
                    mutation_string = 'M({index}):{A}->{B}'.format(index=index, 
                                                                        A=seq_ls[index], B=local_nt)
                    seq_ls[index] = local_nt
                    
                elif (mutation_type == 'deletion'):
                    mutation_string = 'D({index}):-{A}'.format(index=index, 
                                                              A=seq_ls[index])
                    del seq_ls[index]
                
                elif (mutation_type == 'insertion'):
                    local_nt = insertion()
                    mutation_string = 'I({index}):+{A}'.format(index=index, 
                                                              A=seq_ls[index])
                    seq_ls.insert(index, local_nt)
                                 
            # Add the mutated sequence to the pool
            mutation_paths.append(mutations)
            mut_seq = remove_gaps(str().join(seq_ls))
            mutated_variants.append(mut_seq)
            children += 1
        
        # Update children in Tree
        # record_ids, N = add_children_nodes(seqs=mutated_variants, parent=node, 
        #                                   generation=generation, child_count=N,
        #                                   mutations=mutation_paths)
        
        # Append mutated sequences to return List
        return_seqs += mutated_variants
        return_ids += ['%s-%s' % (record.id, j) for j in range(len(mutated_variants))]
    
    return_records = join_to_records(return_seqs, return_ids)
    return return_records


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
    parser.add_option("-t", "--tree-file", action="store", default="", type="string", dest="treefile", help="Json file to write the resulting tree to.")
    parser.add_option("-g", "--generations", action="store", default=def_generations, type="int", dest="iterations", help="Integer defining number of cycles.")
    parser.add_option("-M", "--model", action="store", default='', type="string", dest="cov_model", help="Path to covariance model used as reference. If none is provided, will be calculated from start population.")
    parser.add_option("-m", "--mutationrate", action="store", default=def_mutation_rate, type="int", dest="mutation_rate", help="Integer defining the base mutation per sequence and generation.")
    parser.add_option("-w", "--weights", action="store", default="", type="string", dest="weights", help="Probabilities of different mutation types, in shape P(point mutation),P(deletion),P(insertion),P(crossover), (Default: 0.3,0.3,0.3,0.1).")
    parser.add_option("-s", "--size", action="store", default=def_selection, type="float", dest="selection", help="Generation size to be selected for the next cycle. Can also be a fraction.")
    parser.add_option("-n", "--children", action="store", default=def_children, type="int", dest="children", help="Initial number of children nodes to generate per generation and element (before scoring).")
    parser.add_option("-c", "--crossover-length", action="store", default=4, type="int", dest="crossover_length", help="Length of crossover segments to swap (Default: 4).")
    options, args = parser.parse_args()
    
    # Prepare base parameters
    aln_file = args[-1]
    iterations = options.iterations
    mutation_rate = options.mutation_rate
    n_children = options.children
    fraction = options.selection
    population_file = 'population.fa'
    cm_file = 'cov_model.cm'
    candidates_file = 'candidates.fa'
    score_file = 'cm_scores.txt'
    cm_alignment = 'cm_alignment.fa'


    # TODO: Deal with U/T in input file
    # TODO: Deal with input duplicates !!
    
    if (len(options.weights) == 0):
        weights = default_weights
    else:
        weights = [float(x) for x in options.weights.split(',')]
    
        if (len(weights) != 3):
            print("Exactly 3 base probabilities/weights have to be supplied.")
            sys.exit()
        elif (sum(weights) != 1):
            print("Weights of mutation types have to sum up to 1!")
            sys.exit()

    # Load the base RNA input file
    records = load_seq_records(aln_file)
    
    # Initialize tree & write initial population file
    root = initialize_tree(options, weights)
    ungapped_seqs = [str(r.seq.ungap('-')) for r in records]
    record_ids, _ = add_children_nodes(seqs=ungapped_seqs, parent=root, generation=0, child_count=0)
    
    for i in range(0, len(records)):
        records[i].id = record_ids[i] 

    write_population_file(records, population_file) 
    
    # Alignment of initial population 
    # population_aln = clustalw_alignment(fasta_file=population_file)
    # aligned_records = load_aligned_records(filename=population_aln)
    
    # Build an initial CM model iff none is provided
    if len(options.cov_model) == 0:
        build_cm_model(fasta_file=population_file, cm_file=cm_file)
    else:
        cm_file = options.cov_model
    
    # Score at zero generation (t=0)
    cm_aligned_records = infernal_align(fasta_file=population_file, cm_file=cm_file,
                                        cm_scores=score_file, cm_alignment=cm_alignment)
    scores = parse_cmscores(cm_scores=score_file)
    score_dict = map_seq_to_score(cm_aligned_records, scores=scores)
    
    # Add metrics to tree
    children = get_generation_nodes(root, generation=0)
    assign_features_to_nodes(children, cm_score_dict=score_dict)
    
    # assign_similarity_scores(nodes=children, aligned_records=aligned_records)
    # Initial aligned population
    # records = aligned_records
    
    
    # Start Generation:
    for generation in range(1, iterations+1):
        print('Generation #%s' % generation)
        synthetic_population = 0
        
        # Mutate current population
        mutated_records = mutate(records=records, tree=root, generation=generation,
                              mutation_rate=mutation_rate, weights=weights,
                              n_children=n_children)
        
        # Write candidates to Fasta:
        write_candidate_seqs(records=mutated_records, filename=candidates_file)
        
        # Align mutated entities to CM model and score them
        cm_aligned_records = infernal_align(fasta_file=candidates_file, cm_file=cm_file,
                                         cm_scores=score_file, cm_alignment=cm_alignment)
        scores = parse_cmscores(cm_scores=score_file)
        score_dict = map_seq_to_score(cm_aligned_records, scores=scores)
        
        # Select best candidate(s)
        selected_records = []
        best_candidates = select_best_scoring(score_dict=score_dict, fraction=fraction)
        best_ids = [x[0] for x in best_candidates]
        
        for record in mutated_records:
            if record.id not in best_ids:
                continue
            selected_records.append(record)
        
        # Crossover and append new children to tree
        crossed_records = crossover(selected_records, crossover_length=options.crossover_length)
        parent_nodes = get_generation_nodes(root, generation=generation-1)
        children = update_tree(crossed_records, parent_nodes)
        
        # Write current updated population
        records = []
        
        for child in children:
            sequence=child.seq
            synthetic_population += 1
            synthetic_record = SeqRecord(
            seq=Seq(sequence),
            id=child.id,
            name="",
            description="",
            )
            records.append(synthetic_record)
            
        write_population_file(records, population_file)
        
        # Update CM scores to all children nodes (for reference)
        cm_aligned_records = infernal_align(fasta_file=population_file, cm_file=cm_file,
                                         cm_scores=score_file, cm_alignment=cm_alignment)
        scores = parse_cmscores(cm_scores=score_file)
        score_dict = map_seq_to_score(cm_aligned_records, scores=scores)
        
        # Assign features to child nodes
        assign_features_to_nodes(children, cm_score_dict=score_dict)
    
        # RESTART CYCLE
        # END OF LOOP
        
    
    # Test: render tree
    print(RenderTree(root))
    
    # Export tree as JSON
    exporter = JsonExporter(indent=2, sort_keys=True)
    
    if len(options.treefile) == 0:
        pass
    else:
        with open(options.treefile, 'w') as f:
            f.write(exporter.export(root))
    
    
if __name__ == '__main__':
    main()
