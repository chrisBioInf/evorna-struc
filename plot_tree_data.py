#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 13:54:15 2024

@author: christopher
"""

import sys 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from optparse import OptionParser
from anytree.importer import JsonImporter
from anytree.iterators.levelorderiter import LevelOrderIter
from anytree import RenderTree


# Global plot information
sns.set_theme(style='whitegrid')


def load_tree(filename, options):
    importer = JsonImporter()
    
    with open(filename, 'r') as f:
        data = f.read()    
        tree = importer.import_(data)
    if options.render:
        print(RenderTree(tree))
    return tree


def get_generation_nodes(node, generation):
    iterator = LevelOrderIter(node, filter_=lambda n: n.generation == generation)
    return [node for node in iterator]


def plot_generational_statistics(df, tree, options):
    # root = get_generation_nodes(tree, generation=-1)[0]
    
    title = "N={samples}, generations={generations}, selection={size}, \nmutations={mutation}, children={children}".format(
        samples=df['Samples'].iloc[0],
        generations=df['Generation'].max(),
        size=int(tree.selection),
        mutation=tree.mutations,
        children= int(tree.n_children),
        )
    
    fig, axes = plt.subplots(nrows=2, ncols=2, )
    fig.set_size_inches(15.5, 8.5)

    sns.lineplot(data=df, x='Generation', y='Mean covariance model score', markers=True, hue='Label', ax=axes[0][0])
    axes[0][0].set_ylim((0, 1.1))

    sns.lineplot(data=df, x='Generation', y='Minimum free energy',  hue='Label', ax=axes[0][1])

    sns.lineplot(data=df, x='Generation', y='Mean pairwise identity',  hue='Label', ax=axes[1][0])
    axes[1][0].set_ylim((0, 1.1))

    sns.lineplot(data=df, x='Generation', y='Mean CG content', hue='Label', ax=axes[1][1])
    axes[1][1].set_ylim((0, 1.1))

    if (options.xstep != 0):
        xtick_steps = options.xstep
        axes[0][0].set_xticks([x for x in range(0, max(df['Generation'])+xtick_steps, xtick_steps)])
        axes[0][1].set_xticks([x for x in range(0, max(df['Generation'])+xtick_steps, xtick_steps)])
        axes[1][0].set_xticks([x for x in range(0, max(df['Generation'])+xtick_steps, xtick_steps)])
        axes[1][1].set_xticks([x for x in range(0, max(df['Generation'])+xtick_steps, xtick_steps)])

    fig.suptitle(title)
    
    if (len(options.outfile) > 0):
        plt.savefig(options.outfile, dpi=300, bbox_inches='tight')
    else:
        plt.show()


def main():
    usage = "\nplot_tree_data.py [options] [Json 1] [Json 2] ..."
    parser = OptionParser(usage=usage, version="__version__")
    args = sys.argv
    
    # Define input parameters
    parser.add_option("-g", "--generations", action="store", default=100, type="int", dest="n_generations", help="Number of generations to plot.")
    parser.add_option("-o", "--outfile", action="store", default="", type="string", dest="outfile", help="File to write plot to. Default: Just show.")
    parser.add_option("-l", "--labels", action="store", default="", type="string", dest="labels", help="Comma-separated labels for data sets in trees (Example: tree1,tree2). Default is file names.")
    parser.add_option("-x", "--xsteps", action="store", default=0, type="int", dest="xstep", help="Step size for the x-axis.")
    parser.add_option("-r", "--render-tree", action="store_true", dest="render", help="Should each tree be displayed in terminal?")
    options, args = parser.parse_args()
    
    # Get labels for trees
    if (len(options.labels) == 0):
        file_labels = []
    else:
        file_labels = options.labels.split(',')
        
    label_dict = {}
    for (treefile, label) in zip(args, file_labels):
        label_dict[treefile] = label
    
    # Initialize list of Data Frames
    data_frames = []
    
    for treefile in args:
        root = load_tree(treefile, options)
        
        # Initialize lists for result data frame
        generations = []
        samples = []
        mean_cm_scores = []
        mean_mfes = []
        mean_identity = []
        mean_cg_content = []
        
        for generation in range(0, options.n_generations+1):
            nodes = get_generation_nodes(root, generation)
            
            if (len(nodes) == 0):
                break
            
            # Add averaged parameters for this generation
            generations.append(generation)
            samples.append(len(nodes))
            mean_cm_scores.append(np.mean([n.cm_score for n in nodes]))
            mean_mfes.append(np.mean([n.mfe for n in nodes]))
            mean_identity.append(np.mean([n.mean_hamming_distance for n in nodes]))
            mean_cg_content.append(np.mean([n.CG_content for n in nodes]))
        
        # We will need (at least) 2 generations for the plot to make sense
        if len(generations) < 2:
            print('There is only one generation in the data provided - maybe check tree or --generations parameter?')
            sys.exit()
        
        # Build data frame
        data = {
            'Generation': generations,
            'Samples': samples,
            'Mean covariance model score': mean_cm_scores,
            'Minimum free energy': mean_mfes,
            'Mean pairwise identity': [1-x for x in mean_identity],
            'Mean CG content': mean_cg_content,
            }
        df = pd.DataFrame(data=data)
        df['Label'] = label_dict.get(treefile, treefile)
        data_frames.append(df)
    
    # Connect data frames
    df = pd.concat(data_frames)
    df.index = list(range(0, len(df)))
    print(df)
    
    # Plot the data
    plot_generational_statistics(df, root, options)


main()

