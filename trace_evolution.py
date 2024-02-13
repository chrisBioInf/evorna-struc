#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:48:17 2024

@author: christopher
"""


import sys 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from optparse import OptionParser
from anytree.importer import JsonImporter
from anytree.search import find
from anytree.walker import Walker
from anytree.iterators.levelorderiter import LevelOrderIter
from anytree import RenderTree, AnyNode


# Global plot information
sns.set_theme(style='whitegrid')
sns.color_palette("Set2")


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


def get_minimal_ancestry(ancestor_dict):
    print(ancestor_dict)
    in_order = sorted(list(ancestor_dict.keys()))
    nodes = [AnyNode(id='root', children=[])]
    
    for id_ in in_order:
        if (len(ancestor_dict.get(id_)) == 0):
            node = AnyNode(id=id_, children=[])
            nodes.append(node)
            continue
       
        parent_id = max(ancestor_dict.get(id_))
        parent_node = next(n for n in nodes if n.id == parent_id)
        node = AnyNode(id=id_, children=[], parent=parent_node)
        nodes.append(node)
        
    return nodes[0]


def find_last_common_ancestor(ancestorA, ancestorB):
    while (ancestorA.generation > ancestorB.generation):
        ancestorA = ancestorA.parent
        
    while (ancestorB.generation > ancestorA.generation):
        ancestorB = ancestorB.parent
    
    while (ancestorA != ancestorB):
        ancestorA = ancestorA.parent
        ancestorB = ancestorB.parent
        
    return ancestorA


def find_ancestral_splits(leaves):
    ancestor_dict = {}

    while len(leaves) > 1:
        ancestral_splits = []
        
        for i in range(0, len(leaves)-1):
            node1 = leaves[i]
            for j in range(i+1, len(leaves)):
                node2 = leaves[j]
                lca = find_last_common_ancestor(node1, node2)
                
                for node in (node1, node2):
                    if (node.id not in ancestor_dict):
                        ancestor_dict[node.id] = set()
                    if (lca.id != node.id):
                        ancestor_dict[node.id].add(lca.id)
                
                if (lca not in ancestral_splits):
                    ancestral_splits.append(lca)
                    
        leaves = ancestral_splits
    
    min_ancestry = get_minimal_ancestry(ancestor_dict)
    return min_ancestry
    

def print_lineage(path):
    print(path[0].seq)
    if len(path) == 1:
        return 
    
    for node in path[1:]:
        print('|')
        print(node.seq)


def get_evo_features(path):
    generations = []
    cm_scores = []
    mfes = []
    cg_contents = []
    
    for node in path:
        generations.append(node.generation)
        cm_scores.append(node.cm_score)
        mfes.append(node.mfe)
        cg_contents.append(node.CG_content)
        
    data = {
        'Generation': generations,
        'Covariance model score': cm_scores,
        'Minimum free energy': mfes,
        'CG content': cg_contents,
        }
    df = pd.DataFrame(data=data)
    df['Sequence'] = path[-1].id
    return df

        
def plot_evolution(df, tree, outfile, xtick_steps=0):
    df = df.sort_values(by="Sequence", ascending=True)
    fig, axes = plt.subplots(nrows=2, ncols=2, )
    fig.set_size_inches(15.5, 8.5)
    
    title = "generations={generations}, selection={size}, \nmutations={mutation}, children={children}".format(
        generations=df['Generation'].max(),
        size=int(tree.selection),
        mutation=tree.mutations,
        children= int(tree.n_children),
        )

    sns.lineplot(data=df, x='Generation', y='Covariance model score', markers=True, hue='Sequence', ax=axes[0][0])
    axes[0][0].set_ylim((0, 1.1))

    sns.lineplot(data=df, x='Generation', y='Minimum free energy',  hue='Sequence', ax=axes[0][1])

    sns.lineplot(data=df, x='Generation', y='CG content', hue='Sequence', ax=axes[1][0])
    axes[1][0].set_ylim((0, 1.1))
    
    for ax in (axes[0][0], axes[1][0], axes[1][1]):
        ax.legend([],[], frameon=False)

    sns.move_legend(axes[0][1], "upper left", bbox_to_anchor=(1, 1.02), title='Species')
    fig.suptitle(title)

    if (xtick_steps != 0):
        axes[0][0].set_xticks([x for x in range(0, max(df['Generation'])+xtick_steps, xtick_steps)])
        axes[0][1].set_xticks([x for x in range(0, max(df['Generation'])+xtick_steps, xtick_steps)])
        axes[1][0].set_xticks([x for x in range(0, max(df['Generation'])+xtick_steps, xtick_steps)])
        axes[1][1].set_xticks([x for x in range(0, max(df['Generation'])+xtick_steps, xtick_steps)])
    
    if (len(outfile) > 0):
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
    else:
        plt.show()


def main():
    usage = "\ntrace_evolution.py --tree [Json] [options] [Id1] [Id2] ..."
    parser = OptionParser(usage=usage, version="__version__")
    args = sys.argv
    
    # Define input parameters
    parser.add_option("-t", "--tree", action="store", default="", type="string", dest="treefile", help="The tree file to analyse (Required).")
    parser.add_option("-o", "--outfile", action="store", default="", type="string", dest="outfile", help="File to write plot to. Default: Just show.")
    parser.add_option("-x", "--xsteps", action="store", default=0, type="int", dest="xstep", help="Step size for the x-axis.")
    parser.add_option("-r", "--render-tree", action="store_true", dest="render", help="Set this flag to render tree in the terminal.")
    parser.add_option("-l", "--print-lineage", action="store_true", dest="print", help="Set this flag to print complete lineage for each id.")
    parser.add_option("-p", "--phylogeny", action="store_true", dest="phylogeny", help="Set this flag to display simplified shared phylogeny of given ids.")
    options, args = parser.parse_args()

    if (len(options.treefile) == 0):
        print("We will need a tree file in JSON format at the very least! :)")
        print("I.e.: --tree [JSON]")
        sys.exit()

    root = load_tree(options.treefile, options)
    leaves = []
    dfs = []
    
    for nodearg in args:
        try:
            node = find(root, lambda n: n.id == nodearg)
            leaves.append(node)
        except Exception:
            print("No fitting Node found. Wrong id?")
            continue
        
        walker = Walker()
        path = walker.walk(root, node)[2:][0]
        dfs.append(get_evo_features(path))
        
        if options.print:
            print_lineage(path)
            
    if options.phylogeny:
        minimal_ancestry = find_ancestral_splits(leaves)
        print(RenderTree(minimal_ancestry))
    
    return
    df = pd.concat(dfs)
    df.index = list(range(0, len(df)))
    plot_evolution(df, tree=root, outfile=options.outfile, 
                   xtick_steps=options.xstep)
    

main()
