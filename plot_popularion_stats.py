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

from Bio import SeqIO
from seqgen import rna_mfe, mean_pairwise_identity, CG_content


# Global plot information
sns.set_theme(style='whitegrid')
xtick_steps = 10

param_g = 100
param_c = 2
param_s = 100
param_n = 1000
param_m = 1

title = "N={samples}, generations={generations}, selection={size}, \nmutations={mutation}, children={children}".format(
    samples=param_n,
    generations=param_g,
    size=param_s,
    mutation=param_m,
    children=param_c,
    )


dfs = [pd.read_csv(f, sep='\t') for f in sys.argv[1:]]
df = pd.concat(dfs)

reference_file = 'Models/RF00005.stockholm.txt'
records = [r for r in SeqIO.parse(open(reference_file, 'r'), format='stockholm')]
# print(records)

mean_ref_identity = 0.6142    # mean_pairwise_identity(records)
mean_ref_mfe = -22.3752       # np.mean([rna_mfe(str(r.seq).replace('-', '')) for r in records])
mean_ref_cg = np.mean([CG_content(str(r.seq.ungap('-'))) for r in records])
 
fig, axes = plt.subplots(nrows=2, ncols=2, )
fig.set_size_inches(15.5, 8.5)

sns.lineplot(data=df, x='Generation', y='Mean covariance model score', markers=True, hue='Label', ax=axes[0][0])
# sns.scatterplot(data=df, x='Generation', y='Mean covariance model score', hue='Label', ax=axes[0][0])
axes[0][0].set_ylim((0, 1.1))

sns.lineplot(data=df, x='Generation', y='Minimum free energy',  hue='Label', ax=axes[0][1])
# sns.scatterplot(data=df, x='Generation', y='Minimum free energy', hue='Label', ax=axes[0][1])
axes[0][1].plot([0, max(df['Generation'])], [mean_ref_mfe, mean_ref_mfe], 'r--')
# axes[0][0].set_ylim((min(df['Minimum free energy']) +5, max(df['Minimum free energy'])-5))

sns.lineplot(data=df, x='Generation', y='Mean pairwise identity',  hue='Label', ax=axes[1][0])
# sns.scatterplot(data=df, x='Generation', y='Mean pairwise identity', hue='Label', ax=axes[1][0])
axes[1][0].plot([0, max(df['Generation'])], [mean_ref_identity, mean_ref_identity], 'r--')
axes[1][0].set_ylim((0, 1.1))

sns.lineplot(data=df, x='Generation', y='Mean CG content', hue='Label', ax=axes[1][1])
# sns.scatterplot(data=df, x='Generation', y='Mean CG content', ax=axes[1][1])
axes[1][1].plot([0, max(df['Generation'])], [mean_ref_cg, mean_ref_cg], 'r--')
axes[1][1].set_ylim((0, 1.1))

axes[0][0].set_xticks([x for x in range(0, max(df['Generation'])+xtick_steps, xtick_steps)])
axes[0][1].set_xticks([x for x in range(0, max(df['Generation'])+xtick_steps, xtick_steps)])
axes[1][0].set_xticks([x for x in range(0, max(df['Generation'])+xtick_steps, xtick_steps)])
axes[1][1].set_xticks([x for x in range(0, max(df['Generation'])+xtick_steps, xtick_steps)])

fig.suptitle(title)

plt.savefig("Figures/training_results.pdf", dpi=300, bbox_inches='tight')
plt.savefig("Figures/training_results.svg", dpi=300, bbox_inches='tight')
plt.show()
