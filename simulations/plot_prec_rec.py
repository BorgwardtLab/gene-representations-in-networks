#!/usr/bin/env python3
# @author: Anja Gumpinger


import logging
import os
import ipdb

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set()


def main():

    ROOT = '/links/groups/borgwardt/Projects/master_project_julia'
    OUT_DIR = f'{ROOT}/output/simulation'
    OUT_FILE = f'{OUT_DIR}/simuation_prec_rec.pdf'

    # input the data. 
    skato = [0.06756756756756757, 1.0]
    fisher = [0.4,1.0]
    min_pv = [0.2,1.0]

    fig, ax = plt.subplots()

    ax.set_xlabel('recall')
    ax.set_ylabel('precision')
    ax.plot(skato[1], skato[0], 'o', label='SKAT-O')
    ax.plot(fisher[1], fisher[0], 'o', label='Fisher')
    ax.plot(min_pv[1], min_pv[0], 'o', label='minimum p-value')

    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.grid(True)
    ax.legend(ncol=3, facecolor='1.0')

    fig.savefig(OUT_FILE, bbox_inches='tight')

    pass

if __name__ == "__main__":
    main()

