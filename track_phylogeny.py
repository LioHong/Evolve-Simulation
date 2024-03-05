# -*- coding: utf-8 -*-
"""
Filename: track_phylogeny_01.py
Date created: 2024/03/05, Tue, 21:23:41 (UTC+8)
@author: LioHong
Purpose: Package all attempts at tracking phylogenies.
Steps:

"""
import os
from shutil import copyfile
from math import log10
from datetime import datetime
import pandas as pd
import numpy as np
# This is a borrowed algorithm.
import GlobalAlignment
import genome_handler

# To adjust the dataframe appearance
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option("display.width", 200)
pd.set_option('display.expand_frame_repr', False)


def find_parents(orgid, bol_df_in):
    sporelayer = bol_df_in.loc[orgid, 'Sporelayer']
    quickener = bol_df_in.loc[orgid, 'Quickener']
    if sporelayer == quickener:
        parentage = [sporelayer]
    else:
        parentage = [sporelayer, quickener]
    return parentage


def find_children(orgid, bol_df_in, sex=0):
    s_children = bol_df_in[bol_df_in.Sporelayer == orgid].index
    q_children = bol_df_in[bol_df_in.Quickener == orgid].index
    combo_children = list(set(list(s_children) + list(q_children)))
    combo_children.sort()
    return combo_children

    # s_children = bol_df_in[bol_df_in.Sporelayer == orgid].index
    # if not sex:
    #     combo_children = list(set(list(s_children)))
    # # Quickener short-circuits the generations.
    # else:
    #     q_children = bol_df_in[bol_df_in.Quickener == orgid].index
    #     combo_children = list(set(list(s_children) + list(q_children)))
    # return combo_children


def find_ancestors(orgid, bol_df_in, dist=3):
    gen = bol_df_in.loc[orgid, 'Generation']
    gen_min = gen - dist
    if gen_min < 0:
        gen_min = 0
    # # Slicing breaks on the edge case of a parent from a generation out of bounds.
    # # But also short-circuits the generations in-between.
    # lineage_df = bol_df_in
    lineage_df = bol_df_in[(bol_df_in.Generation >= gen_min) & (bol_df_in.Generation <= gen)]
    # Normally equal to dist but depends on bounds.
    time_jump = gen - gen_min
    return find_lineal_kin(orgid, lineage_df, time_jump, 'up')


def find_descendants(orgid, bol_df_in, dist=3):
    gen = bol_df_in.loc[orgid, 'Generation']
    gen_max = gen + dist
    if gen_max > len(bol_df_in):
        gen_max = len(bol_df_in)
    # # See find_ancestors().
    # lineage_df = bol_df_in
    lineage_df = bol_df_in[(bol_df_in.Generation >= gen+1) & (bol_df_in.Generation <= gen_max)]
    # Normally equal to dist but depends on bounds.
    time_jump = gen_max - gen
    return find_lineal_kin(orgid, lineage_df, time_jump, 'down')


# up_down means ancestor or descendant.
def find_lineal_kin(orgid, lineage_df, time_jump, up_down):
    # Take a snapshot first.
    root = [orgid]
    lineal_kin = []
    outofbounds_kin = []

    for i in range(time_jump):
        rootlets = []
        for r in root:
            # Skip problematic roots.
            try:
                if up_down == 'up':
                    immed = find_parents(r, lineage_df)
                elif up_down =='down':
                    immed = find_children(r, lineage_df)
                rootlets.append(immed)
            # Parent comes from a generation out of bounds.
            except KeyError:
                outofbounds_kin.append(r)
        rootlets = [i for immed in rootlets for i in immed]
        rootlets = list(set(rootlets))
        rootlets.sort()
        lineal_kin.append(rootlets)
        root = rootlets
    lineal_kin = [rl for rootlets in lineal_kin for rl in rootlets]
    lineal_kin = list(set(lineal_kin))
    lineal_kin.sort()

    # # Remove false positives, but not sure why it shows up to begin with.
    # outofbounds_kin = [o for o in outofbounds_kin if o not in lineage_df.index]
    # outofbounds_kin = [o for o in outofbounds_kin if o not in lineal_kin]
    print(outofbounds_kin)
    return lineal_kin


def naive_align(seq1, seq2):
    # Use seq1 as template.
    seq1 = retrieve_aaff(seq1)
    seq2 = retrieve_aaff(seq2)
    seq_diff = []
    # Shortcut: Check if lengths are equal.
    if len(seq1) == len(seq2):
        # Shortcut: Check if last few words are equal ~ 3 words.
        if seq1[-3:] == seq2[-3:]:
            # seq_diff = ['|' if w1 in seq2 else 'X' for w1 in seq1]
            for i in range(len(seq1)):
                if seq1[i] == seq2[i]:
                    seq_diff.append('|')
                else:
                    seq_diff.append('X')
    alignment = pd.DataFrame(data=[seq1,seq2,seq_diff])
    return alignment



# Convert row numbers from format of 'BX_#_' to 'F-'
def compress_row_aaff(aaff_in):
    # Split by 'BX'.
    af_list = geha.aaff_in.split('BX_')
    # Find first '_' in non-main substrings.
    # Find the digit.
    # Replace the digit using ord() and offset.

    # Alt pattern: Numerical.
    # fa_list = ['F' + str(chr(x[x.index('_')]+64)) + x[(x.index('_')+1):] for x in af_list[1:]]
    fa_list = ['F' + str(chr(af_list.index(x)+64))+ x[af_list.index(x)//10+2:] for x in af_list[1:]]
    fa_list.insert(0, af_list[0])
    return ''.join(fa_list)


# Compare genome of parent/child or child/parent, then highlight gen of change.
# Very slow: 36 hr for 500K orgids.
def trace_asexual_lines(bol_df_in):
    lines = {}
    pretracked_orgids = []
    for i in range(len(bol_df_in)-1,-1,-1):
        if i in pretracked_orgids:
            continue
        sex_num = bol_df_in.loc[i, 'Sex_check']
        # If orgid arose sexually, add its line and then move on.
        if bol_df_in.loc[i, 'Sex_check']:
            lines[i] = 0
        else:
            p = i
            # Track within line of descent
            while not sex_num:
                if not p:
                    break
                p = find_parents(p, bol_df_in)[0]
                sex_num += bol_df_in.loc[p, 'Sex_check']
                # Remove orgid from tracking.
                pretracked_orgids.append(p)
                print(p)
                print("Progress update at " + datetime.now().strftime("%H:%M:%S"))
            # Track last asexual ancestor.
            lines[i] = bol_df_in.loc[i, 'Generation'] - bol_df_in.loc[p, 'Generation']

        print(i)
        print("Progress update at " + datetime.now().strftime("%H:%M:%S"))
    return lines


# Receives index of lines.
def sift_lines_by_length(lines_df, bol_in_df, min_len=70):
    sift_df = lines_df[lines_df['0'] > 70]
    changers = {}
    for line in sift_df.index:
        lilen = lines_df.loc[line,'0']
        anc = sbol_df.loc[find_ancestors(line, sbol_df, dist=lilen)]
        anc_counts = len(anc.glen.value_counts())
        if anc_counts < 3:
            continue
        else:
            changers[line] = anc_counts
    return changers


# phylotrackpy takes in Organism objects.
# https://stackoverflow.com/questions/53192602/convert-a-pandas-dataframe-into-a-list-of-objects
class Evorg(object):
    def __init__(self, Orgid, Sporelayer, Quickener, Generation, Birth_step, Death_step, Lifespan, Sex_check, Gnm_Len, Genome):
        self.Orgid = Orgid
        self.Sporelayer = Sporelayer
        self.Quickener = Quickener
        self.Generation = Generation
        self.Birth_step = Birth_step
        self.Death_step = Death_step
        self.Lifespan = Lifespan
        self.Sex_check = Sex_check
        self.Gnm_Len = Gnm_Len
        self.Genome = Genome
        self.taxon = Orgid
        # self.Genome = retrieve_aaff(Genome)
    def __repr__(self):
        return "Evorg object " + self.Genome