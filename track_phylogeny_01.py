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

# sgfif = load_strain_genome(r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Archives\strain_genome_015a.txt")
# sgten = load_strain_genome(r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Archives\strain_genome_010a.txt")
# overlaps = [x for x in sgfif if x in sgten]
# sgall = sgfif | sgten
# len(sgfif) + len(sgten) - len(sgall)


# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.Align import MultipleSeqAlignment
# def align_multiple_ids(orgid_list, strain_genome_dict):
#     seqrec_list = []
#     for oid in orgid_list:
#         aaff_string = strain_genome_dict[oid]
#         # Unspool each aaff_string for each org_id.
#         nt_seq = unspool_aaff_for_msa(aaff_string)
#         # 'Hard' way of writing an alignment: MSA > SeqRecord > Seq
#         seqrec = SeqRecord(Seq(nt_seq), id=str(oid))
#         # Append to a list of SeqRecords.
#         seqrec_list.append(seqrec)
#     # MSA wrap around this list.
#     ms_align = MultipleSeqAlignment(seqrec_list)
#     return ms_align
# # Write to a PHYLIP file.
# AlignIO.write(my_alignments, "my_example.phy", "phylip")

# # Stops line breaks when showing df in console.
# pandas.set_option('display.expand_frame_repr', False)
# overlaps = [x for x in bfif_df.index if x in blife_df.index]
# # Book from run_010 takes precedence.
# bkall_df = blife_df.combine_first(bfif_df)
# len(bkall_df)
# len(blife_df) + len(bfif_df) - len(overlaps)
# # But the death_step for overlapping organisms is updated from run_015.
# bkall_df.loc[overlaps, "Death_step"] = bfif_df.loc[overlaps, "Death_step"]
# bkall_df["Lifespan"] = blife_df.Death_step - blife_df.Birth_step
# bkall_df.loc[overlaps].head(20)
# bkall_df["Lifespan"] = bkall_df.Death_step - bkall_df.Birth_step
# bkall_df.loc[overlaps].head(20)
# # Too many rows to view entirely in spreadsheet.
# bkall_df.to_csv(os.path.join(path_rundir, "book_of_life_015_010.csv"))


# From a list of org_ids from a file, create an ALN file to align with Unipro UGENE.
def align_many_ids(strain_genome_dict, id_list, prefix, path_aln, pad_digits=0):
    # Follow the format: "CLUSTAL W 2.0 multiple sequence alignment\n\n"
    # Blank line.
    aln_list = ["CLUSTAL W 2.0 multiple sequence alignment\n\n"]
    frag_dict = {}
    flen_max = 0
    for orgid in id_list:
        org_genome = unspool_aaff_for_msa(strain_genome_dict[orgid])
        frag_len = 70
        # Split into fragments of length 70.
        frag_list = [org_genome[s:s + frag_len] for s in range(0, len(org_genome), frag_len)
                     if len(org_genome[s:s + frag_len]) > (frag_len-1)]
        # Add the final fragment of length < 70.
        tail_frag = org_genome[len(frag_list) * frag_len:]
        tail_frag = tail_frag + "-" * (frag_len - len(tail_frag))
        frag_list.append(tail_frag)
        frag_dict[orgid] = frag_list
        
        if len(frag_list) > flen_max:
            flen_max = len(frag_list)
    
    # Pad out the frag_list with dash-only fragments.
    for orgid in id_list:

        if len(frag_dict[orgid]) < flen_max:
            for i in range(flen_max - len(frag_dict[orgid])):
                frag_dict[orgid].append("-" * frag_len)
                
    # Decide how many digits each org_id should have.
    padded_ids_dict = {}
    oid_max = len(str(max(id_list)))
    for orgid in id_list:
        if len(str(orgid)) < oid_max:
            padded_ids_dict[orgid] = "0" * (oid_max - len(str(orgid)) + pad_digits) + str(orgid)
        else:
            padded_ids_dict[orgid] = "0" * pad_digits + str(orgid)

    # Add the fragments by org_id and multiple of 70.
    for i in range(flen_max):
        for orgid in id_list:            
            # Sequence name + 4 spaces + 70 nts + 1 space + multiple of 70.
            # line = prefix + "_" + str(padded_ids_dict[orgid]) + " "*4 + frag_dict[orgid][i] + " " + str((i+1)*70) + "\n"
            # Max length of ID is 10 char, so omit prefix for now.
            line = str(padded_ids_dict[orgid]) + " "*4 + frag_dict[orgid][i] + " " + str((i+1)*70) + "\n"
            aln_list.append(line)
        aln_list.append("\n\n")

    # Write to ALN file.
    with open(path_aln+".aln", "wt") as f:
        for line in aln_list:
            f.write(line)

    # Align with Unipro UGENE via batch file.
    path_batrun_ugene = os.path.join(path_rundir, "run_" + run_num + "_" + run_name + "_ugene.bat")
    copyfile(path_bat_ugenetemp, path_batrun_ugene)
    replace_old_with_new(path_batrun_ugene, {"path_in":path_aln, "path_out":path_aln+"a"})
    # Get the filename itself to run.
    os.system(path_batrun_ugene.split('\\')[-1])

    # Switch from UGENE to ClustalW 2.1.


from Bio import Phylo, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
path_aln_in = r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\Runs\Run_015_big_bang\test_a.aln"
path_phy_in = r'C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\Runs\Run_015_big_bang\testb.phy'
def draw_phylos(path_aln_in, path_phy_in):
    # alignment = AlignIO.read(open(r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\Runs\Run_015_big_bang\test.aln"), "clustal")
    alignment = AlignIO.read(open(path_aln_in), "clustal")
    print("Alignment length %i" % alignment.get_alignment_length())
    for record in alignment:
        print(record.seq + " " + record.id)
    # align = AlignIO.read(r'C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\Runs\Run_015_big_bang\testb.phy','phylip')
    align = AlignIO.read(path_phy_in,'phylip')
    print(align)
    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    distMatrix = calculator.get_distance(align)
    print(distMatrix)
    # Create a DistanceTreeConstructor object
    constructor = DistanceTreeConstructor()# Construct the phlyogenetic tree using UPGMA algorithm
    UPGMATree = constructor.upgma(distMatrix)# Construct the phlyogenetic tree using NJ algorithm
    NJTree = constructor.nj(distMatrix)
    # Draw the phlyogenetic tree.
    Phylo.draw(UPGMATree)
    # Draw the phlyogenetic tree using terminal
    Phylo.draw_ascii(NJTree)

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


# # r.ggenealogy works with df containing 'child' and 'parent.
# # But getParent() only returns 1 value, even though it should return 2 values.    
# # Both the code for getChild() and getParent() are very similar: Selection of column in df.
# # So why can there be multiple children but not multiple parents?
# # I added more rows for the second parent and it worked.
# bol_df = pd.read_csv(r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Archives\book_of_life_015_010.csv")
# onebol_df = bol_df.loc[:,['ID','Sporelayer']]
# onebol_df.rename(columns={'ID':'child', 'Sporelayer':'parent'}, inplace=True)
# twobol_df = bol_df.loc[:,['ID','Quickener']]
# twobol_df.rename(columns={'ID':'child', 'Quickener':'parent'}, inplace=True)
# # See how many 2nd parents there are.
# sexbol_df = bol_df.loc[bol_df.Sex_check != 0]
# threebol_df = pd.concat([onebol_df,twobol_df.loc[bol_df.Sex_check != 0]])
# bbbol_df = pd.read_csv(r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Archives\bol_for_r.csv")
# bbbol_df = pd.read_csv(r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Archives\bol_for_r.csv", index_col="Unnamed: 0")
# Actually there is no limit to the number of parents.

# # Which organisms retained the original genome?
# ooo = {}
# for g in range(215):
#     itmd = bgen_df[bgen_df.Generation == g]
#     ooo[g] = list(itmd[bgen_df.Genome == og_gen].index)
# oiu = pd.DataFrame.from_dict(ooo,'index')