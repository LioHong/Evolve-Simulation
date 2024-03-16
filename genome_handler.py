# -*- coding: utf-8 -*-
"""
Filename: genome_handler.py
Date created: 2024/03/05, Tue, 21:04:00 (UTC+8)
@author: LioHong
Purpose: Packages all functions relating to KFORTH genome handling.
Steps:
"""
import numpy as np
from os import system
from json import dump, loads
from pathlib import Path
from pandas import read_csv

# https://stackoverflow.com/questions/28730961/python-slicing-string-in-three-character-substrings
def pair_split(elm):
    return [elm[s:s + 2] for s in range(0, len(elm), 2) if len(elm[s:s + 2]) > 1]

instr_aaff_dict = eval(Path("instr_aaff_dict.txt").read_text(encoding="utf-8"))
aaff_instr_dict = {v: k for k, v in instr_aaff_dict.items()}


# Small but useful for quick manual exports: Copy list to clipboard for pasting elsewhere.
def addToClipBoard(read_list):
    # Use spaces to separate.
    text = '_'.join(read_list)
    command = 'echo ' + text.strip() + '| clip'
    system(command)


# Setup the dicts as basis for comparison.
def get_organics_from_universe(text_phascii):
    # with open(file_phascii, "rt") as fph:
    #     text_phascii = fph.readlines()

    org_dict = {"SPORE": [], "ORGANISM": []}
    for organic in org_dict:
        for line in text_phascii:
            if organic in line:
                if organic == "SPORE":
                    org_dict[organic].append(line)
                elif organic == "ORGANISM":
                    # Must remove the Energy and Age which can change over time.
                    org_dict[organic].append(line.split(' ')[:-2])
    return org_dict


def store_aaff(aaff_genome):
    aaff_string = ''.join(aaff_genome)
    return aaff_string.replace("__","_")


# Convert string to list of instructions/numbers/rows.
def retrieve_aaff(aaff_string,snum=False):
    # Separate out all the numbers.
    aaff_list = aaff_string.split("_")
    # Remove any empty strings.
    aaff_list = [x for x in aaff_list if x != ""]
    aaff_genome = []
    for elm in aaff_list:
        if not (elm.isdigit() or elm[0] == "-"):
            # Split all the remaining AAFFs into pairs.
            # evd_list = [elm[s:s + 2] for s in range(0, len(elm), 2) if len(elm[s:s + 2]) > 1]
            evd_list = pair_split(elm)
            # Flatten the list.
            for evd in evd_list:
                aaff_genome.append(evd)
        # Filter out all the numbers
        else:
            if not snum:
                aaff_genome.append(elm)
            else:
                aaff_genome.append(str(elm))
    return aaff_genome


# Assume input genome is a list.
def shrink_kforth_to_aaff(kforth_genome):
    # Convert instructions.
    aaff_genome = [instr_aaff_dict[i] if i in instr_aaff_dict else i for i in kforth_genome]
    # Segregate positive and negative numbers.
    aaff_genome = ["_" + x + "_" if x.isdigit() else x for x in aaff_genome]
    aaff_genome = ["_" + x + "_" if x[1:].isdigit() else x for x in aaff_genome]
    # Convert rows: Base-26.
    aaff_genome = [chr(70 + int(x[3:])//26) + chr(64+int(x[3:])%26) if x[:3] == 'row' else x for x in aaff_genome]
    return aaff_genome


# To eyeball the genome.
def translate_aaff_to_kforth(aaff_string):
    aaff_genome = retrieve_aaff(aaff_string)
    # Instructions.
    kforth_genome = [aaff_instr_dict[a] if a in aaff_instr_dict.keys() else a for a in aaff_genome]
    # Rows, FA ~ row1.
    kforth_genome = ['row'+str((ord(a[0])-ord('F'))*26+ord(a[1])-ord('A')+1)
                     if a.isupper() and len(a)==2 else a for a in kforth_genome]
    return ' '.join(kforth_genome)


# Split into rows for readability.
def split_kforth_to_read(aaff_string):
    kforth_string = translate_aaff_to_kforth(aaff_string)
    kforth_rows_list = kforth_string.split("row")
    kforth_rows_list = [kforth_rows_list[0]] + ["row"+row for row in kforth_rows_list[1:]]
    return kforth_rows_list


def load_strain_genome(path_sgen):
    # Open the text archive.
    with open(path_sgen, "rt") as f:
        sfgen = f.readlines()
    keys = [int(x.split(" ")[0]) for x in sfgen]
    values = [x.split(":")[1][:-2] for x in sfgen]
    sfgen_dict = dict(zip(keys, values))
    return sfgen_dict


# Combine book_of_life and strain_genome.
def stitch_sgen(sbol_df, path_sgen):
    bgen_df = sbol_df[:]
    sfgen_dict = load_strain_genome(path_sgen)
    bgen_df['Genome'] = bgen_df.index.map(sfgen_dict)
    return bgen_df


def isindict(ipdf, strx):
    try:
        return ipdf[strx]
    except:
        return strx


# Use keys to point to genomes.
# Work with bgen_df first.
def cache_genomes(bgen_df,threshold=5):
    # Store progenitor genome: orgid=0 or smallest.
    # Approach #1: This is ordered chronologically.
    cache_gnm_dict = {1:bgen_df.loc[bgen_df.index[1],'Genome']}
    # Set threshold.
    # Check len of dict and use it as key for next genome.
    # Approach #2: This is ordered by frequency.
    gnm_counts = bgen_df.Genome.value_counts()
    gnm_counts = gnm_counts[gnm_counts>threshold]
    cgd = {k:v for k, v in zip(range(len(gnm_counts)), gnm_counts.index)}
    dgc = {v:k for k, v in cgd.items()}
    bgen_df['cgen'] = bgen_df.Genome.apply(lambda x: isindict(dgc, x))
    return bgen_df, cgd


# Reverse of cache_genomes().
def expand_genomes(cgen_p, cgd_p):
    cgen_df = read_csv(cgen_p)
    clines = cgd_p.read_text(encoding="utf-8")
    cgd = loads(clines)
    # Insert 'Genome' before 'cgd' like in bgen_df.
    cgen_df.insert(8, 'Genome', 0)
    cgen_df['Genome'] = cgen_df.cgen.apply(lambda x: isindict(cgd, x))
    cgen_df.set_index('ID',inplace=True)
    return cgen_df, cgd


# Leave only combined, cached bgen_df - cgen_df.
def compress_book(bgen_df, path_bgen, path_sgen, path_cgen, path_cgd, threshold=5):
    bgen_df = stitch_sgen(bgen_df, path_sgen)
    bgen_df.to_csv(path_bgen)
    cgen_df, cgd = cache_genomes(bgen_df, threshold)
    cgen_df.drop(columns=['Genome']).to_csv(path_cgen)
    with open(path_cgd, "w") as outfile:
        dump(cgd, outfile)
    return cgen_df, cgd


# Another problem: How to stitch book_of_life and strain_genome from consecutive runs?
def collate_books(run_nums_list,grp_num="002"):
    # Assume list of ints.
    run_nums_list.sort()
    r_nstr_list = [f"{x:03}" for x in run_nums_list]
    last = f"{run_nums_list[-1]:03}"
    # Merge with the later run as priority.
    coll_b = []
    coll_s = []
    for r in reversed(r_nstr_list):
        print("Progress update at " + datetime.now().strftime("%H:%M:%S"))
        rdpath = Path(".") / "Runs" / ("Grp_" + grp_num) / ("Run_" + r)
        bkpath = rdpath / ("book_of_life_" + r + ".txt")
        sgpath = rdpath / ("strain_genome_" + r + ".txt")
        bk = bkpath.read_text(encoding="utf-8").splitlines()
        sg = sgpath.read_text(encoding="utf-8").splitlines()
        # Refer from book_of_life because smaller filesize is easier to manage.
        orgids = [b.split(' ')[0] for b in bk]
        if (coll_b):
            og = [b.split(' ')[0] for b in coll_b]
            # Add orgids only if not in later run.
            # add_b = [bk[orgids.index(bb)] for bb in orgids if bb not in og]
            # lost_ones = [orgids.index(bb) for bb in orgids if bb not in og]
            lost_ones = [bb for bb in orgids if bb not in og]
            lost_ones = [orgids.index(lo) for lo in lost_ones]
            add_b = [bk[i] for i in lost_ones]
            add_s = [sg[i] for i in lost_ones]
            # Combine lists.
            coll_b = coll_b + add_b
            coll_s = coll_s + add_s
            coll_b.sort()
            coll_s.sort()
        else:
            coll_b = bk[:]
            coll_s = sg[:]
        (rdpath / ("book_of_life_" + last + "_combo.txt")).write_text("\n".join(coll_b), encoding="utf-8")
        (rdpath / ("strain_genome_" + last + "_combo.txt")).write_text("\n".join(coll_s), encoding="utf-8")


# Measure genome length
# pandas filter genome length

# Branch from execution of next batch
# Input blacklist
# Find organism and all its cells
# Maybe have to find the next orgid, then find the line before that
# Last orgid: Just find the last line.
# Remove current EVOLVE universe. (Is this redundant?)
# Convert PHASCII to new EVOLVE.