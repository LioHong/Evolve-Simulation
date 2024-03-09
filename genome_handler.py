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
from json import dump
from pathlib import Path

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
def stitch_sgen(bgen_df, path_sgen):
    sfgen_dict = load_strain_genome(path_sgen)
    bgen_df['Genome'] = bgen_df.index.map(sfgen_dict)
    return bgen_df


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


    def isindict(strx):
        try:
            return dgc[strx]
        except:
            return strx
    bgen_df['cgen'] = bgen_df.Genome.apply(lambda x: isindict(x))

    return bgen_df, cgd


# Leave only combined, cached bgen_df - cgen_df.
def compress_book(bgen_df, path_bgen, path_sgen, path_cgen, path_cgd, threshold=5):
    bgen_df = stitch_sgen(bgen_df, path_sgen)
    bgen_df.to_csv(path_bgen)
    cgen_df, cgd = cache_genomes(bgen_df, threshold)
    cgen_df.drop(columns=['Genome']).to_csv(path_cgen)
    with open(path_cgd, "w") as outfile:
        dump(cgd, outfile)
    return cgen_df, cgd


# Get genome using key.
