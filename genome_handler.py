# -*- coding: utf-8 -*-
"""
Filename: genome_handler.py
Date created: 2024/03/05, Tue, 21:04:00 (UTC+8)
@author: LioHong
Purpose: Packages all functions relating to KFORTH genome handling.
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

# To adjust the dataframe appearance
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option("display.width", 200)
pd.set_option('display.expand_frame_repr', False)


# Small but useful for quick manual exports: Copy list to clipboard for pasting elsewhere.
def addToClipBoard(read_list):
    # Use spaces to separate.
    text = '_'.join(read_list)
    command = 'echo ' + text.strip() + '| clip'
    os.system(command)


# Setup the dicts as basis for comparison.
def get_organics_from_universe(file_phascii):
    with open(file_phascii, "rt") as fph:
        text_phascii = fph.readlines()

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



# These conversion functions are packaged primarily for readability of wrangle_genome() and also conversion of archives.
# Testing an alternative storage method.

# Formats in use: KFORTH kewyord, base4 evodon, nucleotide evodon, AAFF evodon
# Dicts currently added: KFORTH to base4, base4 to nt (simple), base4 to AAFF
# From KFORTH to base4, convert the raw genome into processed genome: keywords separated by spaces.
# Then convert keywords into base4 evodons and then concat into a single string.
# Base4 and nt are equivalent. Just use string.replace() or list comprehension and string.join(list).
# Intuitively, I think list comprehension is slower because of the need to iterate over all of the elements. Unless string.replace() is already a list comprehension.
# From length-9 evodons to formats of other lengths like KFORTH or AAFF, need to split the genome into substrings of length 9.

# Also consider intervals of variable length since all the simulations are deterministic.


def convert_b4_to_decimal(evd):
    # 1 & 0 are positive.
    if evd[0] == '0':
        dec_num = 0
    if evd[0] == '1':
        dec_num = 2 ** 16
    # 2 & 3 are negative.
    elif evd[0] == '2':
        dec_num = -2 ** 17
    elif evd[0] == '3':
        dec_num = -2 ** 16

    for i in range(1, len(evd)):
        real_num = int(evd[i]) * 4 ** (len(evd) - i - 1)
        dec_num += real_num
    return dec_num


# Add another function around this for spacers.
# Spacers on either side would be safe, but would increase the memory usage by numbers.
# Change to "_" for easier selection of aaff_strings.
# Is there a way to detect the first occurrence of a number after AAFF evodons?
# Maybe slice the genome into two types of sublists: Evodons only and numbers only.
# Sort them into a dict with their relative order as the key, and the sublist as the value.
# KIV as another possible data compression technique.
def convert_b4_to_intstr(evd):
    dec_num = convert_b4_to_decimal(evd)
    # return str(dec_num) + "_"
    return "_" + str(dec_num) + "_"


# Base-4 format used for compatibility with bits and ATGC.
# Consider whether to keep genomes as type list, then ''.join(list) only for the MSA and data storage?
# Most likely usage: Pick genomes from data storage, then convert them into nucleotides for MSA.
# Base-4 format used for compatibility with bits and ATGC.
def convert_kforth_to_base4(kforth_genome):
    def rownum_to_b4(rn):
        elm = np.binary_repr(2**17 - 1 - 48 - int(rn), width=18)
        ppairs_list = pair_split(elm)
        return ''.join([str(int(p,2)) for p in ppairs_list])
        # return ''.join([str(int(elm[s:s + 2],2)) for s in range(0, len(elm), 2) if len(elm[s:s + 2]) > 1])

    # When converting genomes from instructions and ints to base-pairs, convert instructions first.
    b4_genome = [instr_b4_dict[i] if i in instr_b4_dict else i for i in kforth_genome]
    # Then convert rows with numbering.
    # b4_genome = ['F' + str(chr(int(i[3:])+64)) if i[:3] == 'row' else i for i in b4_genome]
    b4_genome = [rownum_to_b4(i[3:]) if i[:3] == 'row' else i for i in b4_genome]
    # Then convert ints.
    b4_genome = [num_b4_dict[keynum] if keynum in num_b4_dict else keynum for keynum in b4_genome]
    return b4_genome


def convert_base4_to_kforth(b4_genome):
    fa_num = 2 ** 17 - 1 - 48
    def test(evnum):
        try:
            if int(evnum,4) > (fa_num-52):
                return int(evnum,4) <= fa_num
            else:
                return False
        except:
            return False
    # Invert dicts
    b4_instr_dict = {v: k for k, v in instr_b4_dict.items()}
    br_num_dict = {v: k for k, v in num_b4_dict.items()}
    kforth_genome = [b4_instr_dict[evd_b4] if evd_b4 in b4_instr_dict else evd_b4 for evd_b4 in b4_genome]
    # Convert row numbering.
    kforth_genome = ['row'+str(fa_num-int(evnum,4)) if test(evnum) else evnum for evnum in kforth_genome]
    kforth_genome = [br_num_dict[evnum] if evnum in br_num_dict else evnum for evnum in kforth_genome]

    return kforth_genome


# Full conversion to ATGC.
def convert_base4_to_nucleotide(b4_genome):
    nt_genome = []
    for elm in b4_genome:
        # nt_genome = [b4_nt_dict[digit] for digit in b4_genome]
        for digit in b4_nt_dict:
            elm = elm.replace(digit, b4_nt_dict[digit])
        nt_genome.append(elm)
    return nt_genome


def join_genome_for_msa(nt_genome):
    return ''.join(nt_genome)


# Not sure why this would ever be needed but just include it for the sake of completeness.
def convert_nucleotide_to_base4(nt_genome):
    b4_genome = nt_genome
    nt_b4_dict = {v: k for k, v in b4_nt_dict.items()}
    for nt in nt_b4_dict:
        b4_genome = b4_genome.replace(nt, nt_b4_dict[nt])
    return b4_genome


# Convert keywords to AA-FF for storage.
def convert_base4_to_aaff(b4_genome):
    fa_num = 2 ** 17 - 1 - 48
    def test(evd):
        try:
            if int(evd,4) > (fa_num-52):
                return int(evd,4) <= fa_num
            else:
                return False
        except:
            return False
    aaff_genome = [b4_aaff_dict[evd] if evd in b4_aaff_dict else evd for evd in b4_genome]
    # Convert row numbering.
    aaff_genome = ['F'+chr(fa_num-int(evd,4)+64) if test(evd) else evd for evd in aaff_genome]
    # Converts to integer string.
    aaff_genome = [str(convert_b4_to_intstr(evd)) if len(evd) > 2 else evd for evd in aaff_genome]
    # aaff_genome = ["_"+str(int(evd,4))+"_" if len(evd) > 2 else evd for evd in aaff_genome]
    return aaff_genome


def store_aaff(aaff_genome):
    aaff_string = ''.join(aaff_genome)
    return aaff_string.replace("__","_")


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


# Data de-compression.
def convert_aaff_to_base4(aaff_genome):
    def fa_to_b4(fa_a):
        elm = np.binary_repr(2**17 - 1 - 48 - ord(fa_a)+64, width=18)
        ppairs_list = pair_split(elm)
        return ''.join([str(int(p,2)) for p in ppairs_list])

    b4_genome = [aaff_b4_dict[evd] if evd in aaff_b4_dict else evd for evd in aaff_genome]
    b4_genome = [fa_to_b4(evd[1]) if evd[0] == 'F' else evd for evd in b4_genome]
    b4_genome = [num_b4_dict[keynum] if keynum in num_b4_dict else keynum for keynum in b4_genome]
    return b4_genome


def convert_aaff_to_xascii(aaff_genome):
    # Converts the double-letter format.
    xa_genome = [aaff_xascii_dict[evd] if evd in aaff_xascii_dict else evd for evd in aaff_genome]
    # Adds spacers to numbers.
    # This accounts for single digits.
    xa_genome = ["_"+x+"_" if x.isdigit() else x for x in xa_genome]
    # This accounts for negative numbers. Single digits not included due to trailing underscore.
    xa_genome = ["_"+x+"_" if x[1:].isdigit() else x for x in xa_genome]
    return xa_genome


# Single function to convert chosen AAFF strings into nucleotide sequences for MSA.
def unspool_aaff_for_msa(aaff_string):
    aaff_genome = retrieve_aaff(aaff_string)
    b4_genome = convert_aaff_to_base4(aaff_genome)
    nt_genome = convert_base4_to_nucleotide(b4_genome)
    nt_seq = join_genome_for_msa(nt_genome)
    return nt_seq


# This function is used for variable inputs.
def check_input_files(path_run_in):
    # Read the filenames of the template files.
    files_list = os.listdir(path_run_in)
    files_list = [x.split('.')[0] for x in files_list]
    # Check that the names of the EVOLVE and PHASCII both match.
    if len(set(files_list)) > 1:
        print("Mismatch in EVOLVE and PHASCII detected.")
    # Then extract the timestep.
    else:
        return int(files_list[0].split('_')[-1])


# To eyeball the genome.
def translate_aaff_to_kforth(aaff_string):
    aaff_genome = retrieve_aaff(aaff_string)
    b4_genome = convert_aaff_to_base4(aaff_genome)
    kforth_genome = convert_base4_to_kforth(b4_genome)
    kforth_string = ' '.join(kforth_genome)
    # kforth_string = kforth_string.replace("row ", "row")
    return kforth_string


# Split into rows for readability.
def split_kforth_to_read(aaff_string):
    kforth_string = translate_aaff_to_kforth(aaff_string)
    kforth_rows_list = kforth_string.split("row")
    kforth_rows_list = [kforth_rows_list[0]] + ["row"+row for row in kforth_rows_list[1:]]
    return kforth_rows_list


# Fix negative numbers in strain_genome_015 as they were far too large.
def fix_negnum_in_aaff(path_in):
    with open(path_in, "rt") as sgen:
        wrong = sgen.readlines()
    rabbit = [x.split("_") for x in wrong]
    turtles = []
    for bun in rabbit:
        ninja = []
        for b in bun:
            if b.isdigit() and len(b) == 5 and b[0] == '3':
                c = int(b) - 2**15
                ninja.append(str(c))
            else:
                ninja.append(b)
        turtles.append(ninja)
        
    torts = ["_".join(t) for t in turtles]
    torts = [t.replace("__", "_") + "\n" for t in torts]
    with open(path_in, "wt") as pstr:
        pstr.truncate(0)
        for line in torts:
            pstr.write(line)


# Take a string and automatically save it.
# File_name is in format "run-num_org-id.seq".
def save_aaff_to_fasta(aaff_string, file_name):
    nt_seq = unspool_aaff_for_msa(aaff_string)
    path_file = os.path.join(path_fasta, file_name + ".seq")
    with open(path_file, "wt") as f:
        f.write(nt_seq)


# Load as org_id:aaff_string. Other metadata isn't strictly required, alr in book_of_life.
def load_strain_genome(path_sgen):
    # Open the text archive.
    with open(path_sgen, "rt") as f:
        sfgen = f.readlines()
    keys = [int(x.split(" ")[0]) for x in sfgen]
    values = [x.split(":")[1][:-2] for x in sfgen]
    sfgen_dict = dict(zip(keys, values))
    return sfgen_dict


# Change data storage format for strain_genome_015.
def fix_aaff_into_xascii(path_in):
    with open(path_in, "rt") as sgen:
        wrong = sgen.readlines()
    org_gnm_pairlist = [x.split(":") for x in wrong]
    orgid_list = [x[0] for x in org_gnm_pairlist]
    genome_list = [x[1] for x in org_gnm_pairlist]
    gxa_list = [''.join(convert_aaff_to_xascii(retrieve_aaff(gnm))) for gnm in genome_list]
    ogxa_pairlist = [':'.join([x, y]).replace("__", "_").replace("\n", "") + "\n" for x, y in zip(orgid_list, gxa_list)]

    with open(path_in, "wt") as pstr:
        pstr.truncate(0)
        for line in ogxa_pairlist:
            pstr.write(line)
