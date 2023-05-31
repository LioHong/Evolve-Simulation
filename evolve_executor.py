# -*- coding: utf-8 -*-
"""
Filename: 
Date created: 2022/08/06, Sat, 12:00:34 (UTC+8)
@author: LioHong
Purpose: Automatically execute an Evolve file (universe) step-by-step for a defined time period and output PHASCII files.
Steps:
1. Set the time period
2. Edit a BATCH file to run the input and output Evolve files.
3. Name the output file based on original simulation and time-step.
4. Export PHASCII for old input.
5. (Operation) Delete the old input file.

"""
import os
from shutil import copyfile
from math import log10
from datetime import datetime
import pandas as pd
import numpy as np

# To adjust the dataframe appearance
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option("display.width", 200)
pd.set_option('display.expand_frame_repr', False)

# ===== PATHS =====
# Edit a BATCH file to run the input and output Evolve files.
path_evodir = r"C:\Program Files (x86)\Evolve"
path_workdir = r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\"
# path_workdir = r"C:\Users\Lio Hong\Documents\LioHong\Evolve-Simulation\\"
path_bat_evotemp = os.path.join(path_workdir, "evo_template.bat")
path_bat_ugenetemp = os.path.join(path_workdir, "ugene_template.bat")
# Eventually can adjust based on user input.
run_num = "026"
# Extract from filename?
run_name = "big_bang"
path_rundir = os.path.join(path_workdir, "Runs", "Run_" + run_num + "_" + run_name)
# Have to change workdir before the batch file can be successfully run.
os.chdir(path_rundir)
# Genome summary.
path_genome = os.path.join(path_rundir, "genomes_over_time_" + run_num + ".txt")
path_strain_genome = os.path.join(path_rundir, "strain_genome_" + run_num + ".txt")
path_book = os.path.join(path_rundir, "book_of_life_" + run_num + ".txt")
path_fasta = os.path.join(path_rundir, "Individual Genomes")
# All genomes present per timestep.
genomes_over_time = {}
# All genomes in the strain over time.
strain_genome = {}
book_of_life = {}


# ===== FUNCTIONS =====
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


def replace_old_with_new(path_input, old_new_dict):
    for oldnew in old_new_dict:
        with open(path_input, "rt") as f:
            text = f.readlines()

        new_text = []
        for line in text:
            if oldnew in line:
                new_text.append(line.replace(oldnew, old_new_dict[oldnew]))
            else:
                new_text.append(line)

        with open(path_input, "wt") as f:
            f.truncate(0)
            for line in new_text:
                f.write(line)


# Function to convert from binary to base-4, so as to retain 2's complement.
binfour_dict = {"00": "0", "01": "1", "10": "2", "11": "3"}
b4_nt_dict = {"0": "A", "1": "T", "2": "G", "3": "C"}
def convert_binary_to_base4(bin_num):
    # Split into pairs.
    # https://stackoverflow.com/questions/28730961/python-slicing-string-in-three-character-substrings
    subs = [bin_num[s:s+2] for s in range(0,len(bin_num),2) if len(bin_num[s:s+2]) > 1]
    subs = [binfour_dict[x] for x in subs]
    return "".join(subs)
# Convert KFORTH genome to DNA: 0123 for ATGC.
with open(os.path.join(path_workdir, "instr_b4_dict.txt"),'r') as filein:
    instr_b4_dict = eval(filein.read())
# AA-FF storage method.
with open(os.path.join(path_workdir, "b4_aaff_dict.txt"),'r') as filein:
    b4_aaff_dict = eval(filein.read())
    aaff_b4_dict = {v: k for k, v in b4_aaff_dict.items()}
# AA-FF storage method. Must have UTF-8 or else UnicodeDecodeError.
with open(os.path.join(path_workdir, "aaff_xascii_dict.txt"),'r', encoding='utf8') as filein:
    aaff_xascii_dict = eval(filein.read())
# Available range of integers: -131022 to 131022. (-2^17 - -49 to 2^17 - 49)
# Useful for generating dict and adjusting based on instructions, but loading from file would be preferable.
bitlength_evodons = 18
num_b4_dict = {}
for posnum in range(0, 131023):
    num_b4_dict[str(posnum)] = convert_binary_to_base4(format(posnum, "#020b")[2:])
for negnum in range(-1, -131023, -1):
    num_b4_dict[str(negnum)] = convert_binary_to_base4(bin(negnum & (2 ** bitlength_evodons - 1))[2:])

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
    # return str(dec_num) + "."
    return "_" + str(dec_num) + "_"


# Base-4 format used for compatibility with bits and ATGC.
# Consider whether to keep genomes as type list, then ''.join(list) only for the MSA and data storage?
# Most likely usage: Pick genomes from data storage, then convert them into nucleotides for MSA.
# Base-4 format used for compatibility with bits and ATGC.
def convert_kforth_to_base4(kforth_genome):
    # When converting genomes from instructions and ints to base-pairs, convert instructions first.
    b4_genome = [instr_b4_dict[instruction] if instruction in instr_b4_dict else instruction for instruction in kforth_genome]
    # Then convert ints.
    b4_genome = [num_b4_dict[keynum] if keynum in num_b4_dict else keynum for keynum in b4_genome]
    return b4_genome


def convert_base4_to_kforth(b4_genome):
    # Invert dicts
    b4_instr_dict = {v: k for k, v in instr_b4_dict.items()}
    br_num_dict = {v: k for k, v in num_b4_dict.items()}
    kforth_genome = [b4_instr_dict[evd_b4] if evd_b4 in b4_instr_dict else evd_b4 for evd_b4 in b4_genome]
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
    aaff_genome = [b4_aaff_dict[evd] if evd in b4_aaff_dict else evd for evd in b4_genome]
    # Converts to integer string.
    aaff_genome = [str(convert_b4_to_intstr(evd)) if len(evd) > 2 else evd for evd in aaff_genome]
    return aaff_genome


def store_aaff(aaff_genome):
    aaff_string = ''.join(aaff_genome)
    return aaff_string.replace("__","_")


def retrieve_aaff(aaff_string):
    # Separate out all the numbers.
    aaff_list = aaff_string.split("_")
    # Remove any empty strings.
    aaff_list = [x for x in aaff_list if x != ""]
    aaff_genome = []
    for elm in aaff_list:
        if not (elm.isdigit() or elm[0] == "-"):
            # Split all the remaining AAFFs into pairs.
            evd_list = [elm[s:s + 2] for s in range(0, len(elm), 2) if len(elm[s:s + 2]) > 1]
            # Flatten the list.
            for evd in evd_list:
                aaff_genome.append(evd)
        # Filter out all the numbers
        else:
            aaff_genome.append(elm)
    return aaff_genome


# Data de-compression.
def convert_aaff_to_base4(aaff_genome):
    b4_genome = [aaff_b4_dict[evd] if evd in aaff_b4_dict else evd for evd in aaff_genome]
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
    kforth_string = kforth_string.replace("row ", "row")    
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
  

# Convert string of numbers into organised df.          
# r"C:\Users\Lio Hong\Documents\LioHong\Evolve-Archives\book_of_life_010.txt"
# r"C:\Users\Lio Hong\Documents\LioHong\Evolve-Archives\strain_genome_010.txt"
# r"C:\Users\Lio Hong\Documents\LioHong\Evolve-Simulation\Runs\Run_015_big_bang\strain_genome_015a.txt"
# r"C:\Users\Lio Hong\Documents\LioHong\Evolve-Simulation\Runs\Run_015_big_bang\book_of_life_015.txt"
def examine_book_of_life(path_book):
    # Open the text archive.
    with open(path_book, "rt") as f:
        bol = f.readlines()
    # Format of an organism's entry: "31 1 1 1:[211, 387]/n"
    # Process the string into lists.
    idnums = [int(x.split(" ")[0]) for x in bol]
    sporelayers = [int(x.split(" ")[1]) for x in bol]
    sporeqks = [int(x.split(" ")[2]) for x in bol]
    generations = [int(x.split(" ")[3].split(":")[0]) for x in bol] 
    lifesteps = [x.split(":")[1][1:-2].split(",") for x in bol] 
    birthsteps = [int(x[0]) for x in lifesteps] 
    # For still-living organisms, set death-step to 0 so that lifespan will become negative.
    deathsteps = [int(x[1]) if len(x)>1 else 0 for x in lifesteps]
    # Combine all the lists into a df.
    zipped = zip(idnums, sporelayers, sporeqks, generations, birthsteps, deathsteps)
    bcols = ["ID", "Sporelayer", "Quickener", "Generation", "Birth_step", "Death_step"]
    blife_df = pd.DataFrame(zipped, columns=bcols)
    blife_df.set_index(["ID"], inplace=True)
    # Sort because organisms aren't added to the archive in order.
    blife_df = blife_df.sort_values(by=["ID"])
    blife_df["Lifespan"] = blife_df.Death_step - blife_df.Birth_step
    blife_df["Sex_check"] = blife_df.Quickener - blife_df.Sporelayer
    
    # ===== Data Handling =====
    # Check which organism lived the longest. 
    blife_df.Lifespan.max()
    # Remove negative lifespans used to rep living organisms. Can raise threshold.
    blife_df.loc[blife_df.Lifespan>-1,"Lifespan"].min()
    # Check how many organisms managed to reproduce.
    len(set(blife_df.Sporelayer))
    len(set(blife_df.Quickener))
    # Check which organism had the most children.
    blife_df.Sporelayer.value_counts()
    blife_df.Quickener.value_counts()
    # Check which organism produced the most children via sexual reproduction.
    blife_df.loc[blife_df.Sex_check!=0, "Sporelayer"].value_counts()
    blife_df.loc[blife_df.Sex_check!=0, "Quickener"].value_counts()
    # Find the oldest still-living organism.
    blife_df.loc[blife_df.Lifespan<0, "Birth_step"].min()
    # Related: Lifespan of oldest still-living organism.
    blife_df.loc[blife_df.Lifespan>0, "Lifespan"].max()

    return blife_df


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
    ogxa_pairlist = [':'.join([x,y]).replace("__","_").replace("\n","") + "\n" for x,y in zip(orgid_list, gxa_list)]
    
    with open(path_in, "wt") as pstr:
        pstr.truncate(0)
        for line in ogxa_pairlist:
            pstr.write(line)


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

       
# # For future formatting of filenames.
# num_lead_zeroes = int(log10(time_period)) + 1
def simulate_universe(time_period, runin_timestep=0, interval=1, express=False):
    # Packaged to ease readability of simulate_universe().
    def wrangle_genome(indiv_genome):
        # Save the stats of the organism while removing from genome.
        vital_stats = indiv_genome.pop(0)
        # Convert indiv_genome from a list of keywords into a single string.
        indiv_genome = ' '.join(indiv_genome)
        # Change of spacer: Convert tabs to spaces.
        indiv_genome = indiv_genome.replace("    ", " ")
        # Separate "rowX" into "row X".
        indiv_genome = indiv_genome.replace("row", "row ")
        # Remove multiple spaces.
        while "  " in indiv_genome:
            indiv_genome = indiv_genome.replace("  ", " ")
        # Split into separate instructions and nums using single spaces for translation of keywords only.
        indiv_genome = indiv_genome.split(" ")
        # Remove blank strings.
        # type(indiv_genome) = list.
        while '' in indiv_genome:
            indiv_genome.remove('')

        # type(evodon_genome) = string, containing 0123.
        b4_genome = convert_kforth_to_base4(indiv_genome)
        # type(nucleotide_genome) = string, containing ATGC.
        nt_genome = convert_base4_to_nucleotide(b4_genome)
        # # type(ab_genome) = string, containing AA-FF.
        aaff_genome = convert_base4_to_aaff(b4_genome)
        aaff_string = store_aaff(aaff_genome)

        # nucleotide_genome = convert_base4_to_nucelotide(evodon_genome)
        # popn_genome[vital_stats] = nucleotide_genome
        # Remove ENERGY and AGE from vital_stats.
        vs_list = vital_stats.split(" ")
        vital_stats = vs_list[1] + " " + " ".join(vs_list[4:7])
        # Can just keep adding genome repeatedly and it'll overwrite.
        # strain_genome[vital_stats] = nucleotide_genome
        popn_genome[vital_stats] = aaff_string
        strain_genome[vital_stats] = aaff_string

        return vital_stats, aaff_string, nt_genome, list(popn_genome.keys())
    
    # Performance metric.
    print("Scraper started at " + datetime.now().strftime("%H:%M:%S"))
    
    for timestep in range(runin_timestep, runin_timestep+time_period, interval):
        # Progress update. Adjust the frequency if time_period becomes larger?
        if timestep % (max(time_period//100,1)) == 0:
            print(timestep)
            print("Progress update at " + datetime.now().strftime("%H:%M:%S"))

        # Step-by-step initialisation.
        if timestep != runin_timestep:
            # Set the text to Find and Replace
            text_find_input = path_input_evo
            text_find_output = path_output_evo
            # Old output becomes input.
            path_input_evo = path_output_evo
            # Name the output file based on original simulation and time-step.
            path_output_evo = os.path.join(path_rundir, run_name + "_" + str(timestep+1))
            # Lives summary from output becomes that for input.
            lives_input = lives_output

        # Timestep equals 1, beginning of simulation. Initialise from template.
        else:
            # Copy EVOLVE and PHASCII from template?
            path_start_evo = os.path.join(path_rundir, run_name + "_" + str(runin_timestep))
            path_output_evo = os.path.join(path_rundir, run_name + "_" + str(runin_timestep+1))
            path_input_evo = path_start_evo
            # Copy the bat file and rename it.
            path_batrun_evo = os.path.join(path_rundir, "run_" + run_num + "_" + run_name + "_evolve.bat")
            copyfile(path_bat_evotemp, path_batrun_evo)
            # Set the text to Find and Replace
            text_find_input = "path_in"
            text_find_output = "path_out"


        # Run the batch file: Update paths in batch file based on timestep.
        replaceds = {text_find_output: path_output_evo, text_find_input: path_input_evo, "1u": str(interval)+"u"}
        replace_old_with_new(path_batrun_evo, replaceds)
        # Get the filename itself to run.
        os.system(path_batrun_evo.split('\\')[-1])

        # Export PHASCII for output: Extract only the ORGANIC section from the PHASCII.
        path_output_phascii = path_output_evo + ".txt"
        # Find the 'ORGANIC' entry and then slice the lines list.
        with open(path_output_phascii, "rt") as phas:
            fas_text = phas.readlines()
        newfas_text = fas_text[fas_text.index("ORGANIC {\n"):]
        with open(path_output_phascii, "wt") as phas:
            phas.truncate(0)
            for line in newfas_text:
                phas.write(line)
        # Update instructions in genome to avoid collisions during genome handling step.
        fix_collisions_dict = {"MAKE-SPORE": "MAKE-SPOR", " - ": " ~ ", "NUM-CELLS": "NUM-CELS"}
        replace_old_with_new(path_output_phascii, fix_collisions_dict)

        # Extract genomes of all organisms in universe.
        with open(path_output_phascii, "rt") as phile:
            raw_phile = phile.readlines()
        org_flag = False
        popn_genome = {}
        organisms_in_timestep = []
        for phline in raw_phile:
            # Add "ORGANISM" line via org_flag=True.
            if "ORGANISM" in phline:
                org_flag = True
                indiv_genome = []
            # Avoid adding the "CELL" line via org_flag=False.
            if "CELL" in phline and org_flag:
                org_flag = False
                vital_stats, indiv_genome, nucleotide_genome, population = wrangle_genome(indiv_genome)
                # Add organism to list of living.
                organisms_in_timestep.append(vital_stats)
                # Add birth-step of organism.
                if vital_stats not in book_of_life:
                    book_of_life[vital_stats] = [timestep]
            if org_flag:
                # Replace everything in indiv_genome that's not a keyword.
                removables = ["\t", "\n", ":", "{", "}", "# program", '"']
                phresh = phline
                for rmvb in removables:
                    phresh = phresh.replace(rmvb, "")
                # Add line.
                indiv_genome.append(phresh)
        # Add death-step of organism.
        for vs_org in book_of_life:
            # Check that organism was still alive.
            if len(book_of_life[vs_org]) < 2:
                # Check that organism is no longer listed among living organisms.
                if vs_org not in population:
                    # Add death-step. Subtract 1 to get final step while alive.
                    book_of_life[vs_org].append(timestep-1)
        
            # For spores, just find the lines between each index of spore, then the final entry.
            # Or pop 1st spore, then find 2nd spore.
            # Then slice the list up until 2nd spore.
            # Repeat until last spore popped, leaving the final entry?

        # Test which data format takes up the least memory: char-num, char-atgc, binary?

        # Operation: Delete the old input evolve file.
        os.remove(path_input_evo + ".evolve")
        # # Add population genome to genome tracking over time.
        # genomes_over_time[timestep] = popn_genome
        # Compare the ORGANISMS section of the input and output PHASCII files.
        lives_output = get_organics_from_universe(path_output_phascii)
        # # If the summaries are identical, then delete the output PHASCII (not lines in console).
        # if lives_input == lives_output and timestep != time_period:
        #     os.remove(path_output_phascii)
        #     # genomes_over_time.pop(timestep)
        # Include a mode which DELETES the intermediate PHASCIIs.
        # Would it be better not to create in the first place? But would require rewrite of the code.
        if timestep == time_period:
            express = False
        if express:
            try:
                os.remove(path_output_phascii)
            except:
                pass

    # Create another alternative dict that only records genomes of organisms. Takes up far less memory.
    file = open(path_strain_genome, "w")
    for key, value in strain_genome.items():
        file.write('%s:%s\n' % (key, value))
    file.close()
    
    # For phylogenetics, record book_of_life.
    file = open(path_book, "w")
    for key, value in book_of_life.items():
        file.write('%s:%s\n' % (key, value))
    file.close()

    # Automate archiving? Store run archive in Evolve-Archives, retain starting files and ending files.


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


# https://johnlekberg.com/blog/2020-10-25-seq-align.html
from itertools import product
from collections import deque
def needleman_wunsch(x, y, score_gap=-1, score_same=1, score_diff=-1):
    """Run the Needleman-Wunsch algorithm on two sequences.
    Code based on pseudocode in Section 3 of:
    Naveed, Tahir; Siddiqui, Imitaz Saeed; Ahmed, Shaftab.
    "Parallel Needleman-Wunsch Algorithm for Grid." n.d.
    https://upload.wikimedia.org/wikipedia/en/c/c4/ParallelNeedlemanAlgorithm.pdf
    """
    N, M = len(x), len(y)
    s = lambda a, b: int(a == b)
    LEFT = -1, 0
    UP = 0, -1
    DIAG = -1, -1

    # Create tables F and Ptr
    F = {}
    Ptr = {}

    F[-1, -1] = 0
    for i in range(N):
        F[i, -1] = -i
    for j in range(M):
        F[-1, j] = -j
        
    print(F)

    option_Ptr = DIAG, LEFT, UP
    for i, j in product(range(N), range(M)):
        option_F = (
            F[i - 1, j - 1] + s(x[i], y[j]),
            F[i - 1, j] - 1,
            F[i, j - 1] - 1,
        )
        F[i, j], Ptr[i, j] = max(zip(option_F, option_Ptr))

    print(F)
    print(Ptr)

    # Work backwards from (N - 1, M - 1) to (0, 0) to find the best alignment.
    alignment = deque()
    i, j = N - 1, M - 1
    while i >= 0 and j >= 0:
        direction = Ptr[i, j]
        if direction == DIAG:
            element = i, j
        elif direction == LEFT:
            element = i, None
        elif direction == UP:
            element = None, j
        alignment.appendleft(element)
        di, dj = direction
        i, j = i + di, j + dj
    while i >= 0:
        alignment.appendleft((i, None))
        i += score_gap
    while j >= 0:
        alignment.appendleft((None, j))
        j += score_gap

    return list(alignment)


def prrrint_align(x, y):
    alignment = needleman_wunsch(x,y)
    print("".join(
        "-" if i is None else x[i] for i, _ in alignment
    ))
    print("".join(
        "-" if j is None else y[j] for _, j in alignment
    ))

# # r.ggenealogy works with df containing 'child' and 'parent.
# # But getParent() only returns 1 value, even though it should return 2 values.    
# # Both the code for getChild() and getParent() are very similar: Selection of column in df.
# # So why can there be multiple children but not multiple parents?
# # I added more rows for the second parent and it worked.
bol_df = pd.read_csv(r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Archives\book_of_life_015_010.csv")
# onebol_df = bol_df.loc[:,['ID','Sporelayer']]
# onebol_df.rename(columns={'ID':'child', 'Sporelayer':'parent'}, inplace=True)
# twobol_df = bol_df.loc[:,['ID','Quickener']]
# twobol_df.rename(columns={'ID':'child', 'Quickener':'parent'}, inplace=True)
# # See how many 2nd parents there are.
sexbol_df = bol_df.loc[bol_df.Sex_check != 0]
# threebol_df = pd.concat([onebol_df,twobol_df.loc[bol_df.Sex_check != 0]])
# bbbol_df = pd.read_csv(r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Archives\bol_for_r.csv")
bbbol_df = pd.read_csv(r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Archives\bol_for_r.csv", index_col="Unnamed: 0")
# Actually there is no limit to the number of parents.


# ===== EXECUTION =====
bgen_df = pd.read_csv(r"C:/Users/Julio Hong/Documents/LioHong/Evolve-Archives/bol_gen_010.csv", index_col="Unnamed: 0")
# Create a version without the genome col.
sbol_df = bgen_df.iloc[:,:-1]

if False:
# if True:
   runin_timestep = check_input_files(path_rundir)
   # simulate_universe(1000, express=True)
   # simulate_universe(100, runin_timestep, express=True)
   # simulate_universe(10000, runin_timestep, 1, express=True)
   # simulate_universe(10000, runin_timestep, 100, express=True)
   simulate_universe(10000, runin_timestep, 1000, express=True)
