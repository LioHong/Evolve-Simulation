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

# ===== PATHS =====
# Edit a BATCH file to run the input and output Evolve files.
path_evodir = r"C:\Program Files (x86)\Evolve"
path_workdir = r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\"
# path_workdir = r"C:\Users\Lio Hong\Documents\LioHong\Evolve-Simulation\\"
path_template_bat = os.path.join(path_workdir, "evo_template.bat")
# Eventually can adjust based on user input.
run_num = "014"
# Extract from filename?
run_name = "big_bang"
path_rundir = os.path.join(path_workdir, "Runs", "Run_" + run_num + "_" + run_name)
# Have to change workdir before the batch file can be successfully run.
os.chdir(path_rundir)
# Genome summary.
path_genome = os.path.join(path_rundir, "genomes_over_time_" + run_num + ".txt")
path_strain_genome = os.path.join(path_rundir, "strain_genome_" + run_num + ".txt")
path_book = os.path.join(path_rundir, "book_of_life_" + run_num + ".txt")
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
        dec_num = -2 ** 15

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
    b4_kword_dict = {v: k for k, v in kword_b4_dict.items()}
    br_num_dict = {v: k for k, v in num_b4_dict.items()}
    kforth_genome = [b4_kword_dict[evd_b4] if evd_b4 in b4_kword_dict else evd_b4 for evd_b4 in b4_genome]
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
    aaff_genome = [str(convert_b4_to_intstr(evd)) if len(evd) > 2 else evd for evd in aaff_genome]
    return aaff_genome


def store_aaff(aaff_genome):
    return ''.join(aaff_genome)


def retrieve_aaff(aaff_string):
    # Separate out all the numbers.
    aaff_list = aaff_string.split("_")
    # Remove any empty strings.
    aaff_list = [x for x in aaff_list if x != ""]
    aaff_genome = []
    for elm in aaff_list:
        if not elm.isdigit():
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
    aaff_b4_dict = {v: k for k, v in b4_aaff_dict.items()}
    b4_genome = [aaff_b4_dict[evd] if evd in aaff_b4_dict else evd for evd in aaff_genome]
    b4_genome = [num_b4_dict[keynum] if keynum in num_b4_dict else keynum for keynum in b4_genome]
    return b4_genome


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

    for timestep in range(runin_timestep, runin_timestep+time_period, interval):
        # Progress update. Adjust the frequency if time_period becomes larger?
        if timestep % (max(time_period//10,1)) == 0:
            print(timestep)

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
            path_run_bat = os.path.join(path_rundir, "run_" + run_num + "_" + run_name + ".bat")
            copyfile(path_template_bat, path_run_bat)
            # Set the text to Find and Replace
            text_find_input = "path_in"
            text_find_output = "path_out"


        # Run the batch file: Update paths in batch file based on timestep.
        replaceds = {text_find_output: path_output_evo, text_find_input: path_input_evo, "1u": str(interval)+"u"}
        replace_old_with_new(path_run_bat, replaceds)
        # Get the filename itself to run.
        os.system(path_run_bat.split('\\')[-1])

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
    

# ===== EXECUTION =====
runin_timestep = check_input_files(path_rundir)
# simulate_universe(1000, express=True)
# simulate_universe(100, runin_timestep, express=True)
# simulate_universe(10000, runin_timestep, 1, express=True)
simulate_universe(10000, runin_timestep, 100, express=True)