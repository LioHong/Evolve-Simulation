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


# Function to convert from binary to base-4, so as to retain 2's complement.
binfour_dict = {"00": "0", "01": "1", "10": "2", "11": "3"}
atgc_dict = {"0": "A", "1": "T", "2": "G", "3": "C"}
def convert_binary_to_base4(bin_num):
    # Split into pairs.
    # https://stackoverflow.com/questions/28730961/python-slicing-string-in-three-character-substrings
    subs = [bin_num[s:s+2] for s in range(0,len(bin_num),2) if len(bin_num[s:s+2]) > 1]
    subs = [binfour_dict[x] for x in subs]
    return "".join(subs)

# Edit a BATCH file to run the input and output Evolve files.
path_evodir = r"C:\Program Files (x86)\Evolve"
# path_workdir = r"C:\Users\Lio Hong\Desktop\Evolve Simulation\\"
path_workdir = r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\"
path_template_bat = os.path.join(path_workdir, "evo_template.bat")
# Eventually can adjust based on user input.
# path_rundir = os.path.join(path_workdir, "Runs", "Run_001_big_bang")
path_rundir = os.path.join(path_workdir, "Runs", "Run_003_big_bang")
# Have to change workdir before the batch file can be successfully run.
os.chdir(path_rundir)
run_name = "big_bang"
path_start_evo = os.path.join(path_rundir, run_name + "_0")

# Convert KFORTH genome to DNA: 0123 for ATGC.
with open(os.path.join(path_workdir, "evodons_dict.txt"),'r') as filein:
    evodons_dict = eval(filein.read())
# Available range of integers: -131022 to 131022. (-2^17 - -49 to 2^17 - 49)
# Useful for generating dict and adjusting based on instructions, but loading from file would be preferable.
bitlength_evodons = 18
evodons_num_dict = {}
for posnum in range(0, 131023):
    evodons_num_dict[str(posnum)] = convert_binary_to_base4(format(posnum, "#020b")[2:])
for negnum in range(-1, -131023, -1):
    evodons_num_dict[str(negnum)] = convert_binary_to_base4(bin(negnum & (2**bitlength_evodons - 1))[2:])

# Set the time period.
time_period = 100
# For future formatting of filenames.
num_lead_zeroes = int(log10(time_period)) + 1
for timestep in range(1, time_period+1):
    # Progress update.
    if timestep % (max(time_period//10,1)) == 0:
        print(timestep)

    if timestep != 1:
        # Set the text to Find and Replace
        text_find_input = path_input_evo
        text_find_output = path_output_evo
        # Old output becomes input.
        path_input_evo = path_output_evo
        # Name the output file based on original simulation and time-step.
        path_output_evo = os.path.join(path_rundir, run_name + "_" + str(timestep))
        # Lives summary from output becomes that for input.
        lives_input = lives_output

        # # Old output becomes input.
        # text_find_input = path_output_evo
        # # Output updated.
        # path_output_evo = os.path.join(path_rundir, run_name + "_" + str(timestep))
        # # Set the text to Find and Replace
        # text_find_output = path_output_evo
    else:
        # Timestep equals 1, beginning of simulation.
        # Copy simulation from template.
        copyfile(os.path.join(path_rundir, "big_bang.evolve"), path_start_evo + ".evolve")
        path_input_evo = path_start_evo
        # Copy the bat file and rename it.
        path_run_bat = os.path.join(path_rundir, "run_003_big_bang.bat")
        copyfile(path_template_bat, path_run_bat)
        # Set the text to Find and Replace
        text_find_input = "path_in"
        text_find_output = "path_out"
        path_output_evo = os.path.join(path_rundir, run_name + "_" + str(timestep))
        # Assume that lives summary for timestep 0 always stays the same.
        path_input_phascii = path_input_evo + ".txt"
        lives_input = get_organics_from_universe(path_input_phascii)


    replaceds = {text_find_output: path_output_evo, text_find_input: path_input_evo}
    # print(replaceds)
    for rp in replaceds:
        with open(path_run_bat, "rt") as bat_file:
            text = bat_file.readlines()


        new_text = []
        for line in text:
            if rp in line:
                new_text.append(line.replace(rp, replaceds[rp]))
            else:
                new_text.append(line)

        with open(path_run_bat, "wt") as bat_file:
            bat_file.truncate(0)
            for line in new_text:
                bat_file.write(line)
    # print(text)

    # Run the bat.
    # Export PHASCII for output.
    # Get the filename itself to run.
    os.system(path_run_bat.split('\\')[-1])
    # Truncate the PHASCII to only ORGANIC onwards.
    path_output_phascii = path_output_evo + ".txt"
    # Open file as lines.
    with open(path_output_phascii, "rt") as phas:
        fas_text = phas.readlines()
    # Find the 'ORGANIC' entry and then slice the lines list.
    newfas_text = fas_text[fas_text.index("ORGANIC {\n"):]
    nffas_text = []
    # Replace "MAKE-SPORE" with "MAKE-SPOR" for later comparison.
    for line in newfas_text:
        if "MAKE-SPORE" in line:
            nffas_text.append(line.replace("MAKE-SPORE", "MAKE-SPOR"))
        else:
            nffas_text.append(line)
    nnfas_text = []
    # Replace "-" because it's part of negative integers and keywords.
    for line in nffas_text:
        if " - " in line:
            nnfas_text.append(line.replace(" - ", " ~ "))
        else:
            nnfas_text.append(line)
    # Overwrite the PHASCII.
    with open(path_output_phascii, "wt") as phas:
        phas.truncate(0)
        for line in nnfas_text:
            phas.write(line)

    # Extract genomes of all organisms in universe.
    with open(path_output_phascii, "rt") as phile:
        raw_phile = phile.readlines()
    org_flag = False
    popn_genome = {}
    for phline in raw_phile:
        if "ORGANISM" in phline:
            org_flag = True
            indiv_genome = []
        # Avoid adding the "CELL" line.
        if "CELL" in phline:
            org_flag = False
            # Save the stats of the organism.
            vital_stats = indiv_genome[0]
            indiv_genome = indiv_genome[1:]
            # Convert the genome into ATGC format.
            # Convert into single string.
            indiv_genome = ' '.join(indiv_genome)
            # Remove tabs.
            indiv_genome = indiv_genome.replace("    ", " ")
            # Remove multiple spaces.
            while "  " in indiv_genome:
                indiv_genome = indiv_genome.replace("  ", " ")
            # Handle occurrence of 'row'.
            indiv_genome = indiv_genome.replace("row", "row ")
            # Split into separate instructions and nums.
            indiv_genome = indiv_genome.split(" ")
            # Remove blank strings.
            while '' in indiv_genome:
                indiv_genome.remove('')

            # for evd in evodons_dict:
            #     indiv_genome = indiv_genome.replace(evd, evodons_dict[evd])
            # for evnum in evodons_num_dict:
            #     indiv_genome = indiv_genome.replace(evnum, evodons_num_dict[evnum])

            # When converting genomes from instructions and ints to base-pairs, convert instructions first.
            indiv_genome = [evodons_dict[evd] if evd in evodons_dict else evd for evd in indiv_genome ]
            indiv_genome = [evodons_num_dict[evnum] if evnum in evodons_num_dict else evnum for evnum in indiv_genome]
            indiv_genome = ' '.join(indiv_genome)
            # Remove all spaces.
            indiv_genome = indiv_genome.replace(" ", "")
            # Convert the list into a dict.
            popn_genome[vital_stats] = indiv_genome

            # Add a converter back from nucleotides to KFORTH.
            # Invert the dicts.
            # Slice the genome into evodons (length=9).
            # List comprehensions with if.

        if org_flag:
            # Replace.
            removables = ["\t", "\n", ":", "{", "}", "# program", '"']
            phresh = phline
            for rmvb in removables:
                phresh = phresh.replace(rmvb, "")
            # Add line.
            indiv_genome.append(phresh)

        # Convert to ATGC.
        # For spores, just find the lines between each index of spore, then the final entry.
        # Or pop 1st spore, then find 2nd spore.
        # Then slice the list up until 2nd spore.
        # Repeat until last spore popped, leaving the final entry?


    # Then convert ints.

    # Finally remove spaces.

    # Test which data format takes up the least memory: char-num, char-atgc, binary?

    # (Operation) Delete the old input evolve file.
    os.remove(path_input_evo + ".evolve")
    # Compare the ORGANISMS section of the input and output PHASCII files.
    lives_output = get_organics_from_universe(path_output_phascii)
    # # If the summaries are identical, then delete the output PHASCII (not lines in console).
    if lives_input == lives_output:
        os.remove(path_output_phascii)