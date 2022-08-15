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

# Edit a BATCH file to run the input and output Evolve files.
path_evodir = r"C:\Program Files (x86)\Evolve"
path_workdir = r"C:\Users\Lio Hong\Desktop\Evolve Simulation\\"
path_template_bat = os.path.join(path_workdir, "evo_template.bat")
# Eventually can adjust based on user input.
path_rundir = os.path.join(path_workdir, "Runs", "Run_001_big_bang")
# Have to change workdir before the batch file can be successfully run.
os.chdir(path_rundir)
run_name = "big_bang"
path_start_evo = os.path.join(path_rundir, run_name + "_0")

# Set the time period.
time_period = 1000
# For future formatting of filenames.
num_lead_zeroes = int(log10(time_period)) + 1
for timestep in range(1, time_period+1):
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
        path_run_bat = os.path.join(path_rundir, "run_001_big_bang.bat")
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
    nnfas_text = []
    # Replace "MAKE-SPORE" with "MAKE-SPOR" for later comparison.
    for line in newfas_text:
        if "MAKE-SPORE" in line:
            nnfas_text.append(line.replace("MAKE-SPORE", "MAKE-SPOR"))
        else:
            nnfas_text.append(line)
    # Overwrite the PHASCII.
    with open(path_output_phascii, "wt") as phas:
        phas.truncate(0)
        for line in nnfas_text:
            phas.write(line)

    # (Operation) Delete the old input evolve file.
    os.remove(path_input_evo + ".evolve")
    # Compare the ORGANISMS section of the input and output PHASCII files.
    lives_output = get_organics_from_universe(path_output_phascii)
    # # If the summaries are identical, then delete the output PHASCII (not lines in console).
    if lives_input == lives_output:
        os.remove(path_output_phascii)


# Convert KFORTH genome to DNA: 0123 for ATGC.



