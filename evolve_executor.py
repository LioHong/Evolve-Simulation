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
# path_workdir = r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\"
path_workdir = r"C:\Users\Lio Hong\Documents\LioHong\Evolve-Simulation\\"
path_template_bat = os.path.join(path_workdir, "evo_template.bat")
# Eventually can adjust based on user input.
# path_rundir = os.path.join(path_workdir, "Runs", "Run_001_big_bang")
run_num = "006"
path_rundir = os.path.join(path_workdir, "Runs", "Run_" + run_num + "_big_bang")
# Have to change workdir before the batch file can be successfully run.
os.chdir(path_rundir)
run_name = "big_bang"
path_start_evo = os.path.join(path_rundir, run_name + "_0")
# Genome summary.
path_genome = os.path.join(path_rundir, "genomes_over_time_" + run_num + ".txt")
path_strain_genome = os.path.join(path_rundir, "strain_genome_" + run_num + ".txt")
# Genomes per timestep.
genomes_over_time = {}
# Genomes per organism.
strain_genome = {}


# ===== FUNCTIONS =====
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
atgc_dict = {"0": "A", "1": "T", "2": "G", "3": "C"}
def convert_binary_to_base4(bin_num):
    # Split into pairs.
    # https://stackoverflow.com/questions/28730961/python-slicing-string-in-three-character-substrings
    subs = [bin_num[s:s+2] for s in range(0,len(bin_num),2) if len(bin_num[s:s+2]) > 1]
    subs = [binfour_dict[x] for x in subs]
    return "".join(subs)
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


# Add a converter back from nucleotides to KFORTH.
def nt_to_kforth(nt_seq):
    conv_evod = ''.join([{v: k for k, v in atgc_dict.items()}[base] for base in nt_seq])
    # Slice the genome into evodons (length=9).
    conv_evod = [conv_evod[s:s + 9] for s in range(0, len(conv_evod), 9) if len(conv_evod[s:s + 9]) > 8]
    # List comprehensions with if and inverted dicts.
    conv_evod = [{v: k for k, v in evodons_dict.items()}[evd] if evd
              in {v: k for k, v in evodons_dict.items()} else evd for evd in conv_evod]
    conv_evod = [{v: k for k, v in evodons_num_dict.items()}[evnum] if evnum
              in {v: k for k, v in evodons_num_dict.items()} else evnum for evnum in conv_evod]
    # Join with spaces.
    lang_genome = ' '.join(conv_evod)
    # Join row with row_number.
    lang_genome = lang_genome.replace("row ", "row")
    return lang_genome


# # For future formatting of filenames.
# num_lead_zeroes = int(log10(time_period)) + 1
def simulate_universe(time_period, express=False):
    # Packaged to ease readability of simulate_universe().
    def wrangle_genome(indiv_genome):
        # Base-4 format used for compatibility with bits and ATGC.
        def convert_keyword_to_evodon(indiv_genome):
            # When converting genomes from instructions and ints to base-pairs, convert instructions first.
            indiv_genome = [evodons_dict[keyword] if keyword in evodons_dict else keyword for keyword in indiv_genome]
            # Then convert ints.
            indiv_genome = [evodons_num_dict[keynum] if keynum in evodons_num_dict else keynum for keynum in indiv_genome]
            indiv_genome = ' '.join(indiv_genome)
            # Remove all spaces.
            indiv_genome = indiv_genome.replace(" ", "")
            return indiv_genome


        # Testing an alternative storage method.
        def convert_evodon_to_decimal(evd):
            # 1 & 0 are positive.
            if evd[0] == '1' or evd[0] == '0':
                dec_num = 2**16
            # 2 & 3 are negative.
            elif evd[0] == '2' or evd[0] == '3':
                dec_num = -2**17

            for i in range(1,len(evd)):
                real_num = int(evd[i]) * 4**(len(evd)-i-1)
                dec_num += real_num
            return dec_num


        # Full conversion to ATGC.
        def convert_base4_to_nucleotide(evodon_genome):
            # Convert to ATGC.
            nucleotide_genome = [atgc_dict[base] for base in evodon_genome]
            nucleotide_genome = ''.join(nucleotide_genome)
            return nucleotide_genome


        # Convert keywords to AA-FF for storage.


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
        evodon_genome = convert_keyword_to_evodon(indiv_genome)
        # type(nucleotide_genome) = string, containing ATGC.
        nucleotide_genome = convert_base4_to_nucleotide(evodon_genome)
        # # type(ab_genome) = string, containing AA-FF.
        # nucleotide_genome = convert_base4_to_nucelotide(evodon_genome)
        popn_genome[vital_stats] = nucleotide_genome
        # Remove ENERGY and AGE from vital_stats.
        vital_stats = " ".join(vital_stats.split(" ")[:-2])
        # Can just keep adding genome repeatedly and it'll overwrite.
        # strain_genome[vital_stats] = nucleotide_genome
        strain_genome[vital_stats] = [convert_evodon_to_decimal(evodon_genome[s:s + 9]) for s in range(0, len(evodon_genome), 9) if len(evodon_genome[s:s + 9]) > 8]

        return vital_stats, indiv_genome, nucleotide_genome

    for timestep in range(1, time_period+1):
        # Progress update. Adjust the frequency if time_period becomes larger?
        if timestep % (max(time_period//10,1)) == 0:
            print(timestep)

        # Step-by-step initialisation.
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

        # Timestep equals 1, beginning of simulation. Initialise from template.
        else:
            # Copy simulation from template.
            copyfile(os.path.join(path_rundir, run_name + ".evolve"), path_start_evo + ".evolve")
            path_input_evo = path_start_evo
            # Copy the bat file and rename it.
            path_run_bat = os.path.join(path_rundir, "run_" + run_num + "_" + run_name + ".bat")
            copyfile(path_template_bat, path_run_bat)
            # Set the text to Find and Replace
            text_find_input = "path_in"
            text_find_output = "path_out"
            path_output_evo = os.path.join(path_rundir, run_name + "_" + str(timestep))
            # Assume that lives summary for timestep 0 always stays the same.
            path_input_phascii = path_input_evo + ".txt"
            lives_input = get_organics_from_universe(path_input_phascii)

        # Run the batch file: Update paths in batch file based on timestep.
        replaceds = {text_find_output: path_output_evo, text_find_input: path_input_evo}
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
        fix_collisions_dict = {"MAKE-SPORE": "MAKE-SPOR", " - ": " ~ "}
        replace_old_with_new(path_output_phascii, fix_collisions_dict)

        # Extract genomes of all organisms in universe.
        with open(path_output_phascii, "rt") as phile:
            raw_phile = phile.readlines()
        org_flag = False
        popn_genome = {}
        for phline in raw_phile:
            # Add "ORGANISM" line via org_flag=True.
            if "ORGANISM" in phline:
                org_flag = True
                indiv_genome = []
            # Avoid adding the "CELL" line via org_flag=False.
            if "CELL" in phline:
                org_flag = False
                vital_stats, indiv_genome, nucleotide_genome = wrangle_genome(indiv_genome)
            if org_flag:
                # Replace everything in indiv_genome that's not a keyword.
                removables = ["\t", "\n", ":", "{", "}", "# program", '"']
                phresh = phline
                for rmvb in removables:
                    phresh = phresh.replace(rmvb, "")
                # Add line.
                indiv_genome.append(phresh)

            # For spores, just find the lines between each index of spore, then the final entry.
            # Or pop 1st spore, then find 2nd spore.
            # Then slice the list up until 2nd spore.
            # Repeat until last spore popped, leaving the final entry?

        # Test which data format takes up the least memory: char-num, char-atgc, binary?

        # Operation: Delete the old input evolve file.
        os.remove(path_input_evo + ".evolve")
        # Add population genome to genome tracking over time.
        genomes_over_time[timestep] = popn_genome
        # Compare the ORGANISMS section of the input and output PHASCII files.
        lives_output = get_organics_from_universe(path_output_phascii)
        # If the summaries are identical, then delete the output PHASCII (not lines in console).
        if lives_input == lives_output and timestep != time_period:
            os.remove(path_output_phascii)
            genomes_over_time.pop(timestep)
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

    # Automate archiving? Store run archive in Evolve-Archives, retain starting files and ending files.
    # For phylogenetics, record organisms_over_time.

# ===== EXECUTION =====
simulate_universe(100, express=True)