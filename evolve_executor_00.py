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
# This is a borrowed algorithm.
import GlobalAlignment
import genome_handler as geha
import track_phylogeny as tphy

# To adjust the dataframe appearance
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 20)
pd.set_option("display.width", 200)
pd.set_option('display.expand_frame_repr', False)

# ===== PATHS =====
# Edit a BATCH file to run the input and output Evolve files.
path_evodir = r"C:\Program Files (x86)\Evolve"
path_workdir = r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\\"
path_bat_evotemp = os.path.join(path_workdir, "evo_template.bat")
path_bat_ugenetemp = os.path.join(path_workdir, "ugene_template.bat")
# Eventually can adjust based on user input.
run_num = "026"
# Extract from filename?
run_name = "big_bang"
# path_rundir = os.path.join(path_workdir, "Runs", "Set_001", "Run_" + run_num + "_" + run_name)
path_rundir = os.path.join(path_workdir, "Runs", "Set_001", "Run_" + run_num)
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


# https://stackoverflow.com/questions/28730961/python-slicing-string-in-three-character-substrings
def pair_split(elm):
    return [elm[s:s + 2] for s in range(0, len(elm), 2) if len(elm[s:s + 2]) > 1]


# Function to convert from binary to base-4, so as to retain 2's complement.
binfour_dict = {"00": "0", "01": "1", "10": "2", "11": "3"}
b4_nt_dict = {"0": "A", "1": "T", "2": "G", "3": "C"}
def convert_binary_to_base4(bin_num):
    # subs = [bin_num[s:s+2] for s in range(0,len(bin_num),2) if len(bin_num[s:s+2]) > 1]
    subs = pair_split(bin_num)
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
        # # Separate "rowX" into "row X".
        # indiv_genome = indiv_genome.replace("row", "row ")
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
        b4_genome = geha.convert_kforth_to_base4(indiv_genome)
        # type(nucleotide_genome) = string, containing ATGC.
        nt_genome = geha.convert_base4_to_nucleotide(b4_genome)
        # # type(ab_genome) = string, containing AA-FF.
        aaff_genome = geha.convert_base4_to_aaff(b4_genome)
        aaff_string = geha.store_aaff(aaff_genome)

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
        lives_output = geha.get_organics_from_universe(path_output_phascii)
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


# Simplify the inputs.
# But need to specify the bol_in_df eventually.
def driver(base, target, gap=-1, match=1, mismatch=-1, debug=False):
    gb = geha.translate_aaff_to_kforth(bgen_df.loc[base, 'Genome']).split(' ')
    gt = geha.translate_aaff_to_kforth(bgen_df.loc[target,'Genome']).split(' ')
    GlobalAlignment.driver(gb, gt, gap, match, mismatch, debug)

# ===== EXECUTION =====
bgen_df = pd.read_csv(r"C:/Users/Julio Hong/Documents/LioHong/Evolve-Archives/bol_gen_010.csv", index_col="Unnamed: 0")
# Create a version without the genome col.
sbol_df = bgen_df.iloc[:,:-1]
# driver(7638, 7854, -1, 1, -1, debug=True)
# bgen_df.loc[find_ancestors(949,sbol_df,5),'Genome']

# Try to fit this into https://github.com/alife-data-standards/alife-std-dev-python/tree/master
digevo_df = bgen_df.copy()
digevo_df['ancestor_list'] = digevo_df.loc[:,['Sporelayer','Quickener']].values.tolist()
# # Remove columns.
# digevo_df = digevo_df.drop(columns=['Sporelayer','Quickener'])
# Set id as index.
digevo_df.set_index('Orgid')
digevo_df.ancestor_list = digevo_df.ancestor_list.apply(lambda x: list(set(x)))
# This works but head() doesn't show that.
digevo_df.ancestor_list = digevo_df.ancestor_list.apply(lambda x: [y for y in x])
digevo_df.rename(columns={'Orgid':'id', 'Birth_step':'origin_time', 'Death_step':'destruction_time', 'Genome':'sequence'}, inplace=True)
# digevo_df.sequence = digevo_df.sequence.apply(lambda x: retrieve_aaff(x),snum=True)

# from ALifeStdDev import phylogeny as asd_phylo
# deeeee = asd_phylo.pandas_df_to_networkx(digevo_df)
# tjhiorpra = asd_phylo.load_phylogeny_to_pandas_df(r"C:\Users\Julio Hong\Documents\LioHong\alife-std-dev-python-master\example_data\asexual_phylogeny_test.csv")
# digevo_df.to_csv(r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Archives\digevo_std_bgen010.csv")

phylogeny = digevo_df.loc[:200,'ancestor_list'].to_dict()

# https://colab.research.google.com/github/emilydolson/alife-phylogeny-tutorial/blob/main/perfect_tracking_final.ipynb#scrollTo=AQleBmbYENpC
# https://deap.readthedocs.io/en/master/api/tools.html?highlight=history#deap.tools.History
from deap import base
from deap import tools
toolbox = base.Toolbox()
history = tools.History()

import matplotlib.pyplot as plt
import networkx
from networkx.drawing.nx_pydot import graphviz_layout

def evalOneMax(individual):
    return sum(individual),
toolbox.register("evaluate", evalOneMax)

graph = networkx.DiGraph(phylogeny)
graph = graph.reverse()     # Make the graph top-down
colors = [toolbox.evaluate(phylogeny[i])[0] for i in graph]
pos = graphviz_layout(graph, prog="dot")
# networkx.draw(graph, node_color=colors, pos=pos)
networkx.draw_networkx(graph, node_color=colors, pos=pos)
plt.show()

evorgs = [tphy.Evorg(**kwargs) for kwargs in bgen_df.to_dict(orient='records')]
from phylotrackpy import systematics
# sys = systematics.Systematics(lambda Evorg: Evorg.Genome)
sys = systematics.Systematics(lambda Evorg: Evorg.__repr__(), True, True, False, False)
for e in evorgs[1:100]:
    e.taxon = sys.add_org(e)
    s_children = bgen_df[bgen_df.Sporelayer == e.Orgid].index
    q_children = bgen_df[bgen_df.Quickener == e.Orgid].index
    # Just ignore q_children for now.
    for s in s_children:
        evorgs[s].taxon = sys.add_org(s, e.taxon)
    if not e.Death_step:
        sys.remove_org(e.taxon)

if False:
# if True:
   runin_timestep = check_input_files(path_rundir)
   # simulate_universe(1000, express=True)
   # simulate_universe(100, runin_timestep, express=True)
   # simulate_universe(10000, runin_timestep, 1, express=True)
   # simulate_universe(10000, runin_timestep, 100, express=True)
   simulate_universe(10000, runin_timestep, 1000, express=True)

