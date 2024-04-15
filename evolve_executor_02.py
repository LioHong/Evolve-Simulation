# -*- coding: utf-8 -*-
"""
Filename: evolve_executor.py
Date created: 2024/03/05, Tue, 21:04:00 (UTC+8)
@author: LioHong
Purpose: Automatically execute an Evolve file (universe) step-by-step for a defined time period and output PHASCII files.
Steps:
1. Set the time period
2. Edit a BATCH file to run the input and output Evolve files.
3. Name the output file based on original simulation and time-step.
4. Export PHASCII for old input.
5. (Operation) Delete the old input file.

"""
from os import system
from shutil import copyfile
from pathlib import Path
from math import log10
from datetime import datetime
import pandas as pd
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
bat_tmpl_path = Path(".") / "evo_template_01.bat"
# Eventually can adjust based on user input.
grp_num =  "002"
run_num = "073"
# run_num = "034_prep"
# Extract from filename?
run_name = "smol02"
# run_name = "need_for_speed"
grp_dirpath = Path(".") / "Runs" / ("Grp_" + grp_num)
run_dirpath = grp_dirpath / ("Run_" + run_num)
rxiv_dirpath = Path.cwd().parent / "Evolve-Archives"
runxiv_dirpath = rxiv_dirpath / "Universes"/ "Grp_002" / ("Run_" + run_num)
# Genome summary.
strain_genome_path = run_dirpath / ("strain_genome_" + run_num + ".txt")
book_path = run_dirpath / ("book_of_life_" + run_num + ".txt")
bgen_path = run_dirpath / ("bgen_" + run_num + ".csv")
# cgen_path = run_dirpath / ("cgen_" + run_num + ".csv")
# cgd_path = run_dirpath / ("cgd_" + run_num + ".csv")
cgen_path = runxiv_dirpath / ("cgen_" + run_num + ".csv")
cgd_path = runxiv_dirpath / ("cgd_" + run_num + ".csv")
# All genomes present per timestep.
genomes_over_time = {}
# All genomes in the strain over time.
strain_genome = {}
book_of_life = {}


def replace_old_with_new(text, old_new_dict):
    for oldnew in old_new_dict:
        text = [line.replace(oldnew, old_new_dict[oldnew]) if oldnew in line else line for line in text]
    return text


# This function is used for variable inputs.
def check_input_files(run_path):
    # Read the filenames of the template files.
    grp_file_list = [x for x in run_path.iterdir() if x.is_file()]
    grp_file_list = [x.name.split('.')[0] for x in grp_file_list]
    # Check that the names of the EVOLVE and PHASCII both match.
    if len(set(grp_file_list)) > 1:
        print("Mismatch in EVOLVE and PHASCII detected.")
    # Then extract the timestep.
    else:
        print("Files matched.")


# Use if I'm testing a bunch of runs with the same universe and PHASCII.
def prep_new_run():
    grp_file_list = [x for x in grp_dirpath.iterdir() if x.is_file()]
    univ = [x for x in grp_file_list if '.evolve' in x.name][0]
    phascii = [x for x in grp_file_list if '.txt' in x.name][0]
    Path(run_dirpath).mkdir(exist_ok=True)
    # Check if template universe and PHASCII already exist.
    copyfile(univ, run_dirpath / univ.name)
    copyfile(phascii, run_dirpath / phascii.name)


# Mainly outputs strain_genome and book_of_life TXTs.
def glue_book(input_path, data_dict):
    lines = []
    for key, value in data_dict.items():
        lines.append('%s:%s\n' % (key, value))
    input_path.write_text("".join(lines), encoding="utf-8")

import cProfile
# # For future formatting of filenames.
# num_lead_zeroes = int(log10(time_period)) + 1
def simulate_universe(time_period, start_step=0, interval=1, delete=False, prep=False, speed=False):
    pr = cProfile.Profile()
    pr.enable()
    # Split and tidy biodata from a single organism.
    def wrangle_biodata(indiv_biodata):
        # Save the stats of the organism while removing from genome.
        vital_stats = indiv_biodata.pop(0)
        # Convert indiv_biodata from a list of keywords into a single string.
        indiv_biodata = ' '.join(indiv_biodata)
        # Change of spacer: Convert tabs to spaces.
        indiv_biodata = indiv_biodata.replace("    ", " ")
        # Remove multiple spaces.
        while "  " in indiv_biodata:
            indiv_biodata = indiv_biodata.replace("  ", " ")
        # Split into separate instructions and nums using single spaces for translation of keywords only.
        indiv_biodata = indiv_biodata.split(" ")
        # Remove blank strings.
        # type(indiv_biodata) = list.
        while '' in indiv_biodata:
            indiv_biodata.remove('')
        aaff_genome = geha.shrink_kforth_to_aaff(indiv_biodata)
        aaff_string = geha.store_aaff(aaff_genome)
        # Remove ENERGY and AGE from vital_stats.
        vs_list = vital_stats.split(" ")
        vital_stats = vs_list[1] + " " + " ".join(vs_list[4:7])
        # Can just keep adding genome repeatedly and it'll overwrite.
        popn_genome[vital_stats] = aaff_string
        strain_genome[vital_stats] = aaff_string
        return vital_stats, aaff_string, list(popn_genome.keys())
    
    # Performance metric.
    print("Scraper started at " + datetime.now().strftime("%H:%M:%S"))
    # Copy EVOLVE and PHASCII from template.
    if prep: prep_new_run()
    lives_output = ""
    for timestep in range(start_step, start_step+time_period, interval):
        # Progress update. Adjust the frequency if time_period becomes larger?
        if (timestep-start_step) % (max(time_period//100,1)) == 0:
            print(timestep)
            print("Progress update at " + datetime.now().strftime("%H:%M:%S"))

        # Step-by-step initialisation.
        if timestep != start_step:
            # Set the text to Find and Replace
            text_find_input = evin_fname
            text_find_output = evout_fname
            # Old output becomes input.
            evin_fname = evout_fname
            # Name the output file based on original simulation and time-step.
            evout_fname = run_name + "_" + str(timestep+interval)
            # Lives summary from output becomes that for input.
            lives_input = lives_output

        # Timestep equals 1, beginning of simulation. Initialise from template.
        else:
            evstart_path = run_name + "_" + str(start_step)
            evout_fname = run_name + "_" + str(start_step+interval)
            evin_fname = evstart_path
            # Copy the bat file and rename it.
            bat_run_path = run_dirpath / ("run_" + run_num + "_" + run_name + "_evolve.bat")
            copyfile(bat_tmpl_path, bat_run_path)
            # Set the text to Find and Replace
            text_find_input = "p_in"
            text_find_output = "p_out"

        # Run the batch file: Update paths in batch file based on timestep.
        br_replaceds = {text_find_output: evout_fname, text_find_input: evin_fname, "1u": str(interval)+"u"}
        text = bat_run_path.read_text(encoding="utf-8").splitlines()
        text = replace_old_with_new(text, br_replaceds)
        bat_run_path.write_text("\n".join(text), encoding="utf-8")
        # Get the filename itself to run.
        system(str(bat_run_path))
        # Export PHASCII for output: Extract only the ORGANIC section from the PHASCII.
        phas_path = run_dirpath / (evout_fname + ".txt")
        if not speed:
            # Find the 'ORGANIC' entry and then slice the lines list.
            phas = phas_path.read_text(encoding="utf-8").splitlines()
            phas = phas[phas.index("ORGANIC {"):]
            # Replace everything in indiv_biodata that's not a keyword.
            removables = ["\t", "\n", ":", "{", "}", "# program", '"']
            for rmvb in removables: phas = [r.replace(rmvb, "") for r in phas]
            # Update instructions in genome to avoid collisions during genome handling step.
            fix_collisions_dict = {"MAKE-SPORE": "MAKE-SPOR", " - ": " ~ ", "NUM-CELLS": "NUM-CELS"}
            phas = replace_old_with_new(phas, fix_collisions_dict)

            # Extract genomes of all organisms in universe.
            org_flag = False
            # Needed for death-step. Actually just need to copy orgids per timestep.
            popn_genome = {}
            organisms_in_timestep = []
            # When see "ORGANISM", start recording the genome.
            for phline in phas:
                if "ORGANISM" in phline:
                    org_flag = True
                    indiv_biodata = []
                # "CELL" indicates end of genome, so store collated biodata.
                if "CELL" in phline and org_flag:
                    org_flag = False
                    vital_stats, indiv_biodata, population = wrangle_biodata(indiv_biodata)
                    # Add organism to list of living.
                    organisms_in_timestep.append(vital_stats)
                    # Add birth-step of organism.
                    if vital_stats not in book_of_life:
                        book_of_life[vital_stats] = [timestep]
                # Line-by-line, record the genome.
                if org_flag:
                    # Add line.
                    indiv_biodata.append(phline)


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

            # # Add population genome to genome tracking over time.
            # genomes_over_time[timestep] = popn_genome
            # Compare the ORGANISMS section of the input and output PHASCII files.
            # Is this any use?
            lives_output = geha.get_organics_from_universe(phas)
            phas_path.write_text("\n".join(phas), encoding="utf-8")
            # # If the summaries are identical, then delete the output PHASCII (not lines in console).
            # if lives_input == lives_output and timestep != time_period:
            #     phas_path.unlink()
            #     # genomes_over_time.pop(timestep)
            # Include a mode which DELETES the intermediate PHASCIIs.
            # Would it be better not to create in the first place? But would require rewrite of the code.
            if (timestep-start_step) == time_period:
                print('(Timestep-start_step) equal to time period.')
                delete = False

        # Operation: Delete the old PHASCII.
        if delete:
            try: phas_path.unlink()
            except FileNotFoundError: print('FileNotFoundError but passing.')
        # Operation: Delete the old input evolve universe.
        (run_dirpath / (evin_fname + ".evolve")).unlink()

    # If not speed, outputs empty files.
    # strain_genome file: Stores genomes only.
    glue_book(strain_genome_path, strain_genome)
    # book_of_life file: Records parentage. Genealogy to phylogeny.
    glue_book(book_path, book_of_life)
    print("Scraper finished at " + datetime.now().strftime("%H:%M:%S"))
    pr.disable()
    pr.print_stats(sort='time')
    # Automate archiving? Store run archive in Evolve-Archives, retain starting files and ending files.


# Convert string of numbers into organised df.
def organise_book_of_life(book_path, save=True):
    # Open the text archive.
    with open(book_path, "rt") as f:
        bol = f.readlines()
    # bol = book_path.read_text(encoding="utf-8").splitlines()

    # Format of an organism's entry: "31 1 1 1:[211, 387]/n"
    # Process the string into lists.
    idnums = [int(x.split(" ")[0]) for x in bol]
    sporelayers = [int(x.split(" ")[1]) for x in bol]
    sporeqks = [int(x.split(" ")[2]) for x in bol]
    generations = [int(x.split(" ")[3].split(":")[0]) for x in bol]
    lifesteps = [x.split(":")[1][1:-2].split(",") for x in bol]
    birthsteps = [int(x[0]) for x in lifesteps]
    # For still-living organisms, set death-step to 0 so that lifespan will become negative.
    deathsteps = [int(x[1]) if len(x) > 1 else 0 for x in lifesteps]
    # Combine all the lists into a df.
    zipped = zip(idnums, sporelayers, sporeqks, generations, birthsteps, deathsteps)
    bcols = ["ID", "Sporelayer", "Quickener", "Generation", "Birth_step", "Death_step"]
    blife_df = pd.DataFrame(zipped, columns=bcols)
    blife_df.set_index(["ID"], inplace=True)
    # Sort because organisms aren't added to the archive in order.
    blife_df = blife_df.sort_values(by=["ID"])
    blife_df["Lifespan"] = blife_df.Death_step - blife_df.Birth_step
    blife_df["Sex_check"] = blife_df.Quickener - blife_df.Sporelayer
    if save: blife_df.to_csv(bgen_path)
    return blife_df


def examine_book_of_life(blife_df):
    # ===== Data Handling =====
    # Check which organism lived the longest.
    print('Longest-lived:')
    print(blife_df.Lifespan.max())
    # Remove negative lifespans used to rep living organisms. Can raise threshold.
    blife_df.loc[blife_df.Lifespan > -1, "Lifespan"].min()
    # Find the oldest still-living organism.
    print('Oldest living:')
    print(blife_df.loc[blife_df.Lifespan < 0, "Birth_step"].min())
    # Related: Lifespan of oldest still-living organism.
    print('Lifespan of oldest living:')
    print(blife_df.loc[blife_df.Lifespan > 0, "Lifespan"].max())
    # Check how many organisms managed to reproduce.
    print('Number of organisms who reproduced (S,Q):')
    print(len(set(blife_df.Sporelayer)))
    print(len(set(blife_df.Quickener)))
    # Check which organism had the most children.
    print('Most children:')
    print(blife_df.Sporelayer.value_counts().head())
    print(blife_df.Quickener.value_counts().head())
    # Check which organism produced the most children via sexual reproduction.
    print('Most children via sex:')
    print(blife_df.loc[blife_df.Sex_check != 0, "Sporelayer"].value_counts().head())
    print(blife_df.loc[blife_df.Sex_check != 0, "Quickener"].value_counts().head())


# Track contents of archive.
def record_archives():
    rxiv_path = Path.cwd().parent / 'Evolve-Archives'
    rlevels = ['*','*/*','*/*/*','*/*/*/*','*/*/*/*/*']
    rxiv_files = [str(rfile).split('Evolve-Archives')[1] for rl in rlevels for rfile in rxiv_path.glob(rl)]
    # rxiv_files = []
    # for rl in rlevels:
    #     for rfile in rxiv_path.glob(rl):
    #         rxiv_files.append(str(rfile))
    # rxiv_files = [r.split('Evolve-Archives')[1] for r in rxiv_files]
    rxiv_files.sort()
    (Path(".") / "archive_records.txt").write_text("\n".join(rxiv_files), encoding="utf-8")


# Another problem: How to stitch book_of_life and strain_genome from consecutive runs?
# Or I can just run this function several times then collate the combos?
def collate_books(run_nums_list,grp_num="002"):
    # Assume list of ints.
    run_nums_list.sort()
    r_nstr_list = [f"{x:03}" for x in run_nums_list]
    coll_b = pd.DataFrame()
    # Merge with the later run as priority.
    for r in reversed(r_nstr_list):
        print("Progress update at " + datetime.now().strftime("%H:%M:%S"))
        rdpath = Path(".") / "Runs" / ("Grp_" + grp_num) / ("Run_" + r)
        bkpath = rdpath / ("book_of_life_" + r + ".txt")
        sgpath = rdpath / ("strain_genome_" + r + ".txt")
        sdf = organise_book_of_life(bkpath, save=False)
        bdf = geha.stitch_sgen(sdf, sgpath)
        if not coll_b.empty:
            # coll_b = pd.concat([coll_b, bdf]).drop_duplicates()
            basket_df = pd.concat([coll_b, bdf])
            # Create new row with earlier birth_step and later death_step.
            dupes = basket_df[basket_df.index.duplicated()].index
            goodlives_df = coll_b.loc[dupes]
            goodlives_df['Birth_step'] = bdf.loc[dupes,'Birth_step']
            # Then re-calculate Lifespan.
            goodlives_df['Lifespan'] = goodlives_df['Death_step'] - goodlives_df['Birth_step']
            # Drop dupes and keep last.
            coll_b = pd.concat([basket_df, goodlives_df])
            coll_b = coll_b[~coll_b.index.duplicated(keep='last')].sort_index()
        else:
            coll_b = bdf.loc[:]
    return coll_b


# Simplify the inputs.
# But need to specify the bol_in_df eventually.
def driver(base, target, gap=-1, match=1, mismatch=-1, debug=False):
    gb = geha.translate_aaff_to_kforth(bgen_df.loc[base, 'Genome']).split(' ')
    gt = geha.translate_aaff_to_kforth(bgen_df.loc[target,'Genome']).split(' ')
    GlobalAlignment.driver(gb, gt, gap, match, mismatch, debug)


# Sanity check 1: Population census.
# All I want is to get list of orgids. It will be a messy function.
# Isn't this redoing the second half of simulate_universe()?
def snapshot_check(xpmt_df, check_p, toggle=True):
    tail = xpmt_df.tail(1).index[0]
    lngh = len(xpmt_df)
    ceil = max(tail, lngh)
    schk = [i for i in range(1,ceil+1) if i not in xpmt_df.index]
    # if tail > lngh:
    #     schk = [i for i in range(1,tail+1) if i not in xpmt_df.index]
    # else:
    #     schk = [i for i in range(1,lngh+1) if i not in xpmt_df.index]
    print(len(schk))

    # Generate reference from scratch using batch utility.
    # OR Provide input file, maybe from the 'parallel' runs.
    # I copied all universes from previous runs into a single folder.
    # Convert EVOLVE univs to PHASCIIs? Copy-paste code from simulate_universe().
    batr02_p = Path(".") / "evo_template_02.bat"
    rlevels = ['*']
    if toggle:
        exper_files = [rfile for rl in rlevels for rfile in check_p.glob(rl)]
        for chp in exper_files:
            # Get the filename itself to run.
            fname = chp.name.split('.')[0]
            batr_p = check_p / (fname + "_evolve.bat")
            br_rplc = {"p_in": fname, "p_out": fname}
            copyfile(batr02_p, batr_p)
            text = batr_p.read_text(encoding="utf-8").splitlines()
            text = replace_old_with_new(text, br_rplc)
            batr_p.write_text("\n".join(text), encoding="utf-8")
            system(str(batr_p))
    # Extract orgids and sgens from universe file.
    exper_files = [rfile for rl in rlevels for rfile in check_p.glob(rl)]
    phas_files = [ef for ef in exper_files if '.txt' in ef.name]
    chk_orgids = {}
    for pf in phas_files:
        phas = pf.read_text(encoding="utf-8").splitlines()
        # phas = pf.read_text(encoding="utf-8")
        phas = phas[phas.index("ORGANIC {"):]
        # # Replace everything in indiv_biodata that's not a keyword.
        # removables = ["\t", "\n", ":", "{", "}", "# program", '"']
        # for rmvb in removables: phas = [r.replace(rmvb, "") for r in phas]
        # # Update instructions in genome to avoid collisions during genome handling step.
        # fix_collisions_dict = {"MAKE-SPORE": "MAKE-SPOR", " - ": " ~ ", "NUM-CELLS": "NUM-CELS"}
        # phas = replace_old_with_new(phas, fix_collisions_dict)
        organisms = geha.get_organics_from_universe(phas, keys=["ORGANISM"])
        orgids = [int(x[1]) for x in organisms["ORGANISM"]]
        # Extract timestep.
        tstep = int(pf.name.split('.')[0].split('_')[1])
        chk_orgids[tstep] = orgids
        # Filter out orgids that are present in xpmt_df.
        cx = [org for org in orgids if org not in schk]
        # Find all orgids within experiment_df and check that they were indeed alive during that timestep.
        mini_df = xpmt_df.loc[cx]
        younglings = mini_df[mini_df.Birth_step > tstep].index
        ghosts = mini_df[mini_df.Death_step < tstep].index
        print(tstep)
        print('young')
        print(younglings)
        print('old')
        print(ghosts)

    # Also check their genomes.
    # Record how many were wrong.
    return chk_orgids


# Sanity check 2: Death check?
# Find all negative lifespans.

# ===== EXECUTION =====
record_archives()
# Create a version without the genome col.
# sbol_df = organise_book_of_life(book_path)
# bgen_df = geha.stitch_sgen(sbol_df, strain_genome_path)
# cgen_df, cgd = geha.compress_book(bgen_df, bgen_path, strain_genome_path, cgen_path, cgd_path, threshold=5)
# driver(7638, 7854, -1, 1, -1, debug=True)
# bgen_df.loc[find_ancestors(949,sbol_df,5),'Genome']

# digevo_df = tphy.fit_phylogeny(bgen_df)
# tphy.draw_phylogeny(digevo_df)

if False:
# if True:
   start_step = check_input_files(run_dirpath)
   # simulate_universe(1000, delete=True)
   # simulate_universe(100, start_step, delete=True)
   # simulate_universe(10000, start_step, 1, delete=True)
   # simulate_universe(10000, start_step, 100, delete=True)
   simulate_universe(10000, start_step, 1000, delete=True)

