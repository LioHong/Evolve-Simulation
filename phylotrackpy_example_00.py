# -*- coding: utf-8 -*-
"""
Filename: phylotrackpy_example.py
Date created: Sat Apr 13 15:48:43 2024
@author: Emily Dolson, Matthew Andres Moreno
Purpose: 
Steps: 
"""
import colorsys
from copy import copy
from dataclasses import dataclass

import alifedata_phyloinformatics_convert as apc
from Bio import Phylo as BioPhylo
from colormath import color_conversions, color_diff, color_objects
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm

from phylotrackpy import systematics

# # Reproducibility.
# np.random.seed(1)

# %load_ext watermark
# %watermark -iwbmuvg -iv

# # Patch numpy for compatibility with colormath package.
# https://github.com/gtaylor/python-colormath/issues/104
import numpy


def patch_asscalar(a):
    return a.item()


setattr(numpy, "asscalar", patch_asscalar)


# Write organism class.
@dataclass
class Organism:
    hue: float = 0.0
    saturation: float = 0.0
    value: float = 0.0

    def mutate(self: "Organism") -> None:
        """Probabilistically tweak stored color information."""
        if np.random.rand() < 0.5:
            self.hue = np.clip(self.hue + np.random.normal(0, 0.05), 0, 1)
        if np.random.rand() < 0.5:
            self.saturation = np.clip(
                self.saturation + np.random.normal(0, 0.02), 0, 1
            )
        if np.random.rand() < 0.5:
            self.value = np.clip(self.value + np.random.normal(0, 0.02), 0, 1)

    def make_offspring(self: "Organism") -> "Organism":
        """Return copy of self with mutation applied."""
        offspring = copy(self)
        offspring.mutate()
        return offspring

    def to_labcolor(self: "Organism") -> color_objects.LabColor:
        """Create colormath `LabColor` object representing stored color data."""
        as_hsv = color_objects.HSVColor(self.hue, self.saturation, self.value)
        return color_conversions.convert_color(as_hsv, color_objects.LabColor)

    def calc_distance(self: "Organism", other: "Organism") -> float:
        """Calculate color-theoretic distance between own color and other
        `Organism`'s color."""
        return color_diff.delta_e_cie1976(
            self.to_labcolor(), other.to_labcolor()
        )
    

# Calculate fitness values for population members, favoring Organisms that contrast other population members.
def calc_fitnesses(organisms: list[Organism]) -> list[float]:
    return [
        np.max(
            [
                Organism.calc_distance(organism, other)
                for other in np.random.choice(organisms, 10)
            ],
        )
        for organism in organisms
    ]


# Set up population tracking infrastructure. 
# Use dummy founder to force common ancestry among all population members.
# Write reproduce and remove helpers to wrap systematics bookkeeping tasks.
def reproduce(parent: Organism) -> Organism:
    offspring = parent.make_offspring()
    parent_taxon = taxa[id(parent)]
    taxa[id(offspring)] = syst.add_org(offspring, parent_taxon)
    return offspring


def remove(org: Organism) -> None:
    taxon = taxa[id(org)]
    syst.remove_org(taxon)
    del taxa[id(org)]
    
# Run rolling evolutionary loop, one tournament at a time.
population = [Organism() for _ in range(50)]

syst = systematics.Systematics(lambda org: str(org.hue))
founder_taxon = syst.add_org(Organism())
taxa = {id(org): syst.add_org(org, founder_taxon) for org in population}
syst.remove_org(founder_taxon)

TOURNAMENT_SIZE = 7
NUM_UPDATES = 2000


def main():
    for update in tqdm(range(NUM_UPDATES)):
        syst.set_update(update)  # track time in systematics manager
    
        # do one tournament
        fitnesses = calc_fitnesses(population)
        target_idx = np.random.randint(len(population))
        selection = max(
            np.random.randint(len(population), size=TOURNAMENT_SIZE),
            key=fitnesses.__getitem__,
        )
    
        # create offspring and replace target index
        offspring = reproduce(population[selection])
        remove(population[target_idx])
        population[target_idx] = offspring
    
    
def draw():
    # Output phylogenetic history, including column storing taxon info (hue values).
    phylo_p = r"C:\Users\Julio Hong\Documents\LioHong\Evolve-Simulation\Working Data\phylo.csv"
    syst.add_snapshot_fun(systematics.Taxon.get_info, "taxinfo")
    syst.snapshot(phylo_p)
    
    # Load from file and convert to BioPython.
    bp_tree = apc.alife_dataframe_to_biopython_tree(
        pd.read_csv(phylo_p),
        setattrs=["taxinfo"],
        setup_branch_lengths=True,
    )
    
    # Set branch colors and draw.
    for node in bp_tree.find_elements():
        if hasattr(node, "taxinfo"):
            rgb_float = colorsys.hsv_to_rgb(node.taxinfo, 1.0, 0.5)
            rgb_int = tuple(int(c * 255) for c in rgb_float)
            node.color = rgb_int
    
    with plt.rc_context({"lines.linewidth": 5}):
        BioPhylo.draw(bp_tree)
        
main()
draw()