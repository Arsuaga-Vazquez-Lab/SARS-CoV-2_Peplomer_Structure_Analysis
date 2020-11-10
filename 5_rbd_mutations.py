"""Studies mutations in the 5 RBD positions specified in Cui et al (2019)
According to that paper, the 5 most important residue locations in the S protein are:
residues 442, 472, 479, 480, and 487
Ideally, we want to test all 20^5 possible mutations in those residues, and dock each mutated S protein
to the ACE2 receptor to test the binding affinity. Sadly, that would take a long long time.
Instead, we can sort the mutations based on how much their most important biochemical properties
differ from the wild type. That gives a prioritized list of mutations to fold and dock to ACE2,
which could help us find important mutations in those 5 residues"""
#TODO: consider alternative folding and docking software -- speed and accuracy are both important

import multiprocessing  # maybe use Slurm instead? idk how that works
import autodock  # look at Ryan's GA code for reference on how to use this
from typing import List
from io import StringIO  # Store the PDB of a protein as a StringIO in python. Distinguishes it from an ordinary string, but it's not a file, so it doesn't get written to the hard drive

def dock_to_ace2(rbd_pdb: StringIO) -> float:
    pass
    # return binding affinity

def fold(5_mutated_residues: List[chr]) -> StringIO:
    pass
    # return a StringIO of the folded protein's PDB file
    # use the code that's already written in modelling.py

def difference_from_wild_type(5_mutated_residues: List[chr]) -> float:
    """ returns a composite score of how different the biochemical properties of the 5 mutated residues
    are from the biochemical properties of those same 5 residues in the wild-type RBD
    Each property has a weight representing how important it is.
    For example, the ability to make hydrogen bonds is one of the most important
    If there are any red flags, so we know the mutations will fail, just return a very big number, like 999"""
    pass

# Begin by creating a giant list of all the mutations we want to test
possible_residues = 'ARNDCEQGHILKMFPSTWYV'  # double check this -- it should be a list of all 20 amino acids
mutations_to_test = []
for residue442 in possible_residues:
    for residue472 in possible_residues:
        for residue479 in possible_residues:
            for residue480 in possible_residues:
                for residue487 in possible_residues:
                    mutations_to_test.append([residue442, residue472, residue479, residue480, residue487])

# Next, weed out some mutations that we are CONFIDENT will fail
# if the difference from the wild type exceeds some threshold (TBD), it's not worth testing at all
# (dunno if that's best to do with a filter, mask, or list comprehension. consider all the options)

# Be sure to use a sorting algorithm that runs in O(nlogn) time
mutations_to_test.sort(key=difference_from_wild_type)
# As we gather more data, we could run some sort of regression to see which biochemical properties
# are most impactful on the binding affinity. Then we continue updating the prioritized list of mutations to test
# as we go by resorting it. This could complicate the code a bit though, and how we keep track of
# which mutations have already been tested

# Use some sort of parallelization to make this next part go faster
# Note: python has multithreading and multiprocessing, but I think these are different
for mutation in mutations_to_test:
    docking_score = dock_to_ace2(fold(mutation))
    
# TODO: finish this up
