""" Defines a function sequence_from_pdb that returns an amino acid sequence
of a given protein by searching through the PDB file

Example usage:

from pdb_to_sequence import sequence_from_pdb
pdb_file = open('6m0j_A.pdb')
print(sequence_from_pdb(pdb_file, chain='A'))

"""
from io import StringIO

amino_acid_letters = {
  # This info is taken from:
  # https://i-base.info/ttfa/hiv-and-drug-resistance/appendix-3-list-of-amino-acids-and-their-abbreviations/
  'ALA': 'A',
  'ARG': 'R',
  'ASN': 'N',
  'ASP': 'D',
  'CYS': 'C',
  'GLU': 'E',
  'GLN': 'Q',
  'GLY': 'G',
  'HIS': 'H',
  'ILE': 'I',
  'LEU': 'L',
  'LYS': 'K',
  'MET': 'M',
  'PHE': 'F',
  'PRO': 'P',
  'SER': 'S',
  'THR': 'T',
  'TRP': 'W',
  'TYR': 'Y',
  'VAL': 'V'
}

def sequence_from_pdb(pdb_filename: str, chain: chr) -> str:
    with open(pdb_filename) as file:
        lines = file.read().split('\n')
    atom_lines = [line for line in lines if line[0:4] == 'ATOM' and line[21] == chain]
    sequence = ''
    for atom in atom_lines:
        current_residue = amino_acid_letters[atom[17:20]]
        current_position = int(atom[22:26])
        if current_position > len(sequence) + 1:
            sequence += '-' * (current_position - len(sequence) - 1)
        if current_position == len(sequence) + 1:
            sequence += current_residue
        if current_residue == len(sequence) and sequence[-1] != current_residue:
            raise Exception(f'Error: multiple amino acids at position {current_position} of chain {chain}')
    # It's okay for a chain to start at a position other than residue 1
    return sequence.lstrip('-')
