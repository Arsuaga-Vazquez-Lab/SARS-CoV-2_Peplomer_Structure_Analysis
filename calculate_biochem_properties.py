""" All values calculated are in line with numbers calculated based on ExPASy numbers/scales

For documentation and more stuff to calculate, see:
https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html

If this fails, be sure to run the following commands in terminal:
pip install propy3
pip install biopython

Example usage:
>>> import json
>>> import calculate_biochem_properties
>>> sequence = 'qwertyipasdfgqqqqq'
>>> print(json.dumps(calculate_biochem_properties.biochemical_properties(sequence), indent=4))
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import MeltingTemp
from Bio.SeqUtils import lcc
from Bio.Seq import Seq
from propy import PyPro
from typing import Dict, Any

def h_bonding_percent(sequence):
    # Returns percentage of AAs in sequence that can form hydrogen bonds
    h_bonding_residues = [aa for aa in sequence if aa in 'CWNQSTYKRHDE']
    # Aboid dividing by zero here:
    if len(h_bonding_residues) == 0: return 0
    return 100 * len(h_bonding_residues) / len(sequence)

def melting_temp(sequence):
    try:
        return MeltingTemp.Tm_GC(sequence, strict=False)
    except ZeroDivisionError: return

def biochemical_properties(sequence: str) -> Dict[str, Any]:
    # Define objects used for calculations
    analysis_object = ProteinAnalysis(sequence)
    descriptor_object = PyPro.GetProDes(sequence)
    sequence_object = Seq(sequence)
    # TODO(Ahmed): Verify that all these calculations are actually returning reasonable values
    # For example, it says the percent composition of every amino acid is zero when I run
    # calculate_biochem_properties.biochemical_properties('qwertyipasdfghklcvnm')
    return {
        'Isoelectric point': analysis_object.isoelectric_point(),
        'Molecular weight': analysis_object.molecular_weight(),  # Daltons? Amu? g/mol?
        'Aromaticity': analysis_object.aromaticity(),
        'Instability index': analysis_object.instability_index(),
        'GRAVY': analysis_object.gravy(),
        'H-bonding percent': h_bonding_percent(sequence),
        'Melting temp': melting_temp(sequence),
        'LCC': lcc.lcc_simp(sequence)
    }


'''
Isoelectric Point --> Matches with https://web.expasy.org/compute_pi/ to 3 decimal places
Molecular Weight --> Is in Daltons and matches with https://web.expasy.org/compute_pi/
Aromaticity --> Calculates aromaticity value of a protein, reference:  https://biopython.org/docs/1.75/api/Bio.sequenceUtils.ProtParam.html#Bio.sequenceUtils.ProtParam.ProteinAnalysis.aromaticity
Instability Index --> Estimate of the stability of a protein in a test tube, above 40 is unstable. Matches with https://web.expasy.org/protparam/
GRAVY --> Gives a value describing the hydropathy of the protein, positive values are hydrophobic. Matches with https://web.expasy.org/protparam/
HBond ratio --> number of amino acids that can form hbonds divded by the total number of amino acids
'''
