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
from propy import PyPro
from typing import Dict, Any

def biochemical_properties(sequence: str) -> Dict[str, Any]:
    # Define objects used for calculations
    analysis_object = ProteinAnalysis(sequence)
    descriptor_object = PyPro.GetProDes(sequence)
    # TODO(Ahmed): Verify that all these calculations are actually returning reasonable values
    # For example, it says the percent composition of every amino acid is zero when I run
    # calculate_biochem_properties.biochemical_properties('qwertyipasdfghklcvnm')
    return {
        'Isoelectric point': analysis_object.isoelectric_point(),
        'AA counts': analysis_object.count_amino_acids(),  # Dict of letters to counts
        'AA percent compositions': descriptor_object.GetAAComp(),  # Gives the percent composition of each amino acid
        'Weight': analysis_object.molecular_weight(),  # Daltons? Amu? g/mol?
        'Aromaticity': analysis_object.aromaticity(),
        'Instability': analysis_object.instability_index(),
        'GRAVY': analysis_object.gravy()
    }
