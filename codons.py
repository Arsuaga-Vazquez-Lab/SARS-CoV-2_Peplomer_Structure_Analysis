import json

def print_json(thingy): print(json.dumps(thingy, indent=4))

# This table below is taken from https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes#SG1
base1       = list('TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG')
base2       = list('TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG')
base3       = list('TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG')
amino_acids = list('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')

# Instead of having nucleotides for each codon in 3 lists, combine them into one list
codons = [''.join(list(codon)) for codon in zip(base1, base2, base3)]
# Convert DNA codons to RNA codons by replacing 'T' with 'U'
codons =  [codon.replace('T', 'U') for codon in codons]
# Create dictionary of all AAs to a list of codons that can map for that AA
aa_to_codons = {}
codon_to_aa = {}
for i in range(64):
    aa_to_codons[amino_acids[i]] = aa_to_codons.get(amino_acids[i], [])
    aa_to_codons[amino_acids[i]].append(codons[i])
    codon_to_aa[codons[i]] = amino_acids[i]

def one_letter_away(word1: str, word2: str) -> bool:
    """ Given two words, return whether they are one point mutation away """
    if len(word1) != len(word2): raise ValueError('Strings must be the same length')
    return sum[letter1 != letter2 for letter1, letter2 in zip(word1, word2)] == 1
