#!/usr/bin/env python
import pdb_to_sequence
import modelling
import os
import csv
import calcs
import calculate_biochem_properties

wuhan_ref_seq = 'MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT'

AA_alphabet = list('GALMFWKQESPVICYHRNDT')
# 6m0j_A.pdb is a file I made by taking 6m0j.pdb and removing the B position of the rotamer at position 493 in the E chain
backbone_of_6m0j = calcs.backbone('6m0j_A.pdb', 'E')  # Used for calculating RMSD later on

mutated_sequences = {}  # dictionary of name to sequence
for substitution477 in AA_alphabet:
    new_seq = list(wuhan_ref_seq)
    new_seq[477-1] = substitution477
    for substitution501 in AA_alphabet:
        new_seq[501-1] = substitution501
        # According to 6m0j.pdb, the RBD is from residues 333-526, inclusive
        # Nicely, the AA sequence of the RBD in 6m0j.pdb perfectly matches the wuhan reference sequence
        new_rbd_seq = ''.join(new_seq[333-1:526])
        mutation_name = f'mutations_S477{substitution477}_and_N501{substitution501}'
        mutated_sequences[mutation_name] = new_rbd_seq

if not os.path.exists('mutations'):
    os.mkdir('mutations')
mutations_already_folded = [filename.split('.pdb')[0] for filename in os.listdir('mutations')]


# Now delete the csv file if it already exists, and create a new one
if os.path.exists('mutations/mutations.csv'):
    os.remove('mutations/mutations.csv')
with open('mutations/mutations.csv', 'w+') as file:
    csv_writer = csv.writer(file)
    fields = ['Name',
        'Sequence',
        'RMSD from 6M0J',
        'AA 477',
        'AA 501',
        'GRAVY @ 477',
        'Isoelectric point @ 477',
        'GRAVY @ 501',
        'Isoelectric point @ 501',
        'GRAVY',
        'Instability index',
        'Aromaticity',
        'Weight',
        'Isoelectric point',
        'Binding energy (calculated by Ahmed using Prodigy)',
        'Binding energy (calculated by Apurva using Prodigy)',
        'Binding energy (calculated by Nathan using Prodigy)',
    ]
    csv_writer.writerow(fields)
    for name, sequence in mutated_sequences.items():
        if name not in mutations_already_folded:
            modelling.fold(name, sequence, '6m0j', 'E')  # TODO: use 6m0j_A instead???
            print(f'Successfully folded {name}')
            # Modeller name folded pdb files weirdly and automatically puts them in the current directory
            os.rename(f'{name}.B99990001.pdb', f'mutations/{name}.pdb')
        alpha_carbons = calcs.backbone(f'mutations/{name}.pdb')
        rmsd_from_6m0j = calcs.rmsd(backbone_of_6m0j, alpha_carbons)
        # The RBD starts at residue 333
        aa477 = sequence[477 - 333]
        aa501 = sequence[501 - 333]
        biochemical_properties_of_aa477 = calculate_biochem_properties.biochemical_properties(aa477)
        biochemical_properties_of_aa501 = calculate_biochem_properties.biochemical_properties(aa501)
        biochemical_properties_of_overall_sequence = calculate_biochem_properties.biochemical_properties(sequence)
        csv_writer.writerow([
            name,
            sequence,
            rmsd_from_6m0j,
            aa477,
            aa501,
            biochemical_properties_of_aa477['GRAVY'],
            biochemical_properties_of_aa477['Isoelectric point'],
            biochemical_properties_of_aa501['GRAVY'],
            biochemical_properties_of_aa501['Isoelectric point'],
            biochemical_properties_of_overall_sequence['GRAVY'],
            biochemical_properties_of_overall_sequence['Instability'],
            biochemical_properties_of_overall_sequence['Aromaticity'],
            biochemical_properties_of_overall_sequence['Weight'],
            biochemical_properties_of_overall_sequence['Isoelectric point']
        ])
        print(f'Finished calcs for {name}.pdb')
