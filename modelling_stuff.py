import re
import os
import modeller as m
from modeller.automodel import *
from io import StringIO

env = m.environ()

def align(target_name: str, target_sequence: str, template_name: str, template_chain: chr) -> StringIO:
    # creates a file called f'alignment_{target_name}_and_{template_name}.pir'
    # assumes a file already exists called f'{template_name}.pdb'
    target_pir = f'>P1;{target_name}\nsequence:{target_name}::::::::\n{target_sequence}*'
    target_pir = StringIO(target_pir)
    alignment_instance = m.alignment(env)
    model_instance = m.model(env)
    model_instance.read(file=template_name, model_segment=(f'FIRST:{template_chain}', f'LAST:{template_chain}'))
    alignment_instance.append_model(model_instance,
                                    align_codes=template_name,
                                    atom_files=template_name)
    alignment_instance.append(file=target_pir, align_codes=target_name)
    alignment_instance.align2d()
    alignment_filename = f'alignment_{target_name}_and_{template_name}.pir'
    alignment_instance.write(file=alignment_filename)
    # We don't actually care about the file, just its contents
    with open(alignment_filename) as file:
        file_contents = file.read()
    os.remove(alignment_filename)
    return StringIO(file_contents)

def fold(target_name: str, target_sequence: str, template_name: str, template_chain: chr) -> StringIO:
    alignment = align(target_name, target_sequence, template_name, template_chain)
    model = automodel(env,
                      alnfile=alignment,
                      knowns=template_name,
                      sequence=target_name,
                      assess_methods=(assess.DOPE, assess.GA341))
    model.make()
    # Delete unnecessary files that that just made
    for filename in os.listdir():
        if re.search(f'^{target_name}.*', filename) and not re.search('.*\.pdb$', filename):
            os.remove(filename)
    # Find the file that we want to return
    for filename in os.listdir():
        if re.search(f'^{target_name}\..*\.pdb$', filename):
            folded_pdb_file = open(filename)
            break
    folded_pdb_file_text = folded_pdb_file.read()
    folded_pdb_file.close()
    return StringIO(folded_pdb_file_text)
