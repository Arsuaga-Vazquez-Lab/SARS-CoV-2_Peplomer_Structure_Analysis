import re
import os
import modeller as m
from modeller.automodel import *

env = m.environ()

def align(target_name: str, target_sequence: str, template_name: str, template_chain: chr) -> None:
    # creates a file called f'alignment_{target_name}_and_{template_name}.pir'
    # assumes a file already exists called f'{template_name}.pdb'
    with open(f'{target_name}.pir', 'w+') as file:
        file.write(f'>P1;{target_name}\nsequence:{target_name}::::::::\n{target_sequence}*')
    alignment_instance = m.alignment(env)
    model_instance = m.model(env)
    model_instance.read(file=template_name, model_segment=(f'FIRST:{template_chain}', f'LAST:{template_chain}'))
    alignment_instance.append_model(model_instance,
                                    align_codes=template_name,
                                    atom_files=template_name)
    alignment_instance.append(file=f'{target_name}.pir', align_codes=target_name)
    alignment_instance.align2d()
    alignment_instance.write(file=f'alignment.pir')

def fold(target_name: str, target_sequence: str, template_name: str, template_chain: chr) -> None:
    alignment = align(target_name, target_sequence, template_name, template_chain)
    model = automodel(env,
                      alnfile=f'alignment.pir',
                      knowns=template_name,
                      sequence=target_name,
                      assess_methods=(assess.DOPE, assess.GA341))
    model.make()
    # Delete unnecessary files that that just made
    for filename in os.listdir():
        if re.search(f'^{target_name}.*', filename) and not re.search('.*\.pdb$', filename):
            os.remove(filename)
