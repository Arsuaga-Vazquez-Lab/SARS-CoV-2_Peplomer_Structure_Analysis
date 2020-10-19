#Will need modifications to automate and make more efficient but I think this should serve as a good basis.
#All values calculated are in line with numbers calculated based on ExPASy numbers/scales
from Bio.SeqUtils.ProtParam import ProteinAnalysis #BioPython Package
from propy import PyPro #ProPy 3 package, pip install propy3


z=input("Enter the name of your sequnce:\n")
a=input("Enter your sequence:\n")

p=ProteinAnalysis(a)#used for the calculations done by Bio.SeqUtils
DesObject=PyPro.GetProDes(a)#used for the calculations by PyPro

pI=(p.isoelectric_point())#calculates Isoelectric point
pI=("Isoelectric point="+str(round(pI,3)))
print(pI)

amino_acid_count=(p.count_amino_acids())#calculates the amino acid count by letter
aac=str(amino_acid_count)
print("Amino Acid count="+aac)
print("Total amount of amino acids in this sequence="+str(len(a)))

amino_acid_percent=DesObject.GetAAComp()#calculates the composition by percentages
aap=str(amino_acid_percent)
print("Amino Acid composition by percentage="+aap)

molec_weight=(p.molecular_weight())#calculates molecular weight
print("Molecular weight="+str(round(molec_weight,3)))

aromaticity=(p.aromaticity())#calculates aromaticity
print("Aromaticity="+str(round(aromaticity,3)))

instability_index=(p.instability_index())#instability_index
print("Instability Index="+str(round(instability_index,3)))

gravy=(p.gravy())#value for gravy
print("Gravy="+str(round(gravy,3)))

file=open(z,"w")
file.write('Sequence:'+z+'\n\n')

file.write('\n'+pI+'\n')

file.write('\n'+"Amino Acid count="+aac+'\n')

file.write(('\n'+"Total amount of amino acids in this sequence="+str(len(a))+"\n"))

file.write('\n'+"Amino Acid composition by percentage="+aap+'\n')

file.write('\n'+"Molecular weight="+str(round(molec_weight,3))+'\n')

file.write('\n'+"Aromaticity="+str(round(aromaticity,3))+'\n')

file.write('\n'+"Instability Index="+str(round(instability_index,3))+'\n')

file.write('\n'+"Gravy="+str(round(gravy,3)))

file.close()

print('Done!')
