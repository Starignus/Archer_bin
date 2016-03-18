#!/usr/bin/env python

#
####################################################################
#
# Sum of Lodwin charges in the Slab: Total Charges, Substrate and adsorbate
#
####################################################################
#Author: Dr. Ariadna Blanca Romero
#        Postdoctoral Research Associate
#        Imperial College London
#        Thomas Young Centre-Chemestry
#        ariadna@starignus.com or starignus@gmail.com
#        https://github.com/Starignus
####################################################################

filename  = raw_input('QE_pdos.out file: ')
adsorbate = raw_input('If there is an adsorbate, True or False: ')
fe_atoms = int(raw_input('Number of Fe atoms in the Slab: '))
atoms = []
cherges_total = []
charges_Fe = []
charges_adsorb = []
with open(filename, "r") as f:
  for line in f:
    if "Lowdin Charges:" in line:
      break
  for line in f:
    line = line.strip() #removing spaces before and after line
    if not line:
      continue 
    if line.startswith("Atom #"):
      atom_num, values = line[len("Atom #"):].split(":", 1) 
      atom_num = int(atom_num)
      values = values.split(",") 
      values.pop()
      values = [value.split("=", 1) for value in values]
      values = [(label.strip(), float(value)) for label, value in values]
      values = dict(values)
      atoms.append(values)
      assert len(atoms) == atom_num


for i, atom in enumerate(atoms, 1):
  print i, atom
print

## Charge Sum
for i, charge in enumerate(atoms):
  charges = atoms[i]["total charge"]
  cherges_total.append(charges)

print "Total Lodwin charge in the system: ", sum(cherges_total)
charges_Fe = cherges_total[: fe_atoms]
print "Total Lodwin charge in Fe: ", sum(charges_Fe)

if adsorbate == "True":
  charges_adsorb = cherges_total [fe_atoms :]
  print "Total Lodwin charge in adsorbate: ", sum(charges_adsorb)
