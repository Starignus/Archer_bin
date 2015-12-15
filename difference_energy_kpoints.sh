#!/bin/sh
#

#Difference from E(N)-E(last)
cat NVsEnergyRy.out | awk '{a[NR]=$1}END {for(i=1;i<NR+1;i++)print a[i]-a[NR]}' | sed '$d' > Difference_Ry_finalenergy.out
cat NVsEnergyeV.out | awk '{a[NR]=$1}END {for(i=1;i<NR+1;i++)print a[i]-a[NR]}' | sed '$d' > Difference_eV_finalenergy.out

###Difference from E(N-1)-E(N)
cat NVsEnergyRy.out | awk '{ print $1 }' | awk 'NR-1{print $0-p}{p=$0}' > Difference_Ry_energy.out
cat NVsEnergyeV.out | awk '{ print $1 }' | awk 'NR-1{print $0-p}{p=$0}' > Difference_eV_energy.out

#Difference wrt slab 10 Ang vacuum
cat NVsEnergyRy.out | awk '{a[NR]=$1}END {for(i=1;i<NR+1;i++)print a[4]-a[i]}'> Difference_Ry_energywrt18kmesh.out
cat NVsEnergyeV.out | awk '{a[NR]=$1}END {for(i=1;i<NR+1;i++)print a[4]-a[i]}'> Difference_eV_energywrt18kmesh.out
