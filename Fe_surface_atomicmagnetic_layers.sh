#!/bin/sh
#
#
####################################################################
#
# Script to get the magnetic moments of the atoms in each Fa slab
#
####################################################################
#Author: Dr. Ariadna Blanca Romero
#        Postdoctoral Research Associate
#        Imperial College London
#        Thomas Young Centre-Chemestry
#        ariadna@starignus.com or starignus@gmail.com
#        https://github.com/Starignus
####################################################################


#### System variables and files
SYSTEM=`pwd`
SYSTEM=${SYSTEM##*/}
VERSION=1.0


###########Options###################
while getopts x:m:f:z:g:h: opt
do
   case "$opt" in
      x) FROM=$OPTARG;;
      m) TO=$OPTARG;;
      f) SEEDFILE=$OPTARG;;
      z) PO="yes";;
      g) FROMDOS="yes";;
      h) SEEDPDOS=$OPTARG;;
      \?) echo """   resend usage:
    -x  From which # of layer you want to start (should be one after the smallest number of layers)
    -m  Up to how many layers
    -f  Inputfile (scf from relax extructure)
    -z  To skipe the analysis of magnetic moment from outpur scf
    -g  To calculate magnetic momnet (polarization) (yes)
    -h  Input file from PDOS calculation
   Example:(110)
           time Fe_surface_atomicmagnetic_layers.sh -x "1" -m "16" -f "Fe_16_sfc.out"
           time Fe_surface_atomicmagnetic_layers.sh -z "yes" -g "yes" -h "Fe_100_l2_pdos.out "

 """
     exit 1  ;;
   esac
done

echo "    Surfaceenergy  Layers (110,100,111) Fe slab." $VERSION
echo "    -------------------" 

if [[ -z $PO ]]; then
  layers_list=`seq $FROM 1 $TO`
  array_layers=(${layers_list// / })
  echo "List of layers" $layers_list
  for j in `seq 0 ${#array_layers[@]}`
    do
    magnetic_layer=`grep "atomic mag" $SEEDFILE |  tail -$TO | awk '{print $5 }'`
    array_magnetic_layer=(${magnetic_layer// / })
    echo ${array_layers[$j]} ${array_magnetic_layer[$j]} >> MagneticMomenteahlayer # Nlayres, M(Mbhor) 
  done

  grep 'magnetization' $SEEDFILE | tail -2 > TMadnAM
  cat TMadnAM | awk '{print $4}' 
  cat MagneticMomenteahlayer | awk '{print $2}'
  grep "Tr" $SEEDFILE | tail -$TO | awk '{print $8, $9, $10}' > OccupUpDnTotal #Occupation 
  cat OccupUpDnTotal
fi

if [[ -n $FROMDOS ]]; then  #If I want to get the data from PDOS
   grep "polarization = " $SEEDPDOS > PDOSpolarizationatomsall
   grep "polarization = " $SEEDPDOS | awk '{print $3}' > PDOSpolarizationtotalatomsatom
   grep "polarization = " $SEEDPDOS | awk '{print $12}' > PDOSpolarizationtotalatomsatom_dstates
cat PDOSpolarizationatomsall
cat PDOSpolarizationtotalatomsatom
cat PDOSpolarizationtotalatomsatom_dstates
fi
echo " - DONE! -" 

