#!/bin/sh
#
# Convergence test for charge density/dual dependence
# once is known the wavefunction cutoff
####################################################################

#### Specify these variables based on your environment
PSEUDO_DIR='/work/e05/e05/ablar/QE_PSEUDO'

#### System variables and files
SYSTEM=`pwd`
SYSTEM=${SYSTEM##*/}
VERSION=1.0
KPOINTS=../$SYSTEM'.kpoints'  #File with the kpoint mesh

###########Options###################
while getopts x:m:a:h:p:f:v:l:r: opt
do
   case "$opt" in
      x) BIN=$OPTARG;;
      m) MPI=$OPTARG;;
      a) NPOOL=$OPTARG;;
      h) NTG=$OPTARG;;
      p) PO="yes";;
      f) SEEDFILE=$OPTARG;;
      v) echo "resend version " $VERSION; exit 1;;
      l) BIN_PROJ=$OPTARG;;
      r) RMV="yes";;
      \?) echo """   resend usage:
    -x  path to the pw.x binary (if not in the \$PATH)
    -m  MPI driver with options (i.e. "mpirun -np 8")
    -a  Number of pools for scaling kpoints
    -h  3D FFT task groups
    -p  Plot Only. It assumes the calculation is already done previously
    -f  Input file without kmesh section
    -v  Print version (i.e. "yes")
    -l  Path to the projwfc.x binary (if not in the \$PATH)
    -r  Delete the folderes "results_" (i.e. -r "yes")
   Example cx1:Fe_Fe_vacuum_surface_layers.sh -m "mpiexec" -a "6" -h "3" -f "seed.input"  > OUT_vaccum_Uoccu_layer_Fe
        archer: Fe_Fe_vacuum_surface_layers.sh -m "aprun -n 336" -a "6" -h "3" -f "seed.input"  > OUT_vaccum_Uoccu_layer_F 
 """
     exit 1  ;;
   esac
done

echo "    Kpoints  test in Fe slab." $VERSION
echo "    -------------------" 

# Check for pw.x binary
PW_BIN=`which pw.x 2>/dev/null`

if [ "$BIN" ]; then
  PW_BIN=$BIN
fi

if [ -z "$PW_BIN" ]; then
   echo " # ERROR: Provide a path to a pw.x binary."
   exit
fi

echo "pw.x executed as: " $MPI $PW_BIN


# Number of pools
if [ "$NPOOL" ]; then
  POOL_C="-npool $NPOOL"
fi

if [ -z "$NPOOL" ]; then
   echo " # If you want to set the npool option  give a number."
   echo " If not leave it this way"
   POOL_C=" "

fi

if [ "$NTG" ]; then
    NTG_C="-ntg $NTG"
fi

if [ -z "$NTG" ]; then
   echo " # If you want to set the task grpups for scaling the 3d FFT  give an integer number."
   echo " If not leave it this way"
   NTG_C=" "

fi

echo "pw.x  with n pools executed as: " $MPI $PW_BIN "< input" $POOL_C  "> output"


# Check for projwfc.x binary
PROJWFC_BIN=`which projwfc.x 2>/dev/null`
if [ "$BIN_PROJ" ]; then  # Checking the path of projwfc.x
   PROJWFC_BIN=$BIN_PROJ
fi

if [ -z "$PROJWFC_BIN" ]; then
    echo " # ERROR: Provide a path to a projwfc.x binary."
    exit
fi

echo "projwfc.x executed as: " $MPI $PROJWFC_BIN




if [[ -z $PO ]]; then    # if PO option present skip this part

 #  Looop for the different files with different kmesh
  V="1" 

  while test 1 -eq `echo " $V < 6" | bc -l`; 
     do
     echo "Vaccum test  Surface: file"$V
     echo "  Temporary directry should be set manually in the input file" $SEEDFILE

  # Seed file
    if [ "$SEEDFILE" ]; then
       SEEDFILENEW=../$SEEDFILE
    fi
    if [ -z "$SEEDFILE" ]; then
       echo " # ERROR: Provide a path to a pw.x binary."
       exit
    fi
  
   # Folder where evrycase is running
     if [[ ! -d ${SYSTEM}_$V ]]; then
        mkdir ${SYSTEM}_$V
     fi 
     cd  ./${SYSTEM}_$V

    LATTICEVEC=../${SYSTEM}_$V.lattice_vectors
 

# self-consistent calculation
   echo 'Setting structure: Seed file do not have the lattice vectors'
   #Checking if the input doesn't already exist
   if [ ! -f ${SYSTEM}_$V.in ]; then

cat  $SEEDFILENEW $LATTICEVEC  >> ${SYSTEM}_$V.in
cat >> ${SYSTEM}_$V.in << EOF
  K_POINTS (automatic)
EOF
cat $KPOINTS >> ${SYSTEM}_$V.in
   fi #End of checking if inputfile areadyexist
 
     # Checking if the .out files already exist
     if [ ! -f ${SYSTEM}_$V.out ]; then
      $MPI $PW_BIN  $POOL_C $NTG_C -inp  ${SYSTEM}_$V.in  > ${SYSTEM}_$V.out 2> /dev/null
     fi # End checking .out exit

     #Saving energies and line of the kpath
     energiesry=$(grep -a "!    total energy              =" ${SYSTEM}_$V.out | awk 'END{print $5}' )
     eneriesev=`echo "$energiesry*13.605691" | bc -l`
     echo $K $energiesry >> ../NVsEnergyRy.out
     echo $K $eneriesev >> ../NVsEnergyeV.out

     #####Getting final cooridnates
     pwd
     CORRD=`grep 'Begin final coordinates'  ${SYSTEM}_$V.out`
     if [ -n "$CORRD" ]; then # relaxed coordinates Enters when string is non-zero
       echo "Checking for the final  structure"
       LASTSTEP=`grep !  ${SYSTEM}_$V.out | awk 'END{print}'`
       grep -A 150 "$LASTSTEP"  ${SYSTEM}_$V.out > TEMP'_fm_' # We have to change the number 150 if we have more than 30 atoms
       ############ Getting atomic coordinates #################################
       nat=`grep nat ${SYSTEM}_$V.in | awk '{ print substr($1,5,2)}'` #If the number of atoms is more of 2 ciphers or digits we most change the "2" to "3" 
       grep -A $nat 'ATOMIC_POSITIONS' TEMP'_fm_'  > ${SYSTEM}_${V}_final_atomic_possitions
       cp ${SYSTEM}_${K}_final_atomic_possitions  ${SYSTEM}_${V}_relax.xyz
       rm  COORD'_fm_' TEMP'_fm_'
       echo "Check the new input with the last structure generated"
     fi #end from "if" when starting looking for the relaxed coordinates 
  V=$[$V+1]  
  rm -r *.save*  
  cd ..
  done # loop while for calculations 
fi # end of PO section

# Delete resutls_ folders 

if [[ -n $RMV ]]; then
   echo $SYSTEM 
    for j in 1 2 3 4 5 
    do
       SYSTEMN=${SYSTEM}_$j;
       echo $SYSTEMN
       rm -r ${SYSTEMN}/${SYSTEMN}.save
    done
fi
echo " - DONE! -"

