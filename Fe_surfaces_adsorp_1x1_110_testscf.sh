#!/bin/sh
#
####################################################################
#
# Script test to launch realaxtion of the Fe(100) surface. There are two
# files needed to start the excecution: The file.temp that contains
# the template of the input file for Quantum Espresso and the file.tempstart
# contains the the same data as file.temp with exeption of the atomic
# postions, lattice vectors, and k-mesh.
#
####################################################################
#Author: Dr. Ariadna Blanca Romero
#        Postdoctoral Research Associate
#        Imperial College London
#        Thomas Young Centre-Chemestry
#        ariadna@starignus.com or starignus@gmail.com
#        https://github.com/Starignus
####################################################################

#### Specify these variables based on your environment
PSEUDO_DIR='/work/e05/e05/ablar/QE_PSEUDO'

#### System variables and files
SYSTEM=`pwd`
SYSTEM=${SYSTEM##*/}
TEMPLATE=$SYSTEM'.temp' # Structure
TEMPLATESTART=$SYSTEM'.tempstart'  # Structure with the start
VERSION=1.0


###########Options###################
while getopts x:m:a:h:p:d:v:l:r: opt
do
   case "$opt" in
      x) BIN=$OPTARG;;
      m) MPI=$OPTARG;;
      a) NPOOL=$OPTARG;;
      h) NTG=$OPTARG;;
      p) PO="yes";;
      d) PDOS="yes";;
      v) echo "resend version " $VERSION; exit 1;;
      l) BIN_PROJ=$OPTARG;;
      r) RMV="yes";;
      \?) echo """   resend usage:
    -x  path to the pw.x binary (if not in the \$PATH)
    -m  MPI driver with options (i.e. "mpirun -np 8")
    -a  Number of pools for scaling kpoints
    -h  3D FFT task groups
    -p  If the relaxation is alredy done, it will skip the first part
    -d  Calculate the PDOS
    -v  Print version (i.e. "yes")
    -l  Path to the projwfc.x binary (if not in the \$PATH)
    -r  Delete the folderes "results_" (i.e. -r "yes")
    -h  Surface arae in Ang^2
   Example:time Fe_surfaces_adsorp_1x1.sh -m "mpiexec" -a "6" -d "yes" > OUT_FeS_FM_RUN_U_occu_layer110Fe (Running PDOS npool6 )
           time Fe_surfaces_adsorp_1x1.sh -m "mpiexec" -a "6" -h "2"  -d "yes" > OUT_FeS_FM_RUN_U_occu_layer110Fe (Running PDOS npool 6 ntg 2)
           time Fe_surfaces_adsorp_1x1.sh -m "mpiexec" -a "6"  > OUT_FeS_FM_RUN_U_occu_layer110Fe (Running just the last scf)
 """
     exit 1  ;;
   esac
done

echo "    Convergence Number of Layers (110) Fe slab." $VERSION
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
   echo " # If you want to set the npool option  give an iteger number which is a common divisor of the number of kpoints and cores."
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

echo "pw.x  with n pools executed as: " $MPI $PW_BIN $POOL_C $NTG_C "-inp input > output"

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

# Checkhing for .temp, .tempstart
if [[ ! -e $TEMPLATE ]] || [[ ! -e $TEMPLATESTART ]]; then
   echo " # One or more files for the input are missing ( .temp, .tempstart)"
   exit
fi



if [[ -z $PO ]]; then    # if PO option present skip this part

############# #Input_file
  #Setting oudir & prefix
  echo "Setting input file ..."
  sed "/outdir=/{s#.*#    outdir='./'#}" $TEMPLATE > ${SYSTEM}.in
  sed "/prefix=/{s#.*#    prefix='${SYSTEM}'#}" -i ${SYSTEM}.in
  sed "/calculation=/{s#.*#    calculation='scf'#}" -i ${SYSTEM}.in
  # Checking if the .out files already exist
  if [ ! -f ${SYSTEM}.out ]; then
     echo "Running:"
     $MPI $PW_BIN  $POOL_C $NTG_C -inp  ${SYSTEM}.in  > ${SYSTEM}.out  2> /dev/null
  fi # End checking .out exit

  pwd
cp ${SYSTEM}.in ${SYSTEM}_scf.in
   rm -r $SYSTEM.save
   rm *.igk* *.satw* *.wfc* *.atwf*
echo "Finish to relax and getting the relax structure"
fi #End of PO

echo "PDOS option:" $PDOS

if [[ -n $PDOS ]]; then   # if PDOS option present it calculates  the PDOS
   # Folder where evrycase is running
   echo "PDOS section active.."
   if [[ ! -d ${SYSTEM}_DOS ]]; then
        mkdir ${SYSTEM}_DOS
   fi
   cp ${SYSTEM}_scf.in ${SYSTEM}_DOS
   cd  ./${SYSTEM}_DOS

   # Settig scratch directory for scf
   if [ ! -d results_${SYSTEM}_dos ]; then
      mkdir  results_${SYSTEM}_dos
      echo "Making"  results_${SYSTEM}_dos
   fi
   SCRATCHDIRDOS=./results_${SYSTEM}_dos
   echo "  Temporary directory="$SCRATCHDIRDOS 
   sed "/outdir=/{s#.*#    outdir='$SCRATCHDIRDOS'#}" -i ${SYSTEM}_scf.in
   if [ ! -f ${SYSTEM}_scf.out ];then
      echo "Calculating scf..."
      $MPI $PW_BIN  $POOL_C $NTG_C -inp   ${SYSTEM}_scf.in   >> ${SYSTEM}_scf.out  2> /dev/null
      EF=`grep 'the Fermi energy' ${SYSTEM}_scf.out | awk '{print $5}'`
      echo "EFermi  from scf=" $EF
      grep 'the Fermi energy' ${SYSTEM}_scf.out | awk '{print $5}' >  ${SYSTEM}_scf_EF
   else
      EF=`grep 'the Fermi energy' ${SYSTEM}_scf.out | awk '{print $5}'`
      echo "EFermi  from scf=" $EF
      grep 'the Fermi energy' ${SYSTEM}_scf.out | awk '{print $5}' >  ${SYSTEM}_scf_EF
   fi # End checking .out exit
  
   if [ ! -f ${SYSTEM}_scf.out ];then
      ########## PDOS  ################
      echo "Post processing PDOS .."
      if [ ! -f ${SYSTEM}_pdos.in ]; then
      cat > ${SYSTEM}_pdos.in << EOF
 &projwfc
    prefix='$SYSTEM'
    outdir='$SCRATCHDIRDOS'
    filpdos='${SYSTEM}'
    DeltaE=0.01
    ngauss=0
    degauss=0.01
    Emin=-50.0, Emax=50.0
 /
EOF
      fi #end checking if ${SYSTEM}_pdos.in  alredy exist

      echo "  running PDOS calculation..."
      if [ ! -f ${SYSTEM}}_pdos.out ]; then
         $MPI $PROJWFC_BIN  $POOL_C -inp  ${SYSTEM}_pdos.in  >> ${SYSTEM}_pdos.out 2> /dev/null
         echo "Finishing PROJWFC  ..."
      fi
      echo "Eding  scf, nscf, and PDOS"

       cat ${SYSTEM}.pdos_tot  | awk -v Ef="$EF" '{print $1-Ef, $2, $3*(-1)}' | sed '1d'  > ${SYSTEM}.pdos.tot.shifted.dat
   fi # If exitst the file _scf.out it will run pdos
   cd ..
  
fi #End PDOS option  

# Delete resutls_ folders
 
echo "RMV option:" $RMV

if [[ -n $RMV ]]; then # If RMV then it deletes the files
   echo " Option to earease .save files" $SYSTEM 
   pwd
   rm -r $SYSTEM.save
   rm *.igk* *.satw* *.wfc* *.atwf* 
   rm -r ${SYSTEM}_DOS/results_${SYSTEM}_dos
   cd ${SYSTEM}_DOS
   rm  *.igk* *.satw* *.wfc* *.atwf*
   cd ..
fi
echo " - DONE! -"

