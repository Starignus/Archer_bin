#!/bin/sh
####################################################################
#
# Convergence test for charge density when increasing nummber of layers
# once is known the wavefunction cutoff (No relaxation
# and nscf for PDOS and with a bigger value of smearing)
# 
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
SPECLIST=../$SYSTEM'.spec'
KPOINTS=../$SYSTEM'.kpoints'  #File with the kpoint mesh
KPOINTSDOS=../$SYSTEM'.kpoints_pdos' #File with the kpoint mesh for PDOS
#LATTICEVEC=../$SYSTEM'.lattice_vectors'
VERSION=1.0


###########Options###################
while getopts x:m:a:p:d:f:v:w:g:l:r:h: opt
do
   case "$opt" in
      x) BIN=$OPTARG;;
      m) MPI=$OPTARG;;
      a) NPOOL=$OPTARG;;
      h) NTG=$OPTARG;;
      p) PO="yes";;
      d) PDOS="yes";;
      f) CUTRHO=$OPTARG;;
      v) echo "resend version " $VERSION; exit 1;;
      w) ECUTWFC=$OPTARG;;
      g) FUNC=$OPTARG;;
      l) BIN_PROJ=$OPTARG;;
      r) RMV="yes";;
      \?) echo """   resend usage:
    -x  path to the pw.x binary (if not in the \$PATH)
    -m  MPI driver with options (i.e. "mpirun -np 8")
    -a  Number of pools for scaling kpoints
    -h  3D FFT task groups
    -p  Plot Only. It assumes the calculation is already done previously
    -d  Calculate the PDOS
    -f  Density cutoff energy in Ry
    -v  Print version (i.e. "yes")
    -w  Cuttoff energy obtained before in Ry.
    -g  Functional (PBE, PBEsol,etc)
    -l  Path to the projwfc.x binary (if not in the \$PATH)
    -r  Delete the folderes "results_" (i.e. -r "yes")
    -h  Surface arae in Ang^2
   Example:time Fe_Fe_converegence_surface_layers.sh -m "mpiexec" -a "6" -w "95" -f "800" -g "PBE" -d "yes" > OUT_FeS_FM_RUN_U_occu_layer110Fe (Running PDOS npool6 )
           time Fe_Fe_converegence_surface_layers.sh -m "mpiexec" -a "6" -h "2" -w "95" -f "800" -g "PBE" -d "yes" > OUT_FeS_FM_RUN_U_occu_layer110Fe (Running PDOS npool 6 ntg 2)
           time Fe_Fe_converegence_surface_layers.sh -m "mpiexec" -a "6" -w "95" -f "800" -g "PBE" > OUT_FeS_FM_RUN_U_occu_layer110Fe (Running just the last scf)
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




if [[ -z $PO ]]; then    # if PO option present skip this part

 #  Looop for the different number of layers
  N="2" #Number of layers (symmetric 110 )

  while test 1 -eq `echo " $N < 33" | bc -l`; 
     do
     echo "Number of Layers Surface (110)" $N
  
     # Coordinates
     COORDS=../${SYSTEM}_l$N.xyz
     LATTICEVEC=../${SYSTEM}_l$N.lattice_vectors
 
#     # Scratch directories for each 
#     if [[ ! -d $TMPDIR/results_$SYSTEM'_'$N ]]; then
#       mkdir  $TMPDIR/results_$SYSTEM'_'$N
#     fi

#     SCRATCHDIR=$TMPDIR/results_$SYSTEM'_'$N

#     echo "  Temporary directory= " $SCRATCHDIR
  
   # Folder where evrycase is running
     if [[ ! -d ${SYSTEM}_l$N ]]; then
        mkdir ${SYSTEM}_l$N
     fi 
     cd  ./${SYSTEM}_l$N
    
     SCRATCHDIR=./
     echo "  Temporary directory= " $SCRATCHDIR
 
   # Number of bands 
    bands=`echo "$N*16" | bc -l`

# self-consistent calculation

   #Checking if the input doesn't already exist
   if [ ! -f ${SYSTEM}_l$N.in ]; then

cat >  ${SYSTEM}_l$N.in << EOF
&control
    calculation='relax'
    restart_mode='from_scratch',    
    pseudo_dir = '$PSEUDO_DIR',
    outdir='$SCRATCHDIR'
    prefix='$SYSTEM',
    wf_collect = .true.
    tprnfor = .true.
    tstress=.true.
    nstep=50000
    etot_conv_thr=1.0D-5
    forc_conv_thr=1.0D-4
    disk_io='default'
    max_seconds=86000

 /
 &system
    ibrav=0,
    celldm(1)=1.889726878
    nat=$N ,
    ntyp= 1,
    ecutwfc = $ECUTWFC, ecutrho = $CUTRHO,
    occupations='smearing', smearing='mv', degauss=0.01,
    nspin=2, 
    starting_magnetization(1)= 0.6
    input_DFT='$FUNC'
    lda_plus_u = .true.,
    Hubbard_U(1) = 1.d-10,


 /
 &electrons
    diagonalization='david'
    electron_maxstep=1000
    conv_thr = 1.0e-8
    mixing_beta = 0.5
    mixing_mode='local-TF'
 /
 &ions
 /
 &cell
   press=0.0
   press_conv_thr=0.1
 /
ATOMIC_SPECIES
EOF
cat $SPECLIST $COORDS $LATTICEVEC >> ${SYSTEM}_l$N.in
cat >> ${SYSTEM}_l$N.in << EOF
  K_POINTS (automatic)
EOF
cat $KPOINTS >> ${SYSTEM}_l$N.in
   fi #End of checking if inputfile areadyexist
 
     # Checking if the .out files already exist
     if [ ! -f ${SYSTEM}_l$N.out ]; then
     $MPI $PW_BIN  $POOL_C $NTG_C -inp  ${SYSTEM}_l$N.in  > ${SYSTEM}_l$N.out  2> /dev/null
     fi # End checking .out exit

     #Saving energies and line of the kpath
     energiesry=$(grep -a "!    total energy              =" ${SYSTEM}_l$N.out | awk 'END{print $5}' )
     eneriesev=`echo "$energiesry*13.605691" | bc -l`
     echo $N $energiesry >> ../NVsEnergyRy.out
     echo $N $eneriesev >> ../NVsEnergyeV.out

     #####Getting final cooridnates
     pwd
     CORRD=`grep 'Begin final coordinates'  ${SYSTEM}_l$N.out`
     if [ -n "$CORRD" ]; then # relaxed coordinates Enters when string is non-zero
       echo "Checking for the final  structure"
       LASTSTEP=`grep !  ${SYSTEM}_l$N.out | awk 'END{print}'`
       grep -A 150 "$LASTSTEP"  ${SYSTEM}_l$N.out > TEMP'_fm_' # We have to change the number 150 if we have more than 30 atoms
       ############ Getting atomic coordinates #################################
       nat=`grep nat ${SYSTEM}_l$N.in | awk '{ print substr($1,5,2)}'` #If the number of atoms is more of 2 ciphers or digits we most change the "2" to "3" 
       grep -A $nat 'ATOMIC_POSITIONS' TEMP'_fm_'  > ${SYSTEM}_l${N}_final_atomic_possitions
       cp ${SYSTEM}_l${N}_final_atomic_possitions  ${SYSTEM}_l${N}_relax.xyz
       rm  COORD'_fm_' TEMP'_fm_'
       echo "Check the new input with the last structure generated"
       # Settig scratch directory for scf
       if [ ! -d results_${SYSTEM}_l${N}_dos ]; then
         mkdir  results_${SYSTEM}_l${N}_dos
       echo "Making"  results_${SYSTEM}_l${N}_dos
       fi
       SCRATCHDIRDOS=./results_${SYSTEM}_l${N}_dos
       echo "  Temporary directory="$SCRATCHDIRDOS 
       
       ##########  Building new input for scf #########################################

cat >  ${SYSTEM}_l${N}_scf.in << EOF
&control
    calculation='scf'
    restart_mode='from_scratch',    
    pseudo_dir = '$PSEUDO_DIR',
    outdir='$SCRATCHDIRDOS'
    prefix='$SYSTEM',
    wf_collect = .true.
    tprnfor = .true.
    tstress=.true.
    nstep=50000
    etot_conv_thr=1.0D-5
    forc_conv_thr=1.0D-4
    disk_io='default'
    max_seconds=86000

 /
 &system
    ibrav=0,
    celldm(1)=1.889726878
    nat=$N ,
    ntyp= 1,
    ecutwfc = $ECUTWFC, ecutrho = $CUTRHO,
    occupations='smearing', smearing='mv', degauss=0.01,
    nspin=2, 
    starting_magnetization(1)= 0.6
    input_DFT='$FUNC'
    lda_plus_u = .true.,
    Hubbard_U(1) = 1.d-10,


 /
 &electrons
    diagonalization='david'
    electron_maxstep=1000
    conv_thr = 1.0e-8
    mixing_beta = 0.5
    mixing_mode='local-TF'
 /
 &ions
 /
 &cell
   press=0.0
   press_conv_thr=0.1
 /
ATOMIC_SPECIES
EOF
cat $SPECLIST ${SYSTEM}_l${N}_relax.xyz $LATTICEVEC >> ${SYSTEM}_l${N}_scf.in
cat >> ${SYSTEM}_l${N}_scf.in << EOF
  K_POINTS (automatic)
EOF
cat $KPOINTS >>  ${SYSTEM}_l${N}_scf.in

         if [ ! -f ${SYSTEM}_l${N}_scf.out ];then
           $MPI $PW_BIN  $POOL_C $NTG_C -inp   ${SYSTEM}_l${N}_scf.in   >> ${SYSTEM}_l${N}_scf.out  2> /dev/null
          EF=`grep 'the Fermi energy' ${SYSTEM}_l${N}_scf.out | awk '{print $5}'`
          echo "EFermi  from scf=" $EF
          grep 'the Fermi energy' ${SYSTEM}_l${N}_scf.out | awk '{print $5}' >  ${SYSTEM}_l${N}_scf_EF
        else
          EF=`grep 'the Fermi energy' ${SYSTEM}_l${N}_scf.out | awk '{print $5}'`
          echo "EFermi  from scf=" $EF
          grep 'the Fermi energy' ${SYSTEM}_l${N}_scf.out | awk '{print $5}' >  ${SYSTEM}_l${N}_scf_EF
        fi # End checking .out exit

     else #from "if" when starting looking for the relaxe coordinates 
        echo "There are not final coordinates in the system with " $N "layers"  
     fi #end from "if" when starting looking for the relaxed coordinates 
  
  if [[ -n $PDOS ]]; then #If one wants to calculate the PDOS

    if [ -f ${SYSTEM}_l${N}_scf.out ]; then
# nscf calcualtion
cat >  ${SYSTEM}_l${N}_nscf.in << EOF
&control
    calculation='nscf'
    restart_mode='from_scratch',    
    pseudo_dir = '$PSEUDO_DIR',
    outdir='$SCRATCHDIRDOS'
    prefix='$SYSTEM',
    wf_collect = .true.
    tprnfor = .true.
    tstress=.true.
    nstep=50000
    etot_conv_thr=1.0D-5
    forc_conv_thr=1.0D-4
    disk_io='default'
    max_seconds=86000

 /
 &system
    ibrav=0,
    celldm(1)=1.889726878
    nat=$N,
    ntyp= 1,
    ecutwfc = $ECUTWFC, ecutrho = $CUTRHO,
    occupations='smearing', smearing='mv', degauss=0.01,
    nspin=2, 
    starting_magnetization(1)= 0.6
    input_DFT='$FUNC'
    lda_plus_u = .true.,
    Hubbard_U(1) = 1.d-10,


 /
 &electrons
    diagonalization='david'
    electron_maxstep=1000
    conv_thr = 1.0e-8
    mixing_beta = 0.5
    mixing_mode='local-TF'
 /
 &ions
 /
 &cell
   press=0.0
   press_conv_thr=0.1
 /
ATOMIC_SPECIES
EOF
cat $SPECLIST ${SYSTEM}_l${N}_relax.xyz $LATTICEVEC >> ${SYSTEM}_l${N}_nscf.in
cat >> ${SYSTEM}_l${N}_nscf.in << EOF
  K_POINTS (automatic)
EOF
cat $KPOINTSDOS >>  ${SYSTEM}_l${N}_nscf.in

      if [ ! -f ${SYSTEM}_l${N}_nscf.out ]; then
        $MPI $PW_BIN  $POOL_C $NTG_C -inp  ${SYSTEM}_l${N}_nscf.in   >> ${SYSTEM}_l${N}_nscf.out 2> /dev/null
        echo "Finishing nscf ..."
      fi

     ########## PDOS  ################
     echo "Post processing PDOS .."
cat > ${SYSTEM}_l${N}_pdos.in << EOF
 &projwfc
    prefix='$SYSTEM'
    outdir='$SCRATCHDIRDOS'
    filpdos='${SYSTEM}_l${N}'
    DeltaE=0.01
    ngauss=0
    degauss=0.01
    Emin=-50.0, Emax=50.0
 /
EOF

echo "  running PDOS calculation..."
      if [ ! -f ${SYSTEM}_l${N}_pdos.out ]; then
      $MPI $PROJWFC_BIN  $POOL_C -inp  ${SYSTEM}_l${N}_pdos.in  >> ${SYSTEM}_l${N}_pdos.out 2> /dev/null
      echo "Finishing PROJWFC  ..."
      fi
      echo "Eding  scf, nscf, and PDOS"

     cat ${SYSTEM}_l${N}.pdos_tot  | awk -v Ef="$EF" '{print $1-Ef, $2, $3*(-1)}' | sed '1d'  > ${SYSTEM}_l${N}.pdos.tot.shifted.dat
    fi # If exitst the file _scf.out it will run nscf, and pdos
  
  fi #End PDOS option  
  rm -r $SYSTEM.save *results_*
  rm *.igk* *.satw* *.wfc* *.atwf* 
  N=$[$N+2]   
  cd ..
  done # loop while for calculations 
fi # end of PO section


# Delete resutls_ folders 

if [[ -n $RMV ]]; then
   echo $SYSTEM 
   layers_list=`seq 2 2 32`
   array_layers=(${layers_list// / })
    for j in `seq 0 ${#arrya_layers[@]}`
    do
       rm -r ${SYSTEM}_l${j}/results_${SYSTEM}_l${j}_dos
       rm -r ${SYSTEM}_l${j}/$SYSTEM.save
    done
fi
echo " - DONE! -"

