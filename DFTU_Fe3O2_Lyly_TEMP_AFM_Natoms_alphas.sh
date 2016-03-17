#!/bin/sh

####################################################
#
# Script to calculate self-consitently the Hubbard
# Parameter by the linear response method. 
#
# Response matrix is calcuated with the
# the resp_mat_2009.f90 (r.x)  
#
#Author: Dr. Ariadna Blanca Romero
#        Postdoctoral Research Associate
#        Imperial College London
#        Thomas Young Centre-Chemestry
#        ariadna@starignus.com or starignus@gmail.com
#        https://github.com/Starignus
####################################################

#### System variables and files
#PSEUDO_DIR='/work/ablancar/QE/PSEUDO'
PSEUDO_DIR='/work/e05/e05/ablar/QE_PSEUDO' # For Archer
SYSTEM=`pwd`
SYSTEM=${SYSTEM##*/}
COORDS_FE=$SYSTEM'_Fe.xyz'  #Fliel with the atomic coordinates in the crystaline system
COORDS_O=$SYSTEM'_O.xyz'  #Fliel with the atomic coordinates in the crystaline system
SPECLISTFE=$SYSTEM'_Fe.spec' #File with the atomic specie, atomic mass, and pseudopotential
SPECLISTO=$SYSTEM'_O.spec' #File with the atomic specie, atomic mass, and pseudopotential
KPOINTS=$SYSTEM'.kpoints'  #File with the kpoint mesh
LATTICEVECT=$SYSTEM'.lattice_invect' #Check this file is correct i.e. already in bohrs and multiplied by alat
TMPDIR=./
###########Options###################
while getopts x:y:a:p:d:m:g:w:f:n:s:z:r: opt
do
   case "$opt" in
      x) BIN=$OPTARG;;
      y) RESP=$OPTARG;;
      a) NPOOL=$OPTARG;;
      p) PO="yes";;
      d) LRPU="yes";;
      m) MPI=$OPTARG;;
      g) FUNC=$OPTARG;;
      w) ECUTWFC=$OPTARG;;
      f) ECUTRHO=$OPTARG;;
      n) NPA=$OPTARG;;
      s) DIFFATOM=$OPTARG;;
      z) ALPHA+=("$OPTARG");;
      r) RMV="yes";;
      \?) echo """   resend usage:
    -x  path to the pw.x binary (if not in the \$PATH)
    -y  path to the r.x binary (if not in the \$PATH)
    -a  Number of npools
    -p  It assumes the non-perturbed calculation was already done (i.e. -p "yes")
    -d  If present it does not calculate the response matrix with r.x (i.e. -d "yes")
    -m  MPI driver with options (i.e. "mpirun -np 8")
    -g  Functional: PBE,PBEsol, etc.     -w  Cutoff energy obtained before in Ry.
    -w  Ecutoff 
    -f  Cutoffrho.
    -n  Number of atoms in the system
    -s  Number of type of aotms (ntyp)
    -z  List of perturbation alpha
    -r  Delete the folderes "results_" (i.e. -r "yes")
    -h  print this help 
   Example:time DFTU_Fe_TEMP.sh -m "mpiexec" -w 95 -f 800 -g "PBE" -n "2" -s "4" -z "-0.05 0.00 0.05" (cx1)
           
 """
     exit 1  ;;
   esac
done

echo "    DFT+U Fe"
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

# Chcking if we have NPOOLS
if [ "$NPOOL" ]; then
  POOL_C="-npool $NPOOL"
fi

if [ -z "$NPOOL" ]; then
   echo " # If you want to set the npool option  give a number."
   echo " If not leave it this way"
   POOL_C=" "
fi

echo "pw.x  with n pools executed as: " $MPI $PW_BIN "< input" $POOL_C  "> output"

# Cheching for .spect, .xyz. .kpoints .kpoint_pdos
if [[ ! -e $SPECLISTFE ]] || [[ ! -e $SPECLISTO ]] || [[ ! -e $COORDS_FE ]] || [[ ! -e $COORDS_O ]] || [[ ! -e $LATTICEVECT ]] || [[ ! -e $KPOINTS ]]; then
   echo " # One or more files for the input are missing ( .spec, .xyz. .lattice_invect)"
   exit
fi

#############DFT+U Section ############
NPA2=$(( $NPA * 2))

############ Input for response matirx: resp_mat.f90 -> r.x  #########

# Check for r.x binary
RESP_BIN=`which r.x 2>/dev/null`

if [ "$RESP" ]; then
  PW_BIN=$RESP
fi

if [ -z "$RESP_BIN" ]; then
   echo " # ERROR: Provide a path to a r.x binary."
fi

# Settig outdir directory (cx1)
OUTDIR=$TMPDIR/DFTU/$SYSTEM

 if [ ! -d  $OUTDIR ]; then
   mkdir -p $OUTDIR
 fi

rm -rf $OUTDIR/*

if [[ -z $PO ]]; then    # if PO option present skip this part

########### Calculating unperturbed system
# Check if your label in the atomic positions (Crystal, Angstroms,etc) is correct. 
# Also check if you used ibrav or ibrav=0 (in this case you need to add the lattice vectors)
# 5.42 exp conventional lattice parameter, 5.3526 from DFT calculations (PAW)

cat > ${SYSTEM}.unperturbed_scf_Fe.in << EOF
&control
    calculation='scf'
    verbosity = 'high'
    restart_mode='from_scratch'
    prefix='$SYSTEM'
    pseudo_dir = '$PSEUDO_DIR'
    outdir='$OUTDIR/'
    wf_collect = .true.
    tstress=.true.
    tprnfor = .true.
    etot_conv_thr=1.0D-4
    forc_conv_thr=1.0D-3
    nstep=50000
 /
&system
    ibrav=0,
    celldm(1)=1.889726878,
    nat= $NPA, ntyp = $DIFFATOM, nspin = 2
    nbnd = 70,
    ecutwfc = $ECUTWFC, ecutrho = $ECUTRHO,
    occupations='smearing', smearing='mv', degauss=0.01,
    starting_magnetization(1)=0.6
    starting_magnetization(2)=-0.6
    starting_magnetization(4)=0.6
    lda_plus_u = .true.,
    U_projection_type = 'atomic',
    Hubbard_U(1)= 1.d-20
    Hubbard_U(2)= 1.d-20
    Hubbard_U(3)= 1.d-20
    Hubbard_U(4)= 1.d-20
 /
&electrons
    mixing_beta = 0.3 
    diagonalization='david'
    electron_maxstep=1000
    conv_thr =  1.0d-9
 /
 &ions
 /
 &cell
   press=0.0
   press_conv_thr=0.1
 /

ATOMIC_SPECIES
EOF
cat $SPECLISTFE $COORDS_FE $LATTICEVECT >> ${SYSTEM}.unperturbed_scf_Fe.in
cat >> ${SYSTEM}.unperturbed_scf_Fe.in << EOF
K_POINTS (automatic)
EOF
cat $KPOINTS >> ${SYSTEM}.unperturbed_scf_Fe.in

  # Checking if the .out files already exist
  if [ ! -f ${SYSTEM}.unperturbed_scf_Fe.out ]; then
     $MPI $PW_BIN < ${SYSTEM}.unperturbed_scf_Fe.in  $POOL_C >> ${SYSTEM}.unperturbed_scf_Fe.out 2> /dev/null
  fi # End checking .out exit

  echo "Ending unperturbed calculation...i"
fi # end of PO section

  
####Copying Potential from OUTDIR
UNPETURBPOT=UPOTENTIAL
if [ ! -d  $UNPETURBPOT ]; then
   mkdir -p $UNPETURBPOT
   rm -rf $UNPETURBPOT/*
   cp -r $OUTDIR/* $UNPETURBPOT 
fi
########### Calculating perturbed system ###########

# Settig  directory for pertirbed calculation

PERTURBED_DIRECT=PERTURBED
 if [ ! -d  $PERTURBED_DIRECT ]; then
   mkdir -p $PERTURBED_DIRECT
 fi

###Inputs for reposne Matrix
# pos.in file
sed '1d' $COORDS_FE > tempxyz
cat tempxyz | awk '{print $2, $3, $4}' > tempxyz2
sed '1d' $LATTICEVECT > templatvect
cat templatvect tempxyz2 > ./${PERTURBED_DIRECT}/pos.in
rm tempxyz tempxyz2 templatvect

# r.x input
cat > ./${PERTURBED_DIRECT}/resp_mat.in << EOF
&input_mat
 ntyp = 2
 na(1) =4
 na(2) =6 
 nalfa = ${#ALPHA[@]}
 magn = .true.
 filepos = "pos.in"
 back = "no"
 filednda = "dnda.in"
 n1 = 5
 n2 = 5
 n3 = 5 
&end
EOF

### Description input
#&input_mat
# ntyp = 2    Nuber of type of atoms
# na(1) = 8 Number of type 1 
# na(2) = 8 Number of Type 2
# nalfa = 3
# magn = .true. True for magnetic systems, false otherwise
# filepos = "pos.in" File containing atomic positions 
# back = "no" To add a neutralizing backgroud see PRB 71-35105
# filednda = "dnda.in" File containing the names of dn*mat
# n1 = 5  #Number of unit cells used for extrapolating response 
# n2 = 5  # matrices to larged supercells 
# n3 = 5 
#&end
### End inputs response Matrix 

#Preparing input for perturbed system
ETHR=`grep ethr ${SYSTEM}.unperturbed_scf_Fe.out | tail -1 |awk '{print $3}'`
echo "ETHR" $ETHR
cp  ${SYSTEM}.unperturbed_scf_Fe.in ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe.in

# add these lines before conv_thr
#line=`cat -n ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe.in | grep conv_thr | head -1 | awk '{print $1+1}'`#Active if input we dont have forc and etot_thr
line=`cat -n ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe.in | grep 'conv_thr' | tail -2 | head -1 |  awk '{print $1+1}'` #When we have in input orc and etot_thr
sed ''${line}'i\    diago_thr_init = '$ETHR'' -i ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe.in
sed ''${line}'i\    startingpot = 'file',' -i ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe.in
sed ''${line}'i\    startingwfc = 'file',' -i ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe.in
# add perturbation lines after Hubbard_U line 
line=`cat -n ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe.in | grep Hubbard_U | tail -1 | awk '{print $1+1}'`
sed ''${line}'i\    Hubbard_alpha(1)= 0.00' -i ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe.in
# Settig outdir directory perturbed (cx1)
OUTDIR2=$TMPDIR/DFTU/${SYSTEM}_2
if [ ! -d  $OUTDIR2 ]; then
   mkdir -p $OUTDIR2
fi
# output directory setting
sed "/outdir=/{s#.*#    outdir= '$OUTDIR2/'#}" -i ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe.in
# perturbed run
for a in ${ALPHA[@]}; do
 rm -rf $OUTDIR2
 cp -r $UNPETURBPOT/* $OUTDIR2
 sed '/Hubbard_alpha(1)/{s/.*/    Hubbard_alpha(1)= '$a'/}' -i ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe.in
 cp ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe.in ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe_$a.in
 echo "${SYSTEM}: Calculating response to alpha=$a"
 #Checking if the .out files already exist
 if [ ! -f ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe_$a.out ]; then
    $MPI $PW_BIN < ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe_$a.in  $POOL_C >> ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe_$a.out 2> /dev/null
 fi # End checking .out exit
done
#Starting unperturbed calculatin for Oxygen
echo "Starting unperturbed calculation O..."

cat > ${SYSTEM}.unperturbed_scf_O.in << EOF
&control
    calculation='scf'
    verbosity = 'high'
    restart_mode='from_scratch'
    prefix='$SYSTEM'
    pseudo_dir = '$PSEUDO_DIR'
    outdir='$OUTDIR2/'
    wf_collect = .true.
    tprnfor = .true.
    tstress=.true.
    etot_conv_thr=1.0D-4
    forc_conv_thr=1.0D-3
    nstep=50000
 /
&system
    ibrav=0,
    celldm(1)=1.889726878,
    nat= $NPA, ntyp = $DIFFATOM, nspin = 2
    nbnd = 70,
    ecutwfc = $ECUTWFC, ecutrho = $ECUTRHO,
    occupations='smearing', smearing='mv', degauss=0.01,
    starting_magnetization(1)=0.6
    starting_magnetization(4)=-0.6
    lda_plus_u = .true.,
    U_projection_type = 'atomic',
    Hubbard_U(1)= 1.d-20
    Hubbard_U(2)= 1.d-20
    Hubbard_U(3)= 1.d-20
    Hubbard_U(4)= 1.d-20

 /
&electrons
    mixing_beta = 0.3 
    diagonalization='david'
    startingwfc='atomic+random',
    electron_maxstep=1000
    conv_thr =  1.0d-9
 /
 &ions
 /
 &cell
   press=0.0
   press_conv_thr=0.1
 /
ATOMIC_SPECIES
EOF
cat $SPECLISTO $COORDS_O $LATTICEVECT >> ${SYSTEM}.unperturbed_scf_O.in
cat >> ${SYSTEM}.unperturbed_scf_O.in << EOF
K_POINTS (automatic)
EOF
cat $KPOINTS >> ${SYSTEM}.unperturbed_scf_O.in

  # Checking if the .out files already exist
  if [ ! -f ${SYSTEM}.unperturbed_scf_O.out ]; then
     $MPI $PW_BIN < ${SYSTEM}.unperturbed_scf_O.in  $POOL_C >> ${SYSTEM}.unperturbed_scf_O.out 2> /dev/null
  fi # End checking .out exit

  echo "Ending unperturbed calculation..."

####Copying Potential from OUTDIR
UNPETURBPOTO=UPOTENTIALO
if [ ! -d  $UNPETURBPOTO ]; then
   mkdir -p $UNPETURBPOTO
   rm -rf $UNPETURBPOTO/*
   cp -r $OUTDIR2/* $UNPETURBPOTO
fi
########### Calculating perturbed system ###########

# Settig  directory for pertirbed calculationi Oxyge
#Preparing input for perturbed system
ETHR=`grep ethr ${SYSTEM}.unperturbed_scf_O.out | tail -1 |awk '{print $3}'`
echo "ETHR" $ETHR
cp  ${SYSTEM}.unperturbed_scf_O.in ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O.in
# Remove one line in input
sed '/startingwfc/d' -i ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O.in
# add these lines before conv_thr
#line=`cat -n ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O.in | grep conv_thr | head -1 | awk '{print $1+1}'`#Active if input we dont have forc and etot_thr
line=`cat -n ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O.in | grep 'conv_thr' | tail -2 | head -1 |  awk '{print $1+1}'` #When we have in input orc and etot_thr
sed ''${line}'i\    diago_thr_init = '$ETHR'' -i ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O.in
sed ''${line}'i\    startingpot = 'file',' -i ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O.in
sed ''${line}'i\    startingwfc = 'file',' -i ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O.in
# add perturbation lines after Hubbard_U line 
line=`cat -n ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O.in | grep Hubbard_U | tail -1 | awk '{print $1+1}'`
sed ''${line}'i\    Hubbard_alpha(2)= 0.00' -i ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O.in
# Settig outdir directory perturbed (cx1)
OUTDIR3=$TMPDIR/DFTU/${SYSTEM}_3
if [ ! -d  $OUTDIR3 ]; then
   mkdir -p $OUTDIR3
fi
# output directory setting
sed "/outdir=/{s#.*#    outdir= '$OUTDIR3/'#}" -i ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O.in
# perturbed run
for a in ${ALPHA[@]}; do
 rm -rf $OUTDIR3
 cp -r $UNPETURBPOTO/* $OUTDIR3
 sed '/Hubbard_alpha(2)/{s/.*/    Hubbard_alpha(2)= '$a'/}' -i ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O.in
 cp ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O.in ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O_$a.in
 echo "${SYSTEM}: Calculating response to alpha=$a"
 #Checking if the .out files already exist
 if [ ! -f ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O_$a.out ]; then
    $MPI $PW_BIN < ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O_$a.in  $POOL_C >> ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O_$a.out 2> /dev/null
 fi # End checking .out exit

# generate dnda files for reponsematrix 
#For Fe
 n0totFe=(`grep "Tr" ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe_$a.out | head -$NPA2 | tail -$NPA | awk '{print $NF}'`)
 ntotFe=(`grep "Tr"  ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_Fe_$a.out | tail -$NPA | awk '{print $NF}'`)
 np=1
 for (( idx=0; idx<$NPA; idx++)); do
   ni=$(( $idx + 1 ))
   echo $a ${ntotFe[idx]} >> ./${PERTURBED_DIRECT}/dn_${ni}_da_${np}.dat
   echo $a ${n0totFe[idx]} >> ./${PERTURBED_DIRECT}/dn0_${ni}_da_${np}.dat
   if [[ "$a" == "${ALPHA[0]}" ]]; then
     echo " dn_${ni}_da_${np}.dat  dn0_${ni}_da_${np}.dat" >>  ./${PERTURBED_DIRECT}/dnda.in
   fi
 done
 #For O
 n0totO=(`grep "Tr" ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O_$a.out | head -$NPA2 | tail -$NPA | awk '{print $NF}'`)
 ntotO=(`grep "Tr"  ./${PERTURBED_DIRECT}/${SYSTEM}.perturbed_scf_O_$a.out | tail -$NPA | awk '{print $NF}'`)
 np=17
 for (( idx=0; idx<$NPA; idx++)); do
   ni=$(( $idx + 1 ))
   echo $a ${ntotO[idx]} >> ./${PERTURBED_DIRECT}/dn_${ni}_da_${np}.dat
   echo $a ${n0totO[idx]} >> ./${PERTURBED_DIRECT}/dn0_${ni}_da_${np}.dat
   if [[ "$a" == "${ALPHA[0]}" ]]; then
     echo " dn_${ni}_da_${np}.dat  dn0_${ni}_da_${np}.dat" >>  ./${PERTURBED_DIRECT}/dnda.in
   fi
 done

done #End for looping the perturbaions

echo "Ending perturbed calculations ..."

if [[ -z $LRPU ]]; then    # if LRPU option present skip this part
 # calculate U
 sed '/filednda/{s/.*/ filednda = "dnda.in"/}' -i ./${PERTURBED_DIRECT}/resp_mat.in
 sed '/filepos/{s/.*/ filepos = "pos.in"/}' -i ./${PERTURBED_DIRECT}/resp_mat.in
 for ni in 1 2 3 4 5;do #Extrapolation to superlattices 1x1x1 2x2x2 3x3x3 4x4x4 5x5x5
   cp ./${PERTURBED_DIRECT}/resp_mat.in ./${PERTURBED_DIRECT}/resp_mat_$ni.in
   sed '/n1/{s/.*/ n1 = '$ni'/}' -i ./${PERTURBED_DIRECT}/resp_mat_$ni.in
   sed '/n2/{s/.*/ n2 = '$ni'/}' -i ./${PERTURBED_DIRECT}/resp_mat_$ni.in
   sed '/n3/{s/.*/ n3 = '$ni'/}' -i ./${PERTURBED_DIRECT}/resp_mat_$ni.in
   cd ${PERTURBED_DIRECT}
   if [ ! -f Umat_$ni.out]; then
      r.x < resp_mat_$ni.in
      nat=`grep 'number of atoms in the supercell'  Umat.out |tail -1|awk '{print $7}'`
      u1=`grep U1 Umat.out | awk '{print $NF}' | head -1`
      u2=`grep U1 Umat.out | awk '{print $NF}' | tail -1`
      echo "$nat $u1 $u2" >> Usc1_pbe.out
      mv Umat.out Umat_$ni.out
      cd ..
   else
      nat=`grep 'number of atoms in the supercell'  Umat_$ni.out |tail -1|awk '{print $7}'`
      u1=`grep U1 Umat_$ni.out | awk '{print $NF}' | head -1`
      u2=`grep U1 Umat_$ni.out | awk '{print $NF}' | tail -1`
      echo "$nat $u1 $u2" >> Usc1_pbe.out
      cd ..
   fi  
 done
fi # End oif LRPU

############ Delete resutls_ folders ######################

if [[ -n $RMV ]]; then # if RMV option is present it removes the temporary files
   echo "Removing .save from" $TMPDIR
   rm -r $TMPDIR #archer
fi
