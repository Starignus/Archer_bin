#!/bin/sh
#
# Runing for relaxing a structure and PDOS
####################################################################

#### Specify these variables based on your environment
PSEUDO_DIR='/work/e05/e05/ablar/QE_PSEUDO'


#### System variables and files
SYSTEM=`pwd`
SYSTEM=${SYSTEM##*/}
SPECLIST=$SYSTEM'.spec' #File with the atomic specie, atomic mass, and pseudopotential
COORDS=$SYSTEM'.xyz'    #Fliel with the atomic coordinates in the crystaline system
KPOINTS=$SYSTEM'.kpoints'  #File with the kpoint mesh
KPOINTSDOS=$SYSTEM'.kpoints_pdos' #File with the kpoint mesh for PDOS
KPOINTSBAND=$SYSTEM'.kpoints_bands' #File with the high symmetry kpoint 
LATTICEVECT=$SYSTEM'.lattice_invect' # Used insteed ibrav
VERSION=1.0

###########Options###################
while getopts x:a:m:p:k:c:g:d:l:s:q:o:u:v:w:f:e:y:t:r:z:b: opt
do
   case "$opt" in
      x) BIN=$OPTARG;;
      a) NPOOL=$OPTARG;;      
      m) MPI=$OPTARG;;
      p) PO="yes";;
      k) RSTR="yes";;
      c) CAL=$OPTARG;;
      g) FUNC=$OPTARG;;
      d) PDOS="yes";;
      l) BIN_PROJ=$OPTARG;;
      s) BIN_DOS=$OPTARG;;
      q) OCCUP=$OPTARG;;
      o) BIN_BAND=$OPTARG;;
      u) BIN_BANDPLOT=$OPTARG;;
      v) echo "resend version " $VERSION; exit 1;;
      w) ECUTWFC=$OPTARG;;
      f) ECUTRHO=$OPTARG;;
      e) RMV="yes";;
      y) CPWF="yes";;
      t) WALLT=$OPTARG;;
      r) REST=$OPTARG;;
      z) PLOTEV="yes";;
      b) SUMDOS="yes";;
      \?) echo """   resend usage:
    -x  path to the pw.x binary (if not in the \$PATH)
    -a  Number of npools
    -m  MPI driver with options (i.e. "mpirun -np 8")
    -p  It assumes the calculation (vc-relax, scf, relax) is already done previously  (i.e. -p "yes")
    -k  If it is present, it checks for the output file to get the relaxed or the last structure obtained  
        (i.e. -k "yes"). If we don't want to run relxation again also put -p "yes" option
    -c  Type of calculation:scf, vc-relax ,relax
    -g  Functional: PBE,PBEsol, etc.
    -d  Calculate PDOS  (i.e. -r "yes")
    -l  Path to the projwfc.x binary (if not in the \$PATH)
    -s  Path to the dos.x binary (if not in the \$PATH)
    -q  Kind of occupation in PDOS:smearing (metals),tetrahedra (for DOS), fixed (for insulators) and from_input (occupation are read from input:go to manual)
    -o  Path to the bands.x binary (if not in the \$PATH)
    -u  Path to the plotband.x  binary (if not in the \$PATH)
    -v  echo "resend version " $VERSION; exit 1;;
    -w  Cutoff energy obtained before in Ry.
    -f  Cutoffrho.
    -r  Delete the folderes "results_" (i.e. -r "yes")
    -y  Copy the "$SYSTEM".save from the temporary folder to "$SYSTEM"
    -t  Wall time in seconds (Default 24 hours = 86400 sec)
    -r  Restarting mode: restart, from_scratch (default)
    -z  Option to plot E Vs. volume curve (i.e. -z "yes")
    -b  Option to use the sum_states.py for DOS
    -h  print this help 
   Example:time Fe_FM_RUN -x /prog/x86_64/espresso/espresso-5.0.2/bin/pw.x -m "mpirun -np 2" -f 200-w 30 -c scf -g PBE
           time Fe_FM_RUN -x /prog/x86_64/espresso/espresso-5.0.2/bin/pw.x -m "mpirun -np 2" -f 200 -w 30 -p "yes" -r "yes" -c scf -g PBE
           time Fe_FM_RUN -x /prog/x86_64/espresso/espresso-5.0.2/bin/pw.x  -m "mpirun -np 2" -w 50 -f 200 -c scf -g PBE (when running locally)
           In my desktop:
           time Fe_FM_RUN -x "\$espresso" -m "mpirun -np 4" -w 50 -f 200 -c vc-relax -g PBE -d yes -l "\$projwfc" -s "\$dos" -q smearing
           time Fe_FM_RUN.sh -x \$espresso -m "mpirun -np 4" -w 50 -f 500 -c "vc-relax" -g PBE -z "yes" -d "yes" -l \$projwfc -s \$dos -q "smearing" -p "yes" -o \$bands -u \$plotbands
           time Fe_FM_RUN.sh -x \$espresso -m "mpirun -np 4" -w 50 -f 500 -c vc-relax -g PBE -p "yes" -d yes -l "\$projwfc" -s "\$dos" -q smearing -o \$bands -u \$plotbands (for runing just PDOS part)
"""
     exit 1  ;;
   esac
done  

######################## End Options #######################

echo "    S_FM_primitive_RUN.sh.v." $VERSION
echo "    -------------------" 

##################### System setings #######################

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

echo "pw.x  with n pools executed as: " $MPI $PW_BIN "< input" $POOL_C  "> output"


# Checking wall time

if [ "$WALLT" ]; then
    SECOND=$WALLT
  elif [ -z  "$WALLT" ]; then
    SECOND="86400"   
fi

# Restart mode 

if [ "$REST" ]; then
    RESTMOD=$REST
  elif [ -z "$REST" ]; then
    RESTMOD="from_scratch"
fi   

# Cheching for .spect, .xyz. .kpoints .kpoint_pdos
if [[ ! -e $SPECLIST ]] || [[ ! -e $COORDS ]] || [[ ! -e $KPOINTS ]] || [[ ! -e $KPOINTSDOS ]] || [[ ! -e $KPOINTSBAND ]] || [[ ! -e $LATTICEVECT ]]; then
   echo " # One or more files for the input are missing ( .spec, .xyz. .kpoints .kpoints_pdos .kpoints_bands)"
   exit
fi

if [[ -z $PO ]]; then    # if PO option present skip this part
echo "Starting" $CAL "..."

# Settig scratch directory for relaxing
# if [ ! -d  $TMPDIR/results_$SYSTEM'_fm_'$CAL ]; then
#   mkdir   $TMPDIR/results_$SYSTEM'_fm_'$CAL
# fi

#SCRATCHDIR=$TMPDIR/results_$SYSTEM'_fm_'$CAL
SCRATCHDIR=./


echo "  Temporary directory= " $SCRATCHDIR
echo "  calculating FM=" " ... "
if [[ -e $SYSTEM'_fm_2'.in ]]; then  #Checking if there is already a relaxation
    echo "Perfoming again the relaxation:"  $SYSTEM'_fm_2'.in 
    $MPI $PW_BIN < $SYSTEM'_fm_2'.in $POOL_C >> $SYSTEM'_fm_2'.out 2> /dev/null
    echo "copy" $SYSTEM'_fm_2'.in "to another file to keep the record of structures"
    cp $SYSTEM'_fm_2'.in $SYSTEM'_fm_2_1'.in
    echo "Copying" $SYSTEM'_fm_'.out "to" $SYSTEM'_fm_1'.out 
    cp $SYSTEM'_fm_'.out $SYSTEM'_fm_1'.out
    echo "Moving" $SYSTEM'_fm_2'.out "to" $SYSTEM'_fm_'.out
    cp $SYSTEM'_fm_2'.out $SYSTEM'_fm_'.out
    echo "copying" $SYSTEM'_fm_2'.in "in" $SYSTEM'_fm_'.in
    cp  $SYSTEM'_fm_2'.in $SYSTEM'_fm_'.in
   else

###################  Type of CAL=relax,scf, vc-relax  ################
cat >  $SYSTEM'_fm_'.in << EOF
&control
    calculation='$CAL'
    restart_mode='$RESTMOD',    
    pseudo_dir ='$PSEUDO_DIR',
    outdir='$SCRATCHDIR'
    prefix='$SYSTEM',
    tprnfor = .true.
    tstress=.true.
    nstep=50000
    etot_conv_thr=1.0D-5
    forc_conv_thr=1.0D-4
    disk_io='default'
    max_seconds=$SECOND 

 /
 &system
    ibrav=0, 
    celldm(1)=1,
    nat=128 , ntyp= 1, ecutwfc = $ECUTWFC, ecutrho = $ECUTRHO,
    occupations='smearing', smearing='mv', degauss=0.01,
    input_DFT='$FUNC'
    london=.true.
    london_rcut = 200
    london_s6 = 0.75

 /
 &electrons
    diagonalization='david'
    electron_maxstep=1000
    conv_thr = 1.0e-8
    mixing_beta = 0.2
 /
 &ions
 /
 &cell
   press=0.0
   press_conv_thr=0.1
 /
ATOMIC_SPECIES
EOF
cat $SPECLIST $COORDS $LATTICEVECT >>$SYSTEM'_fm_'.in
cat >> $SYSTEM'_fm_'.in << EOF
K_POINTS (automatic)
EOF
cat $KPOINTS >> $SYSTEM'_fm_'.in

# Checking if the .out files already exist
if [ ! -f $SYSTEM'_fm_'.out ]; then
   $MPI $PW_BIN < $SYSTEM'_fm_'.in $POOL_C >> $SYSTEM'_fm_'.out 2> /dev/null
fi # End checking .out exit

echo "Ending" $CAL "..."
fi # If for second relaxation
fi # end of PO section

if [ -n "$RSTR" ]; then # In case we want to do this without runing again the previous scf

echo "Cheking if there is a relaxed structure ..."
############ Checking if there is  a relaxed structure ########################

if [[ -e $SYSTEM'_fm_'.out ]] &&  [[ "$CAL" == *vc-relax* ]]; then  # Here I will check if I got a relaxed structure
   CORRD=`grep 'Begin final coordinates'  $SYSTEM'_fm_'.out` 
   if [ -z "$CORRD" ]; then # If we don't have relaxed structure but we want to get the last coordinates
       RESTMOD="restart"
       echo "The final structure was not found, check your ouput. "
       echo "Check the new input with the last structure generated "   
       cp $SYSTEM'.xyz' $SYSTEM'_1.xyz' 
       LASTSTEP=`grep !  $SYSTEM'_fm_'.out | awk 'END{print}'`
       grep -A 150 "$LASTSTEP"  $SYSTEM'_fm_'.out > TEMP'_fm_' # We have to change the number 150 if we have more than 30 atoms
       ########### Getting lattice parameters ##################################
       grep -A 3 "CELL"  TEMP'_fm_' > $SYSTEM'_fm_latvec' # Print just the lattice vectors
       grep "new unit-cell volume" TEMP'_fm_' > $SYSTEM'_fm_volume'
       alat=`grep 'alat' $SYSTEM'_fm_latvec' | awk 'NR=1{ print substr($3,1,10)}'` #to get initial lattice parameters  
       Normvec_angle.py $SYSTEM'_fm_latvec' $alat > $SYSTEM'_fm_lattice'
       nalat=`awk 'NR==2{print $1}' $SYSTEM'_fm_lattice' ` # new lattice parameters in Bhors to put in celldm()
       nblat=`awk 'NR==2{print $2}' $SYSTEM'_fm_lattice' `
       nclat=`awk 'NR==2{print $3}' $SYSTEM'_fm_lattice' `  
       ####for the hexagonal case ###
       #nac=$(echo "$nclat/$nalat" | bc -l)    
       ###### Here this section won't be used because this is a convetional cell ##############
       #volbohrp=`grep 'volume' $SYSTEM'_fm_volume_p' | awk '{print $5}'` # Primtive volume bhor
       #volangp=`grep 'volume' $SYSTEM'_fm_volume_p' | awk '{print $8}'` # Primitive volume angs
       #volbohrc=$(echo "2*$volbohrp" | bc -l ) # Conventional volume bohr 
       #volangc=$(echo "2*$volangp" | bc -l) #Conventional volume angs
       #echo "Conventional cell (Bohr):" $volbohrc "Conventional cell (angs):" $volangc > $SYSTEM'_fm_volume_c'
       #nalatbhorc=$(echo "e(l($volbohrc)/3)" | bc -l ) # new lattice parameters in bhor of conventional  V_prim= (a_con)^3/2
       #nalatagnsc=$(echo "$nalatbhorc*0.529177" | bc -l)  # new lattice parameters in angs
       #echo "Conventional a (Bohr):" $nalatbhorc "Conventional a (angs):" $nalatagnsc >  $SYSTEM'_fm_alattice_c'       
       ############ Getting atomic coordinates #################################
       nat=`grep nat $SYSTEM'_fm_'.in | awk '{ print substr($1,5,2)}'` #If the number of atoms is more than to 2 ciphers or digits we most change the "2" to "3" 
       natline=$(echo "$nat+5" | bc)
       grep -A $natline "CELL" TEMP'_fm_' > COORD'_fm_'
       grep -A $nat 'ATOMIC_POSITIONS' COORD'_fm_'  > $SYSTEM'_fm_atomic_possitions'  
       cp $SYSTEM'_fm_atomic_possitions'  $SYSTEM'.xyz' # Coping the new atomic possitions
       rm  COORD'_fm_' TEMP'_fm_'
       ##########  Building new input #########################################
cat >  $SYSTEM'_fm_2'.in << EOF
&control
    calculation='$CAL'
    restart_mode='$RESTMOD',    
    pseudo_dir = '$PSEUDO_DIR',
    outdir='$SCRATCHDIR'
    prefix='$SYSTEM',
    tprnfor = .true.
    tstress=.true.
    nstep=50000
    etot_conv_thr=1.0D-5
    forc_conv_thr=1.0D-4
    disk_io='default'
    max_seconds=$SECOND 

 /
 &system
    ibrav=0, 
    celldm(1)=1,
    nat=128 , ntyp= 1, ecutwfc = $ECUTWFC, ecutrho = $ECUTRHO,    
    occupations='smearing', smearing='mv', degauss=0.01,
    input_DFT='$FUNC'
    london=.true.
    london_rcut = 200
    london_s6 = 0.75

 /
 &electrons
    diagonalization='david'
    electron_maxstep=1000
    conv_thr = 1.0e-8
    mixing_beta = 0.5
 /
 &ions
 /
 &cell
   press=0.0
   press_conv_thr=0.1
 /
ATOMIC_SPECIES
EOF
cat $SPECLIST $COORDS $SYSTEM'_fm_latvec' >>$SYSTEM'_fm_2'.in
cat >> $SYSTEM'_fm_2'.in << EOF
K_POINTS (automatic)
EOF
cat $KPOINTS >> $SYSTEM'_fm_2'.in
       exit
    else #if we obtain a final structure
       RESTMOD="from_scratch"
       echo "The final structure was found, check if the last scf cycle was performed. "  
       cp $SYSTEM'.xyz' $SYSTEM'_1.xyz'
       ########### Getting lattice parameters ##################################
       grep -A 1 'Begin final coordinates' $SYSTEM'_fm_'.out  > $SYSTEM'_fm_final_volume' #Print the final volume
       grep -A 6 'Begin final coordinates' $SYSTEM'_fm_'.out  > $SYSTEM'_fm_CELL'  
       grep -A 3 'CELL' $SYSTEM'_fm_CELL' > $SYSTEM'_fm_latvec_final'  # Print just the lattice vectors
       alat=`grep 'alat' $SYSTEM'_fm_latvec_final' | awk 'NR=1{ print substr($3,1,10)}'` #to get initial lattice parameters given in input as well
       Normvec_angle.py $SYSTEM'_fm_latvec_final' $alat > $SYSTEM'_fm_lattice_final'
       nalat=`awk 'NR==2{print $1}' $SYSTEM'_fm_lattice_final' ` # new lattice parameters in Bhors of the primitive cell
       nblat=`awk 'NR==2{print $2}' $SYSTEM'_fm_lattice_final' `
       nclat=`awk 'NR==2{print $3}' $SYSTEM'_fm_lattice_final' `
       ########for the hexagonal case ####
       #nac=$(echo "$nclat/$nalat" | bc -l)
       ###### Here this section won't be used because this is a convetional cell ##############
       #volbohrp=`grep 'volume' $SYSTEM'_fm_final_volume' | awk '{print $5}'` # Primtive volume bhor
       #volangp=`grep 'volume' $SYSTEM'_fm_final_volume' | awk '{print $8}'` # Primitive volume angs
       #volbohrc=$(echo "2*$volbohrp" | bc -l ) # Conventional volume bohr 
       #volangc=$(echo "2*$volangp" | bc -l) #Conventional volume angs
       #echo "Conventional cell (Bohr):" $volbohrc "Conventional cell (angs):" $volangc > $SYSTEM'_fm_final_volume'
       #nalatbhorc=$(echo "e(l($volbohrc)/3)" | bc -l ) # new lattice parameters in bhor of conventional  V_prim= (a_con)^3/2
       #nalatagnsc=$(echo "$nalatbhorc*0.529177" | bc -l )  # new lattice parameters in angs
       #echo "Conventional a (Bohr):" $nalatbhorc "Conventional a (angs):" $nalatagnsc >  $SYSTEM'_fm_final_alattice'
       rm $SYSTEM'_fm_CELL'
       ############ Getting atomic coordinates #################################
       nat=`grep nat $SYSTEM'_fm_'.in | awk '{ print substr($1,5,2)}'` #If the number of atoms is more of 2 ciphers or digits we most change the "2" to "3" 
       natline=$(echo "$nat+8" | bc)
       grep -A $natline 'Begin final coordinates' $SYSTEM'_fm_'.out > COORD'_fm_'
       grep -A $nat 'ATOMIC_POSITIONS' COORD'_fm_'  > $SYSTEM'_fm_final_atomic_possitions'  
       cp $SYSTEM'_fm_final_atomic_possitions'  $SYSTEM'.xyz'
       rm  COORD'_fm_'
cat >  $SYSTEM'_fm_2'.in << EOF
&control
    calculation='$CAL'
    restart_mode='$RESTMOD',    
    pseudo_dir = '$PSEUDO_DIR',
    outdir='$SCRATCHDIR'
    prefix='$SYSTEM',
    tprnfor = .true.
    tstress=.true.
    nstep=50000
    etot_conv_thr=1.0D-5
    forc_conv_thr=1.0D-4
    disk_io='default'
    max_seconds=$SECOND 

 /
 &system
    ibrav=0, 
    celldm(1)=1,
    nat=128 , ntyp= 1, ecutwfc = $ECUTWFC, ecutrho = $ECUTRHO,        
    occupations='smearing', smearing='mv', degauss=0.01,
    input_DFT='$FUNC'
    london=.true.
    london_rcut = 200
    london_s6 = 0.75
 /
 &electrons
    diagonalization='david'
    electron_maxstep=1000
    conv_thr = 1.0e-8
    mixing_beta = 0.5
 /
 &ions
 /
 &cell
   press=0.0
   press_conv_thr=0.1
 /
ATOMIC_SPECIES
EOF
cat $SPECLIST $COORDS  $SYSTEM'_fm_latvec_final'   >>$SYSTEM'_fm_2'.in
cat >> $SYSTEM'_fm_2'.in << EOF
K_POINTS (automatic)
EOF
echo "Copying original kpoints into a new file structure ..." $SYSTEM'_fm_2'.in
cat $KPOINTS >> $SYSTEM'_fm_2'.in
       #################### List of energies and volumes#############
       energies=$(grep -a "!    total energy              =" $SYSTEM'_fm_'.out | awk '{print $5}' )  
       volumescbohr=$(grep -a "new unit-cell volume" $SYSTEM'_fm_'.out | awk '{print $5*2}' ) # volumes conventional bohr
       volumescangsr=$(grep -a "new unit-cell volume" $SYSTEM'_fm_'.out | awk '{print $8*2}' ) # volumes primitive angrom
       ######### arrays ####
       declare -a array_volumebohr
       declare -a array_volumangs
       declare -a array_ener
       array_volumebohr=(${volumescbohr// / })
       array_volumeangrs=(${volumescangsr// / }) 
       array_ener=(${energies// / }) 
       rm -f energies.$SYSTEM'_fm_angs' energies.$SYSTEM'_fm_bohr'
       for i in `seq 0 ${#array_ener[@]}`; do
        echo ${array_ener[$i]} ${array_volumeangrs[$i]}  >>  energies.$SYSTEM'_fm_angs'
        echo ${array_ener[$i]} ${array_volumebohr[$i]}  >>  energies.$SYSTEM'_fm_bohr'  
       done
       ############## Plot Energy Vs. volume #############
       GP_COMMAND=`which gnuplot 2>/dev/null`  
       if [[ -n "$PLOTEV" ]] && [[ -e energies.$SYSTEM'_fm_bohr' ]] && [[ -e energies.$SYSTEM'_fm_angs' ]]; then
          echo "Plotting E Vs. V ..."
          if [ "$GP_COMMAND" != "" ]; then
             $GP_COMMAND -persist -e "set xlabel 'Volume (Ang^3)'; set ylabel 'E (Ry)'; p 'energies.${SYSTEM}_fm_angs' u 1:2 w linespoints"
             $GP_COMMAND -e "set terminal png; set output 'fig.energies.${SYSTEM}_fm_angs.png'; set xlabel 'Volume (Ang^3)'; set ylabel 'E (Ry)'; p 'energies.${SYSTEM}_fm_angs' u 1:2 w linespoints"
             $GP_COMMAND -persist -e "set xlabel 'Volume (u.a.^3)'; set ylabel 'E (a.u.)'; p 'energies.${SYSTEM}_fm_bohr' u 1:2 w linespoints"
             $GP_COMMAND -e "set terminal png; set output 'fig.energies.${SYSTEM}_fm_bohr.png'; set xlabel 'Volume (a.u.^3)'; set ylabel 'E (Ry)'; p 'energies.${SYSTEM}_fm_bohr' u 1:2 w linespoints"
          else
            d        $ECHO "No gnuplot in PATH. Results not plotted."
          fi
        else
          echo ' ' 
          echo ' # ERROR: There is no energy file (called energies. $SYSTEM_fm_angs (bohr). Unable to plot.' 
       fi # End Plotting E Vs. V 
  fi # End for obtaining new input with final relaxed strucure or last structure
fi # end if we get a vc-relax structure or not

####### CAL relax #####
if [[ -e $SYSTEM'_fm_'.out ]] && [[  "$CAL" == "relax" ]]; then  # Here I will check if I got a relaxed structure
   CORRD=`grep 'Begin final coordinates'  $SYSTEM'_fm_'.out`
   if [ -z "$CORRD" ]; then # If we don't have relaxed structure but we want to get the last coordinates
       RESTMOD="restart"
       echo "The final structure was not found, check your ouput. "
       echo "Check the new input with the last structure generated"
       cp $SYSTEM'.xyz' $SYSTEM'_1.xyz'   
       LASTSTEP=`grep !  $SYSTEM'_fm_'.out | awk 'END{print}'`
       grep -A 150 "$LASTSTEP"  $SYSTEM'_fm_'.out > TEMP'_fm_' # We have to change the number 150 if we have more than 30 atoms
       ############ Getting atomic coordinates #################################
       nat=`grep nat $SYSTEM'_fm_'.in | awk '{ print substr($1,5,2)}'` #If the number of atoms is more of 2 ciphers or digits we most change the "2" to "3" 
       natline=$(echo "$nat+2" | bc)
       grep -A $natline "CELL" TEMP'_fm_' > COORD'_fm_'
       grep -A $nat 'ATOMIC_POSITIONS' COORD'_fm_'  > $SYSTEM'_fm_final_atomic_possitions'
       cp $SYSTEM'_fm_final_atomic_possitions'  $SYSTEM'.xyz'
       rm  COORD'_fm_' TEMP'_fm_'
       cp $LATTICEVECT $SYSTEM'_fm_latvec'
       ##########  Building new input #########################################
cat >  $SYSTEM'_fm_2'.in << EOF
&control
    calculation='$CAL'
    restart_mode='$RESTMOD',    
    pseudo_dir = '$PSEUDO_DIR',
    outdir='$SCRATCHDIR'
    prefix='$SYSTEM',
    tprnfor = .true.
    tstress=.true.
    nstep=50000
    etot_conv_thr=1.0D-5
    forc_conv_thr=1.0D-4
    disk_io='default'
    max_seconds=$SECOND 

 /
 &system
    ibrav=0, 
    celldm(1)=1,
    nat=128 , ntyp= 1, ecutwfc = $ECUTWFC, ecutrho = $ECUTRHO,    
    occupations='smearing', smearing='mv', degauss=0.01,
    input_DFT='$FUNC'
    london=.true.
    london_rcut = 200
    london_s6 = 0.75

 /
 &electrons
    diagonalization='david'
    electron_maxstep=1000
    conv_thr = 1.0e-8
    mixing_beta = 0.5
 /
 &ions
 /
 &cell
   press=0.0
   press_conv_thr=0.1
 /
ATOMIC_SPECIES
EOF
cat $SPECLIST $COORDS $LATTICEVECT >>$SYSTEM'_fm_2'.in
cat >> $SYSTEM'_fm_2'.in << EOF
K_POINTS (automatic)
EOF
cat $KPOINTS >> $SYSTEM'_fm_2'.in
       exit
    else   #if we obtain a final structure
       RESTMOD="from_scratch"
       echo "The final structure (relax) was found, check if the last scf cycle was performed. "  
       cp $SYSTEM'.xyz' $SYSTEM'_1.xyz'
        ############ Getting atomic coordinates #################################
       nat=`grep nat $SYSTEM'_fm_'.in | awk '{ print substr($1,5,2)}'` #If the number of atoms is more of 2 ciphers or digits we most change the "2" to "3" 
       natline=$(echo "$nat+2" | bc)
       grep -A $natline 'Begin final coordinates' $SYSTEM'_fm_'.out > COORD'_fm_'
       grep -A $nat 'ATOMIC_POSITIONS' COORD'_fm_'  > $SYSTEM'_fm_final_atomic_possitions'
       cp $SYSTEM'_fm_final_atomic_possitions'  $SYSTEM'.xyz'
       rm  COORD'_fm_'
       cp $LATTICEVECT $SYSTEM'_fm_latvec'
cat >  $SYSTEM'_fm_2'.in << EOF
&control
    calculation='$CAL'
    restart_mode='$RESTMOD',    
    pseudo_dir = '$PSEUDO_DIR',
    outdir='$SCRATCHDIR'
    prefix='$SYSTEM',
    tprnfor = .true.
    tstress=.true.
    nstep=50000
    etot_conv_thr=1.0D-5
    forc_conv_thr=1.0D-4
    disk_io='default'
    max_seconds=$SECOND 

 /
 &system
    ibrav=0, 
    celldm(1)=1,
    nat=128 , ntyp= 1, ecutwfc = $ECUTWFC, ecutrho = $ECUTRHO,    
    occupations='smearing', smearing='mv', degauss=0.01,
    input_DFT='$FUNC'
    london=.true.
    london_rcut = 200
    london_s6 = 0.75

 /
 &electrons
    diagonalization='david'
    electron_maxstep=1000
    conv_thr = 1.0e-6
    mixing_beta = 0.5
 /
 &ions
 /
 &cell
   press=0.0
   press_conv_thr=0.1
 /
ATOMIC_SPECIES
EOF
cat $SPECLIST $COORDS $LATTICEVECT >>$SYSTEM'_fm_2'.in
cat >> $SYSTEM'_fm_2'.in << EOF
K_POINTS (automatic)
EOF
cat $KPOINTS >> $SYSTEM'_fm_2'.in
   fi # End for obtaining new input with final relaxed strucure or last structure
fi # end if we have relax
echo "Ending section for checking if there is a relaxed structure ..."
fi # End of $RSTR
   

######################## Calculation of PDOS ###############################
if [[ -n $PDOS ]]; then   # if PDOS option present it calculates  the PDOS
   echo "Enterign to PDOS"
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
 
   # Check for dos.x binary 
   DOSWFC_BIN_=`which dos.x 2>/dev/null`
   if [ "$BIN_DOS" ]; then  # Checking the path of dos.x
      DOSWFC_BIN=$BIN_DOS
   fi

   if [ -z "$DOSWFC_BIN" ]; then
       echo " # ERROR: Provide a path to a dos.x binary."
       exit
   fi

   echo "dos.x executed as: " $MPI $DOSWFC_BIN 
   
   # Check for bands.x binary
   BANDS_BIN=`which bands.x 2>/dev/null`
   if [ "$BIN_BAND" ]; then  # Checking the path of projwfc.x
      BANDS_BIN=$BIN_BAND
   fi

   if [ -z "$BANDS_BIN" ]; then
       echo " # ERROR: Provide a path to a bands.x binary."
       exit
   fi

   echo "bands.x executed as: " $MPI $BANDS_BIN
   
   # Check for plotband.x binary 
   BANDPLOT_BIN=`which plotband.x 2>/dev/null`
   if [ "$BIN_BANDPLOT" ]; then  # Checking the path of dos.x
      BANDPLOT_BIN=$BIN_BANDPLOT
   fi

   if [ -z "$BANDPLOT_BIN" ]; then
       echo " # ERROR: Provide a path to a plotband.x binary."
       exit
   fi

   echo "plotband.x executed as: " $MPI $DOSWFC_BIN 

   # Creating PDOS directory
   if [ ! -d $SYSTEM'_fm_PDOSBAND' ]; then
       mkdir  $SYSTEM'_fm_PDOSBAND'
       echo "Making" $SYSTEM'_fm_PDOSBAND'
      # mkdir  $SYSTEM'_fm_PBANDS'
       #echo "Making" $SYSTEM'_fm_PDOSBAND'
   fi
      
   # Occupations for PDOS
   if [[ "$OCCUP" == *tetrahedra* ]]; then
      OCC="tetrahedra"
      echo "type occupation for DOS:" $OCC
   elif [[ "$OCCUP" == *smearing* ]]; then
      OCC="smearing"
      OPTIONS="smearing='mv', degauss=0.01,"
      echo "type occupation for DOS:" $OCC
   fi
    
# Settig scratch directory for scf
#   if [ ! -d results_$SYSTEM'_fm_scf' ]; then
#      mkdir  results_$SYSTEM'_fm_scf'
#      echo "Making" results_$SYSTEM'_fm_scf'
#   fi
#SCRATCHDIRDOS=./results_$SYSTEM'_fm_scf'
SCRATCHDIRDOS=./

echo "  Temporary directory=" $SCRATCHDIRDOS 

# alat value
   if [[ -e $SYSTEM'_fm_lattice_final' ]] &&  [[ "$CAL" == "vc-relax" ]]; then
       nalat=`awk 'NR==2{print $1}' $SYSTEM'_fm_lattice_final' ` # new lattice parameters in Bhors of the primitive cell
       nblat=`awk 'NR==2{print $2}' $SYSTEM'_fm_lattice_final' `
       nclat=`awk 'NR==2{print $3}' $SYSTEM'_fm_lattice_final' `
       #nac=$(echo "$nclat/$nalat" | bc -l)
       #volbohrp=`grep 'volume' $SYSTEM'_fm_final_volume_p' | awk '{print $5}'` # Primtive volume bhor
       #volbohrc=$(echo "2*$volbohrp" | bc -l ) # Conventional volume bhor 
       #nalatbhorc=$(echo "e(l($volbohrc)/3)" | bc -l ) # new lattice parameters in bhor of conventional  V_prim= (a_con)^3/2
     #elif [ "$CAL" == "relax" ]; then
       #nalatbhor=`grep "(alat)" $SYSTEM'_fm_'.out  | awk '{print $5}'` #to get initial lattice parameters 
   fi   
# self-consistent calculation
     echo "Setting up a scf with the relaxed structure; a_primitive="$nalat " a_conventional="$nalatbhorc
cat >  ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_scf'.in << EOF
&control
    calculation='scf'
    restart_mode='from_scratch',    
    pseudo_dir = '$PSEUDO_DIR',
    outdir='$SCRATCHDIRDOS'
    prefix='$SYSTEM',
    tprnfor = .true.
    tstress=.true.
    nstep=50000
    etot_conv_thr=1.0D-5
    forc_conv_thr=1.0D-4
    verbosity='high'
    disk_io='default'

 /
 &system
    ibrav=0, 
    celldm(1)=1,
    nat=128 , ntyp= 1, ecutwfc = $ECUTWFC, ecutrho = $ECUTRHO,        
    occupations='smearing', smearing='mv', degauss=0.01,
    nbnd = 88,
    input_DFT='$FUNC'
    london=.true.
    london_rcut = 200
    london_s6 = 0.75

 /
 &electrons
    diagonalization='david'
    electron_maxstep=1000
    conv_thr = 1.0e-8
    mixing_beta = 0.5
 /
 &ions
 /
 &cell
   press=0.0
   press_conv_thr=0.1
 /
ATOMIC_SPECIES
EOF
cat $SPECLIST $COORDS $SYSTEM'_fm_latvec_final'  >> ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_scf'.in
cat >> ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_scf'.in << EOF
K_POINTS (automatic) 
EOF
cat $KPOINTS >> ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_scf'.in
  if [ ! -f ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_scf'.out ]; then
    $MPI $PW_BIN < ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_scf'.in $POOL_C >> ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_scf'.out 2> /dev/null
        EF=`grep 'the Fermi energy' ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_scf'.out | awk '{print $5}'`
        echo "EFermi  from scf=" $EF
        grep 'the Fermi energy' ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_scf'.out | awk '{print $5}' > $SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_EF'
        cp ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_scf'.out ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM.scf
    else 
        EF=`grep 'the Fermi energy' ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_scf'.out | awk '{print $5}'`
        echo "EFermi  from scf=" $EF
        grep 'the Fermi energy' ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_scf'.out | awk '{print $5}' > $SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_EF'
        cp ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_scf'.out ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM.scf
  fi
       
################ BAND structure calculation along high-symmetry lines##################
    echo "Seeting structure for bands calculation ..."
cat >  $SYSTEM'_fm_bands'.in << EOF
&control
    calculation='bands'
    restart_mode='from_scratch',    
    pseudo_dir = '$PSEUDO_DIR',
    outdir='$SCRATCHDIRDOS'
    prefix='$SYSTEM',
    tprnfor = .true.
    tstress=.true.
    nstep=50000
    etot_conv_thr=1.0D-5
    forc_conv_thr=1.0D-4
    verbosity='high'
    disk_io='default'

 /
 &system
    ibrav=0, 
    celldm(1)=1,
    nat=128 , ntyp= 1, ecutwfc = $ECUTWFC, ecutrho = $ECUTRHO,        
    occupations='$OCC', $OPTIONS
    nbnd = 88,
    input_DFT='$FUNC'
    london=.true.
    london_rcut = 200
    london_s6 = 0.75

 /
 &electrons
    diagonalization='david'
    electron_maxstep=1000
    conv_thr = 1.0e-8
    mixing_beta = 0.5
 /
 &ions
 /
 &cell
   press=0.0
   press_conv_thr=0.1
 /
ATOMIC_SPECIES
EOF
cat $SPECLIST $COORDS $SYSTEM'_fm_latvec_final'  $KPOINTSBAND  >> $SYSTEM'_fm_bands'.in
  if [ ! -f $SYSTEM'_fm_bands'.out ]; then
    $MPI $PW_BIN < $SYSTEM'_fm_bands'.in $POOL_C >> $SYSTEM'_fm_bands'.out 2> /dev/null
    echo "Finished to run bands..."
    cp $SYSTEM'_fm_bands'.out ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM.pwscf
  fi
############## BANDS ############
#### New bands bands_inter_2.x  ####
if [ -e $SYSTEM.inp ]; then
   cp $SYSTEM.inp ./$SYSTEM'_fm_PDOSBAND'
   echo "Copying "  $SYSTEM.inp "to plot later the bands with bands_inter_2.x"  
  else
   echo "Remember to generate the file" $SYSTEM.inp " for ploting the bands"
fi

##################################
### Post processig data for using  plotband.x ####
## This method doesn't give good data for plotting #### 
####UP####
echo "Post procesing bands up ..."
cat >  $SYSTEM'_fm_bands_up'.in << EOF
 &bands
    outdir='$SCRATCHDIRDOS'
    prefix='$SYSTEM'    
    filband = '${SYSTEM}_fm_sup.bands.dat'
    spin_component=1,
    no_overlap=.true.,

 /
EOF
echo "  running BADS_UP calculation..."

     if [ ! -f ${SYSTEM}_fm_sup.bands.dat ]; then
       $MPI $BANDS_BIN < $SYSTEM'_fm_bands_up'.in >> $SYSTEM'_fm_bands_up'.out
     fi

###### PLOT BANDS UP  ##############
### Example input #######
# ibands.dat   # name_of_file_produced_by_bands.dat
# -6.0 10  Emin, Emax
# sibands.xmgr name_of_grace_output.xmgr
# sibands.ps name_of_ps_output.ps
# 6.255  (EF)
# 1.0, 6.255   (deltaE, reference E (for tics), FE)
echo "where it is running:"
ls ${SYSTEM}_fm_sup.bands.dat
 
echo "Ploting bands up ..."
cat > $SYSTEM'_fm_bands_up'.plotband.in << EOF
${SYSTEM}_fm_sup.bands.dat  
-1, 30
${SYSTEM}_fm_sup.bands.xmgr
${SYSTEM}_fm_sup.bands.ps
$EF
5.0, $EF
EOF

$BANDPLOT_BIN < $SYSTEM'_fm_bands_up'.plotband.in  > $SYSTEM'_fm_bands_up'.plotband.out

#### DOWN #####
echo "Post procesing bands dn ..."
cat >  $SYSTEM'_fm_bands_dn'.in << EOF
 &bands
    outdir='$SCRATCHDIRDOS'
    prefix='$SYSTEM'    
    filband = '${SYSTEM}_fm_sdn.bands.dat'
    spin_component=2,
    no_overlap=.true.,
 /
EOF

      if [ ! -f ${SYSTEM}_fm_sdn.bands.dat ]; then
       $MPI $BANDS_BIN < $SYSTEM'_fm_bands_dn'.in >> $SYSTEM'_fm_bands_dn'.out
      fi

###### PLOT BANDS UP  ##############
echo "Ploting bands down ..."
echo "where it is running:"
ls ${SYSTEM}_fm_sdn.bands.dat

cat > $SYSTEM'_fm_bands_dn'.plotband.in << EOF
${SYSTEM}_fm_sdn.bands.dat  
-1.0, 30
${SYSTEM}_fm_sdn.bands.xmgr
${SYSTEM}_fm_sdn.bands.ps
$EF
5.0, $EF
EOF

$BANDPLOT_BIN < $SYSTEM'_fm_bands_dn'.plotband.in  > $SYSTEM'_fm_bands_dn'.plotband.out


############## nscf for DOS ######################
echo "Setting nscf calculation for PDOS ..."
# non self-consistent calculation
cat >  ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_nscf'.in << EOF
&control
    calculation='nscf'
    restart_mode='from_scratch',    
    pseudo_dir = '$PSEUDO_DIR',
    outdir='$SCRATCHDIRDOS'
    prefix='$SYSTEM',
    tprnfor = .true.
    tstress=.true.
    nstep=50000
    etot_conv_thr=1.0D-5
    forc_conv_thr=1.0D-4
    verbosity='high'
    disk_io='default'

 /
 &system
    ibrav=0, 
    celldm(1)=1,
    nat=128 , ntyp= 1, ecutwfc = $ECUTWFC, ecutrho = $ECUTRHO,        
    occupations='$OCC', $OPTIONS
    nbnd = 88,
    input_DFT='$FUNC'
    london=.true.
    london_rcut = 200
    london_s6 = 0.75

 /
 &electrons
    diagonalization='david'
    electron_maxstep=1000
    conv_thr = 1.0e-8
    mixing_beta = 0.5
 /
 &ions
 /
 &cell
   press=0.0
   press_conv_thr=0.1
 /
ATOMIC_SPECIES
EOF
cat $SPECLIST $COORDS $SYSTEM'_fm_latvec_final'   >> ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_nscf'.in
cat >> ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_nscf'.in << EOF
K_POINTS (automatic) 
EOF
cat $KPOINTSDOS >> ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_nscf'.in

      if [ ! -f ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_nscf'.out ]; then
        $MPI $PW_BIN < ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_nscf'.in $POOL_C >> ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_nscf'.out 2> /dev/null
        echo "Finishing nscf ..."
      fi
##########DOS ################
echo "Post processing DOS .."
cat >  ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_dos'.in << EOF
 &dos
    outdir='$SCRATCHDIRDOS'
    prefix='$SYSTEM'
    fildos='${SYSTEM}_fm.dos'
    DeltaE=0.1,
    ngauss=0
    degauss=0.01
    Emin=-50.0, Emax=50
 /
EOF

echo "  running DOS calculation..."

      if [ ! -f ./$SYSTEM'_fm_PDOSBAND'/${SYSTEM}_fm.dos ]; then
       $MPI $DOSWFC_BIN < ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_dos'.in >> ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_dos'.out
      fi

echo "copying" ${SYSTEM}_fm.dos "..."
if [ -e ${SYSTEM}_fm.dos ]; then
    cp ${SYSTEM}_fm.dos ./$SYSTEM'_fm_PDOSBAND'
    rm ${SYSTEM}_fm.dos
fi

########## PDOS  ################
echo "Post processing PDOS .."
cat > ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_pdos'.in << EOF
 &projwfc
    prefix='$SYSTEM'
    outdir='$SCRATCHDIRDOS'
    filpdos='${SYSTEM}_fm.pdos'
    DeltaE=0.01
    ngauss=0
    degauss=0.01
    Emin=-50.0, Emax=50.0
 /
EOF

echo "  running PDOS calculation..."

$MPI $PROJWFC_BIN < ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm_pdos'.in >> ./$SYSTEM'_fm_PDOSBAND'/$SYSTEM'_fm'.pdos.out 2> /dev/null

cp *.pdos_tot* *.pdos_atm#* ./$SYSTEM'_fm_PDOSBAND'
rm *.pdos_atm#*

echo "Ending DOS & PDOS & BANDS"
fi # End of DOS &  PDOS & BANDS

if [[ -n $SUMDOS ]]; then   # if we want to sum PDOS with sum_states.py
    cd ./$SYSTEM'_fm_PDOSBAND'
    cat  *.dos*  | awk '{print $1-$EF, $2, $3*(-1)}' | sed '1d'  > ${SYSTEM}_fm.dos.shifted.dat
    cat *pdos_tot*  | awk '{print $1-$EF, $2, $3*(-1)}' | sed '1d'  > ${SYSTEM}_fm.pdos.tot.shifted.dat
#   cat ./$SYSTEM'_fm_PDOSBAND'/*pdos_tot*  | awk '{print $1-$EF, $4, $5}' | sed '1d'  > ./$SYSTEM'_fm_PDOSBAND'/${SYSTEM}_fm.pdos.tot.shifted_2.dat
    sum_states.py -o *scf.out* -s *pdos_atm#1\(Fe\)_wfc#3\(d\)* -t "${SYSTEM}_fm" -xr -10 10
    cp sum_dos.out $SYSTEM.3d.sum_dos.out
    cat $SYSTEM.3d.sum_dos.out | awk '{print $1, $2, $3*(-1)}'  > ${SYSTEM}_fm.pdos.3d.sum_dos.dat
    sum_states.py -o *scf.out* -s *pdos_atm#1\(Fe\)_wfc#4\(s\)* -p sum.Fe.d.dos.out -t "BCC Bulk Fe" -xr -10 10
    cp sum_dos.out $SYSTEM.4s.sum_dos.out  
    cat $SYSTEM.4s.sum_dos.out | awk '{print $1, $2, $3*(-1)}'  > ${SYSTEM}_fm.pdos.4s.sum_dos.dat
fi



############ Delete resutls_ folders ######################

if [[ -n $RMV ]]; then # if RMV option is present it removes the temporary files
   echo "Removing" $SYSTEM 
   rm -r *results_$SYSTEM* 
fi

########### Copying the .save file from relaxation here #####

if [[ -n $CPWF ]]; then  # If CPWF is present it copies the $SYSYEM.save file in the current working directory
   echo "Copying .save to the work directory "$SYSTEM results_$SYSTEM'_fm_' 
   cp -r results_$SYSTEM'_fm_'/$SYSYEM.save .
fi
echo " - DONE! -"
