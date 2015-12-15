#!/bin/sh
SYSTEM=`pwd`
SYSTEM=`basename "$SYSTEM"`
echo $SYSTEM
LENGTH=`echo -n $SYSTEM | wc -c`
echo $LENGTH
LENGTHNEW=`echo $(($LENGTH - 4))`
echo $LENGTHNEW
echo $SYSTEM
SYSTEM_PATH=`echo $SYSTEM | cut -c -${LENGTHNEW}`
echo $SYSTEM_PATH

for i in 1 2 3 4 5 A
  do
    if [[ ! -d L$i ]]; then
       echo "Creating folders..."
       mkdir L$i
   fi
done

# # look for empty dir 
if [ -z "`find L1 -type f`" ]; then 
  echo "Copying dos to L1..."
  cp ${SYSTEM_PATH}.pdos_atm#1\(* L1/
  cp ${SYSTEM_PATH}.pdos_atm#2\(* L1/
  cp ${SYSTEM_PATH}.pdos_atm#3\(* L1/
  cp ${SYSTEM_PATH}.pdos_atm#4\(* L1/
fi

if [ -z "`find L2 -type f`" ]; then
  echo "Copying dos to L2..."
  cp ${SYSTEM_PATH}.pdos_atm#5\(* L2/
  cp ${SYSTEM_PATH}.pdos_atm#6\(* L2/
  cp ${SYSTEM_PATH}.pdos_atm#7\(* L2/
  cp ${SYSTEM_PATH}.pdos_atm#8\(* L2/
fi

if [ -z "`find L3 -type f`" ]; then
  echo "Copying dos to L3..."
  cp ${SYSTEM_PATH}.pdos_atm#9\(* L3/
  cp ${SYSTEM_PATH}.pdos_atm#10\(* L3/
  cp ${SYSTEM_PATH}.pdos_atm#11\(* L3/
  cp ${SYSTEM_PATH}.pdos_atm#12\(* L3/
fi

if [ -z "`find L4 -type f`" ]; then
  echo "Copying dos to L4..."
  cp ${SYSTEM_PATH}.pdos_atm#13\(* L4/
  cp ${SYSTEM_PATH}.pdos_atm#14\(* L4/
  cp ${SYSTEM_PATH}.pdos_atm#15\(* L4/
  cp ${SYSTEM_PATH}.pdos_atm#16\(* L4/
fi

if [ -z "`find L5 -type f`" ]; then
  echo "Copying dos to L5..."
  cp ${SYSTEM_PATH}.pdos_atm#17\(* L5/
  cp ${SYSTEM_PATH}.pdos_atm#18\(* L5/
  cp ${SYSTEM_PATH}.pdos_atm#19\(* L5/
  cp ${SYSTEM_PATH}.pdos_atm#20\(* L5/
fi

if [ -z "`find LA -type f`" ]; then
  echo "Copying dos to LA..."
  cp ${SYSTEM_PATH}.pdos_atm#21\(* LA/
  cp ${SYSTEM_PATH}.pdos_atm#22\(* LA/
  cp ${SYSTEM_PATH}.pdos_atm#23\(* LA/
fi

# Fermi energy
EF=`cat ${SYSTEM_PATH}_scf_EF`
SCF_out="../${SYSTEM_PATH}_scf.out"
echo $SCF_out

cd ./L1
echo "PDOS L1..."
sum_states.py -s "*Fe*wfc*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_tot_L1.dat

sum_states.py -s "*Fe*wfc*d*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_d_tot_L1.dat

sum_states.py -s "*Fe*wfc*s*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_s_tot_L1.dat

sum_states.py -s "*Fe*wfc*p*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_p_tot_L1.dat

cat Fe_tot_L1.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_tot_L1_shifted.dat
cat Fe_d_tot_L1.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_d_tot_L1_shifted.dat
cat Fe_s_tot_L1.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_s_tot_L1_shifted.dat
cat Fe_p_tot_L1.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_p_tot_L1_shifted.dat
cd ..

cd ./L2
pwd
echo "PDOS L2..."
sum_states.py -s "*Fe*wfc*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_tot_L2.dat

sum_states.py -s "*Fe*wfc*d*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_d_tot_L2.dat

sum_states.py -s "*Fe*wfc*s*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_s_tot_L2.dat

sum_states.py -s "*Fe*wfc*p*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_p_tot_L2.dat

cat Fe_tot_L2.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_tot_L2_shifted.dat
cat Fe_d_tot_L2.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_d_tot_L2_shifted.dat
cat Fe_s_tot_L2.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_s_tot_L2_shifted.dat
cat Fe_p_tot_L2.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_p_tot_L2_shifted.dat
cd ..

cd ./L3
pwd
echo "PDOS L3..."
sum_states.py -s "*Fe*wfc*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_tot_L3.dat

sum_states.py -s "*Fe*wfc*d*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_d_tot_L3.dat

sum_states.py -s "*Fe*wfc*s*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_s_tot_L3.dat

sum_states.py -s "*Fe*wfc*p*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_p_tot_L3.dat

cat Fe_tot_L3.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_tot_L3_shifted.dat
cat Fe_d_tot_L3.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_d_tot_L3_shifted.dat
cat Fe_s_tot_L3.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_s_tot_L3_shifted.dat
cat Fe_p_tot_L3.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_p_tot_L3_shifted.dat
cd ..

cd ./L4
echo "PDOS L4..."
sum_states.py -s "*Fe*wfc*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_tot_L4.dat

sum_states.py -s "*Fe*wfc*d*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_d_tot_L4.dat

sum_states.py -s "*Fe*wfc*s*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_s_tot_L4.dat

sum_states.py -s "*Fe*wfc*p*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_p_tot_L4.dat

cat Fe_tot_L4.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_tot_L4_shifted.dat
cat Fe_d_tot_L4.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_d_tot_L4_shifted.dat
cat Fe_s_tot_L4.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_s_tot_L4_shifted.dat
cat Fe_p_tot_L4.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_p_tot_L4_shifted.dat
cd ..

cd ./L5
echo "PDOS L5..."
sum_states.py -s "*Fe*wfc*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_tot_L5.dat

sum_states.py -s "*Fe*wfc*d*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_d_tot_L5.dat

sum_states.py -s "*Fe*wfc*s*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_s_tot_L5.dat

sum_states.py -s "*Fe*wfc*p*" -o $SCF_out -xr -100 100
mv sum_dos.out Fe_p_tot_L5.dat

cat Fe_tot_L5.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_tot_L5_shifted.dat
cat Fe_d_tot_L5.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_d_tot_L5_shifted.dat
cat Fe_s_tot_L5.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_s_tot_L5_shifted.dat
cat Fe_p_tot_L5.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > Fe_p_tot_L5_shifted.dat
cd ..

cd ./LA
echo "PDOS LA..."

if [[ ! -d C_states ]]; then
   mkdir C_states
fi

mv ${SYSTEM_PATH}.pdos_atm#21\(* C_states
cd C_states

sum_states.py -s "*C*wfc*" -o ../$SCF_out -xr -100 100
mv sum_dos.out C_tot_LA.dat

sum_states.py -s "*C*wfc*s*" -o ../$SCF_out -xr -100 100
mv sum_dos.out C_s_tot_LA.dat

sum_states.py -s "*C*wfc*p*" -o ../$SCF_out -xr -100 100
mv sum_dos.out C_p_tot_LA.dat

cat C_tot_LA.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > C_tot_LA_shifted.dat
cat C_s_tot_LA.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > C_s_tot_LA_shifted.dat
cat C_p_tot_LA.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > C_p_tot_LA_shifted.dat
cd ..

sum_states.py -s "*O*wfc*" -o $SCF_out -xr -100 100
mv sum_dos.out O_tot_LA.dat

sum_states.py -s "*O*wfc*s*" -o $SCF_out -xr -100 100
mv sum_dos.out O_s_tot_LA.dat

sum_states.py -s "*O*wfc*p*" -o $SCF_out -xr -100 100
mv sum_dos.out O_p_tot_LA.dat

cat O_tot_LA.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > O_tot_LA_shifted.dat
cat O_s_tot_LA.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > O_s_tot_LA_shifted.dat
cat O_p_tot_LA.dat | awk '{print $1, $2, $3*(-1)}' |  sed '1d' > O_p_tot_LA_shifted.dat
cd ..

echo "Finish.."
