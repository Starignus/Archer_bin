#!/bin/sh
#
####################################################################
#
# Script that calcuates the relaxation steps, the total time in 
# seconds spend in the whole relaxation, and the average tiem in 
# seconds for relaxations step.
#
####################################################################
#Author: Dr. Ariadna Blanca Romero
#        Postdoctoral Research Associate
#        Imperial College London
#        Thomas Young Centre-Chemestry
#        ariadna@starignus.com or starignus@gmail.com
#        https://github.com/Starignus
####################################################################

grep 'total cpu time spent' *out* | awk '{ print $9 }' | awk 'NR-1{print $0-p}{p=$0}' > TIME
SCFsteps=`grep 'acc'  *out* | wc -l`
sumtimesteps=`cat TIME | awk '{ sum+=$1} END {print sum}'`
averagedtime=$(echo " $sumtimesteps/$SCFsteps" | bc -l)
echo "Steps: $SCFsteps"
echo "Sum of seconds: $sumtimesteps"
echo "Average time in the scf per step:$averagedtime"
