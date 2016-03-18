#!/bin/sh

#
####################################################################
#
# Script that check the status of a job in the queue in Archer 
#
####################################################################
#Author: Dr. Ariadna Blanca Romero
#        Postdoctoral Research Associate
#        Imperial College London
#        Thomas Young Centre-Chemestry
#        ariadna@starignus.com or starignus@gmail.com
#        https://github.com/Starignus
####################################################################

while getopts x: opt
do
   case "$opt" in
      x) search=$OPTARG;;
      \?) echo """   resend usage:
    -x  qstat ID
    -h  print this help 
   Example:time  qstatlsit.sh -x "2342342" 
 """
     exit 1  ;;
   esac
done

date +"%D"
date +"%T"
pwd
qstat | awk '{print $1, $5}' > list 
tail -n+3 list >list2
rm list
g=`grep "$search"  list2`  
S=`echo $g | awk '{print $2}'`

if  [[ -z $g ]]; then 
   echo "there is not such job in the queue"
 else
   echo "the job" $g "is in the queue with status" $S 
fi
rm list2
