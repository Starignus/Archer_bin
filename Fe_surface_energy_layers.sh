#!/bin/sh
#
# Convergence test for charge density/dual dependence
# once is known the wavefunction cutoff
####################################################################


#### System variables and files
SYSTEM=`pwd`
SYSTEM=${SYSTEM##*/}
VERSION=1.0


###########Options###################
while getopts x:m:a:p:g:b:c:d: opt
do
   case "$opt" in
      x) FROM=$OPTARG;;
      m) TO=$OPTARG;;
      a) AREA=$OPTARG;;
      p) PO="yes";;
      g) POP="yes";;
      b) POPP="yes";;
      c) SLOPERY=$OPTARG;;
      d) POPLOT="yes";;
      \?) echo """   resend usage:
    -x  From which # of layer you want to start (should be one after the smallest number of layers)
    -m  Up to how many layers
    -a  Surface area in Ang^2
    -p  Plot data if it is activated ("yes")
    -g  This option will not run the first part ("yes")
    -b  Surface energy calculaterd when the bulk energy comes from the slope NlayersVsE ("yes")
    -c  Slope in Ry
    -d  Plot data with Ebulk from fit
   Example:(110)
           time Fe_surface_energy_layers.sh -x "4" -m "16" -a"5.673182" -p "yes"
           time Fe_surface_energy_layers.sh -x "4" -m "16" -a"5.673182" -g "yes" -b "yes" -c "-329.2715" (for running Ebulk from slpe)
           time Fe_surface_energy_layers.sh -d "yes" -g "yes" (Option just to plot Ry results with E fit and no efit together)
           (100)
           time Fe_surface_energy_layers.sh -x "4" -m "16" -a "8.02309024" -p "yes"
           time Fe_surface_energy_layers.sh -x "4" -m "16" -a "8.02309024" -g "yes" -b "yes" -c "-329.2717" (for running Ebulk from slpe)
           time Fe_surface_energy_layers.sh -d "yes" -g "yes"
           (111)
           time Fe_surface_energy_layers.sh -x "4" -m "16" -a "13.8964" -p "yes"
           time Fe_surface_energy_layers.sh -x "4" -m "16" -a "13.8964" -g "yes" -b "yes" -c "-329.2712" (for running Ebulk from slpe)
           time Fe_surface_energy_layers.sh -d "yes" -g "yes"

 """
     exit 1  ;;
   esac
done

echo "    Surfaceenergy  Layers (110,100,111) Fe slab." $VERSION
echo "    -------------------" 

if [[ -z $POP ]]; then #If this option is given this part will not run
  ## Surface energies##
  #Area=5.67331815 Ang^2
  cat  NVsEnergyRy.out | awk '{a[NR]=$2}{b[NR]=$1}{c=NR}END {for(i=1;i<=c;i++) print (a[i+1]-b[i+1]*(a[i+1]-a[i])/2)/2}' | sed '$d'  > SurfaceenergyRy
  cat  NVsEnergyeV.out  | awk '{a[NR]=$2}{b[NR]=$1}{c=NR}END {for(i=1;i<=c;i++) print (a[i+1]-b[i+1]*(a[i+1]-a[i])/2)/2}' | sed '$d' > SurfaceenergyeV
  cat  NVsEnergyRy.out |  awk -v AR="$AREA" '{a[NR]=$2}{b[NR]=$1}{c=NR}END {for(i=1;i<=c;i++) print ((a[i+1]-b[i+1]*(a[i+1]-a[i])/2)/(2*AR))}' | sed '$d'  > SurfaceenergyRyAngs2
  cat  NVsEnergyeV.out  | awk  -v AR="$AREA" '{a[NR]=$2}{b[NR]=$1}{c=NR}END {for(i=1;i<=c;i++) print ((a[i+1]-b[i+1]*(a[i+1]-a[i])/2)/(2*AR))}' | sed '$d'  > SurfaceenergyeVAngs2

  layers_list=`seq $FROM 2 $TO`
  array_layers=(${layers_list// / })
  echo "List of layers" $layers_list
  for j in `seq 0 ${#array_layers[@]}`
    do
    enrgy_layerRy=`cat SurfaceenergyRy | awk '{print $1 }'`
    array_energy_layerRy=(${enrgy_layerRy// / })
    #echo ${array_layers[$j]} ${array_energy_layerRy[$j]} 
    #echo "uno "
    echo ${array_layers[$j]} ${array_energy_layerRy[$j]} >> LayersVsSurfaceenergyRy # Nlayres, E(Ry) 
    enrgy_layereV=`cat SurfaceenergyeV | awk '{print $1 }'`
    array_energy_layereV=(${enrgy_layereV// / })
    #echo ${array_layers[$j]} ${array_energy_layereV[$j]}
    #echo "dos "
    echo ${array_layers[$j]} ${array_energy_layereV[$j]} >> LayersVsSurfaceenergyeV  # Nlayres, E(eV)
    enrgy_layerRyarea=`cat SurfaceenergyRyAngs2 | awk '{print $1 }'`
    array_energy_layerRyarea=(${enrgy_layerRyarea// / })
    # echo ${array_layers[$j]} ${array_energy_layerRyarea[$j]}
    # echo  "tres  "
    echo ${array_layers[$j]} ${array_energy_layerRyarea[$j]} >> LayersVsSurfaceenergyRyarea  # Nlayres, E(Ry) area Ang (for J/m2)
    enrgy_layereVarea=`cat SurfaceenergyeVAngs2 | awk '{print $1 }'`
    array_energy_layereVarea=(${enrgy_layereVarea// / })
    #echo ${array_layers[$j]} ${array_energy_layereV[$j]} 
    #echo "cuatro "
    echo ${array_layers[$j]} ${array_energy_layereV[$j]} >> LayersVsSurfaceenergyeVarea   # Nlayres, E(eV) area Ang
  done
fi #END POP 

if [[ -n $POPP ]]; then # option to estimate surface energy once we know the slope of N Vs E for estimating the Bulk energy
  echo $SLOPERY $SLOPEEV
  cat  NVsEnergyRy.out | awk -v EBulkRy="$SLOPERY" '{a[NR]=$2}{b[NR]=$1}{c=NR}END {for(i=1;i<=c;i++) print (a[i+1]-b[i+1]*EBulkRy)/2}' | sed '$d'  > SurfaceenergyRy_Ebulkslope  #Energy in Ry with slope
#  cat  NVsEnergyeV.out  | awk -v  EBulkeV="$SLOPEEV" '{a[NR]=$2}{b[NR]=$1}{c=NR}END {for(i=1;i<=c;i++) print (a[i+1]-b[i+1]*EBulkeV)/2}' | sed '$d' >> SurfaceenergyeV_Ebulkslope
  cat  NVsEnergyRy.out |  awk -v AR="$AREA" -v  EBulkRy="$SLOPERY" '{a[NR]=$2}{b[NR]=$1}{c=NR}END {for(i=1;i<=c;i++) print ((a[i+1]-b[i+1]*EBulkRy)/(2*AR))}' | sed '$d'  > SurfaceenergyRyAngs2_Ebulkslope #Energy in Ry with slope area Ang
#  cat  NVsEnergyeV.out  | awk  -v AR="$AREA" -v  EBulkeV="$SLOPEEV" '{a[NR]=$2}{b[NR]=$1}{c=NR}END {for(i=1;i<=c;i++) print ((a[i+1]-b[i+1]*EBulkeV)/(2*AR))}' | sed '$d' | sed '$d' >> SurfaceenergyeVAngs2_Ebulkslope

  layers_list=`seq $FROM 2 $TO`
  array_layers=(${layers_list// / })
  echo "List of layers" $layers_list
  for j in `seq 0 ${#array_layers[@]}`
    do
    enrgy_layerRy=`cat SurfaceenergyRy_Ebulkslope | awk '{print $1 }'`
    array_energy_layerRy=(${enrgy_layerRy// / })
    echo ${array_layers[$j]} ${array_energy_layerRy[$j]} >> LayersVsSurfaceenergyRy_Ebulkslope # Nlayres, E(Ry) with slope
    cat LayersVsSurfaceenergyRy_Ebulkslope | awk '{print $1, $2*13.605691}' | sed '$d' > LayersVsSurfaceenergyeV_Ebulkslope  # Nlayres, E(eV) with slope
    #enrgy_layereV=`cat SurfaceenergyeV_Ebulkslope | awk '{print $1 }'`
    #array_energy_layereV=(${enrgy_layereV// / })
    #echo ${array_layers[$j]} ${array_energy_layereV[$j]} >> LayersVsSurfaceenergyeV_Ebulkslope
    enrgy_layerRyarea=`cat SurfaceenergyRyAngs2_Ebulkslope | awk '{print $1 }'`  
    array_energy_layerRyarea=(${enrgy_layerRyarea// / })
    echo ${array_layers[$j]} ${array_energy_layerRyarea[$j]} >> LayersVsSurfaceenergyRyarea_Ebulkslope # Nlayres, E(Ry) with slope area Ang (for J/m2)
    cat LayersVsSurfaceenergyRyarea_Ebulkslope  | awk '{print $1, $2*13.605691}' | sed '$d' >  LayersVsSurfaceenergyeVarea_Ebulkslope  # Nlayres, E(eV) with slope area Ang
    # enrgy_layereVarea=`cat SurfaceenergyeVAngs2_Ebulkslope | awk '{print $1 }'`
    # array_energy_layereVarea=(${enrgy_layereVarea// / })
    # echo ${array_layers[$j]} ${array_energy_layereV[$j]} >> LayersVsSurfaceenergyeVarea_Ebulkslope
   
    #Tranformation to J/m2
    cat LayersVsSurfaceenergyRyarea | awk '{print $1, $2*217.98741}' > LayersVsSurfaceenergyRyarea_Jm2
    cat LayersVsSurfaceenergyRyarea_Ebulkslope | awk '{print $1, $2*217.98741}' > LayersVsSurfaceenergyRyarea_Ebulkslope_Jm2 
  done
  cat LayersVsSurfaceenergyRyarea_Jm2 | awk '{ print $2 }' | awk 'NR-1{print $0-p}{p=$0}' > Surface.energy.diferenceJm2
  cat LayersVsSurfaceenergyRyarea_Ebulkslope_Jm2 | awk '{ print $2 }' | awk 'NR-1{print $0-p}{p=$0}' > Surface.energy.diference_Ebulkslope_Jm2
  cat LayersVsSurfaceenergyRyarea_Jm2 | awk '{a[NR]=$2}END {for(i=1;i<NR-1;i++)print a[i]-a[NR-1]}' > Surface.energy.diference_finalavalue_Jm2
  cat LayersVsSurfaceenergyRyarea_Ebulkslope_Jm2 | awk '{a[NR]=$2}END {for(i=1;i<NR-1;i++)print a[i]-a[NR-1]}' > Surface.energy.diference_finalavalue_Ebulkslope_Jm2

fi #end POPP 

if [[ -n $PO ]]; then  
  # Plotting results with gnuplot

  gnuplot_bin=`which gnuplot 2>/dev/null`
  if [ "$gnuplot_bin" != "" ]; then
     $gnuplot_bin -persist -e "set xlabel 'Number of Layers'; set ylabel 'Surface Energy (Ry)'; p 'LayersVsSurfaceenergyRy' u 1:2 w linespoints"
     $gnuplot_bin -e "set terminal png; set output 'fig.surface.energyRy.out.png'; set xlabel 'Number of Layers'; set ylabel 'Surface Energy (Ry)'; p 'LayersVsSurfaceenergyRy' u 1:2 w linespoints"
     $gnuplot_bin -persist -e "set xlabel 'Number of Layers'; set ylabel 'Surface Energy (eV)'; p 'LayersVsSurfaceenergyeV' u 1:2 w linespoints"
     $gnuplot_bin -e "set terminal png; set output 'fig.surface.energyeV.out.png'; set xlabel 'Number of Layers'; set ylabel 'Surface Energy (eV)'; p 'LayersVsSurfaceenergyeV' u 1:2 w linespoints"
     $gnuplot_bin -persist -e "set xlabel 'Number of Layers'; set ylabel 'Surface Energy (Ry/Ang^2)'; p 'LayersVsSurfaceenergyRyarea' u 1:2 w linespoints"
     $gnuplot_bin -e "set terminal png; set output 'fig.surface.energyRyarea.out.png'; set xlabel 'Number of Layers'; set ylabel 'Surface Energy (Ry/Ang^2)'; p 'LayersVsSurfaceenergyRyarea' u 1:2 w linespoints"
     $gnuplot_bin -persist -e "set xlabel 'Number of Layers'; set ylabel 'Surface Energy (eV/Ang^2)'; p 'LayersVsSurfaceenergyeVarea' u 1:2 w linespoints"
     $gnuplot_bin -e "set terminal png; set output 'fig.surface.energyeVarea.out.png'; set xlabel 'Number of Layers'; set ylabel 'Surface Energy (eV/Ang^2)'; p 'LayersVsSurfaceenergyeVarea' u 1:2 w linespoints"
  else
  d        $ECHO "No gnuplot in PATH. Results not plotted."
  fi
fi #End PO

if [[ -n $POPLOT ]]; then
  # Plotting Ebulk with gnuplot

  gnuplot_bin=`which gnuplot 2>/dev/null`
  if [ "$gnuplot_bin" != "" ]; then
     $gnuplot_bin -persist -e "set xlabel 'Number of Layers'; set ylabel 'Surface Energy (Ry)'; p 'LayersVsSurfaceenergyRy' u 1:2 w linespoints,'LayersVsSurfaceenergyRy_Ebulkslope' u 1:2 w linespoints"
     $gnuplot_bin -e "set terminal png; set output 'fig.surface.energyRy_Ebulk.out.png'; set xlabel 'Number of Layers'; set ylabel 'Surface Energy (Ry)'; p 'LayersVsSurfaceenergyRy' u 1:2 w linespoints,'LayersVsSurfaceenergyRy_Ebulkslope' u 1:2 w linespoints"
     $gnuplot_bin -persist -e "set xlabel 'Number of Layers'; set ylabel 'Surface Energy (Ry/Ang^2)';  p 'LayersVsSurfaceenergyRyarea' u 1:2 w linespoints,'LayersVsSurfaceenergyRyarea_Ebulkslope' u 1:2 w linespoints"
     $gnuplot_bin -e "set terminal png; set output 'fig.surface.energyRyarea.out_Ebulk.png'; set xlabel 'Number of Layers'; set ylabel 'Surface Energy (Ry/Ang^2)';  p 'LayersVsSurfaceenergyRyarea' u 1:2 w linespoints,'LayersVsSurfaceenergyRyarea_Ebulkslope' u 1:2 w linespoints "
     $gnuplot_bin -persist -e "set xlabel 'Number of Layers'; set ylabel 'Surface Energy (J/m^2)'; p 'LayersVsSurfaceenergyRyarea_Jm2' u 1:2 w linespoints,'LayersVsSurfaceenergyRyarea_Ebulkslope_Jm2' u 1:2 w linespoints"
     $gnuplot_bin -e "set terminal png; set output 'fig.surface.energyRy_Ebulk_Jm2.out.png'; set xlabel 'Number of Layers'; set ylabel 'Surface Energy (J/m^2)'; p 'LayersVsSurfaceenergyRyarea_Jm2' u 1:2 w linespoints,'LayersVsSurfaceenergyRyarea_Ebulkslope_Jm2' u 1:2 w linespoints"
  else
  d        $ECHO "No gnuplot in PATH. Results not plotted."
  fi
fi #End POPLOT


echo " - DONE! -" 

