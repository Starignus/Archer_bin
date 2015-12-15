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
      g) NLAYER=$OPTARG;;
      p) PO="yes";;
      b) POPP="yes";;
      c) SLOPERY=$OPTARG;;
      d) POPLOT="yes";;
      \?) echo """   resend usage:
    -x  Set alwasy to 1
    -m  Up to how many calculations of kpoints or vaccum (e.g 3,4) 
    -a  Surface area in Ang^2
    -g  Number of layers in our slab (ATOMS)
    -p  Plot data if it is activated ("yes")
    -b  Surface energy calculaterd when the bulk energy comes from the slope NlayersVsE ("yes")
    -c  Slope in Ry
    -d  Plot data with Ebulk from fit
   Example:(110)
           time Fe_surface_energy_layers_vacc_kpoints.sh -x "1" -m "6" -a"5.673182" -g "20" -b "yes" -c "-329.2715" (for running Ebulk from slpe)
           time Fe_surface_energy_layers.sh -d "yes" -g "yes" (Option just to plot Ry results with E fit and no efit together)
           (100)
           time Fe_surface_energy_layers.sh -x "1" -m "6" -a "8.02309024" -g "20" -b "yes" -c "-329.2717" (for running Ebulk from slpe)
           time Fe_surface_energy_layers.sh -d "yes" -g "yes"
           (111)
           time Fe_surface_energy_layers.sh -x "1" -m "5" -a "13.8964" -g "18" -b "yes" -c "-329.2712" (for running Ebulk from slpe)
           time Fe_surface_energy_layers.sh -d "yes" -g "yes"

 """
     exit 1  ;;
   esac
done

echo "    Surfaceenergy  Layers (110,100,111) Fe slab." $VERSION
echo "    -------------------" 

if [[ -n $POPP ]]; then # option to estimate surface energy once we know the slope of N Vs E for estimating the Bulk energy
  echo $SLOPERY $SLOPEEV
   cat  NVsEnergyRy.out | awk -v EBulkRy="$SLOPERY" -v Nlay="$NLAYER" '{a[NR]=$1}{c=NR}END {for(i=1;i<=c;i++) print (a[i]-Nlay*EBulkRy)/2}'  > SurfaceenergyRy_Ebulkslope  #Energy in Ry with slope
#  cat  NVsEnergyeV.out  | awk -v  EBulkeV="$SLOPEEV" '{a[NR]=$2}{b[NR]=$1}{c=NR}END {for(i=1;i<=c;i++) print (a[i+1]-b[i+1]*EBulkeV)/2}' | sed '$d' >> SurfaceenergyeV_Ebulkslope
  cat  NVsEnergyRy.out |  awk -v AR="$AREA" -v  EBulkRy="$SLOPERY" -v Nlay="$NLAYER" '{a[NR]=$1}{c=NR}END {for(i=1;i<=c;i++) print ((a[i]-Nlay*EBulkRy)/(2*AR))}'  > SurfaceenergyRyAngs2_Ebulkslope #Energy in Ry with slope area Ang
#  cat  NVsEnergyeV.out  | awk  -v AR="$AREA" -v  EBulkeV="$SLOPEEV" '{a[NR]=$2}{b[NR]=$1}{c=NR}END {for(i=1;i<=c;i++) print ((a[i+1]-b[i+1]*EBulkeV)/(2*AR))}' | sed '$d' | sed '$d' >> SurfaceenergyeVAngs2_Ebulkslope

  vackpints_list=`seq $FROM 2 $TO`
  vackpints_layers=(${layers_list// / })
  echo "List of layers" $vackpints_layers
  for j in `seq 0 ${#vackpints_layers[@]}`
    do
    enrgy_layerRy=`cat SurfaceenergyRy_Ebulkslope | awk '{print $1 }'`
    array_energy_layerRy=(${enrgy_layerRy// / })
    echo ${vackpints_layers[$j]} ${array_energy_layerRy[$j]} >> LayersVsSurfaceenergyRy_Ebulkslope # Nlayres, E(Ry) with slope
    cat LayersVsSurfaceenergyRy_Ebulkslope | awk '{print $1, $2*13.605691}' | sed '$d' > LayersVsSurfaceenergyeV_Ebulkslope  # Nlayres, E(eV) with slope
    enrgy_layerRyarea=`cat SurfaceenergyRyAngs2_Ebulkslope | awk '{print $1 }'`  
    array_energy_layerRyarea=(${enrgy_layerRyarea// / })
    echo ${vackpints_layers[$j]} ${array_energy_layerRyarea[$j]} >> LayersVsSurfaceenergyRyarea_Ebulkslope # Nlayres, E(Ry) with slope area Ang (for J/m2)
    cat LayersVsSurfaceenergyRyarea_Ebulkslope  | awk '{print $1, $2*13.605691}' | sed '$d' >  LayersVsSurfaceenergyeVarea_Ebulkslope  # Nlayres, E(eV) with slope area Ang
    # enrgy_layereVarea=`cat SurfaceenergyeVAngs2_Ebulkslope | awk '{print $1 }'`
    # array_energy_layereVarea=(${enrgy_layereVarea// / })
    
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

