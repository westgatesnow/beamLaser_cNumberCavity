#This bash file is to run multiple simulations with respect to 
#--different transit time tau
#--different average atom number nAtom
#for the beam laser
#using cNumberCavity theory.

iFile=input.txt
nMax1=10
nMax2=10
init_tau=0.5
interval_tau=0.1
init_nAtom=50
interval_nAtom=50

for ((i=0; i<nMax1; i+=1)) do 
  for ((j=0; j<nMax2; j+=1)) do 
    tau=$(echo "$init_tau + $interval_tau * $i" | bc -l)
    nAtom=$(echo "$init_nAtom + $interval_nAtom * $j" | bc -l)
    dens=$(echo "$nAtom / $tau" | bc -l)
    if (( $(bc <<< "$tau < 1") )) 
    then
      printf "dt 0.005
      tmax 20
      nStore 1000
      nTrajectory 1000
      nBin 20
      yWall 1.0
      lambda 1.0
      deltaZ 0.5
      deltaPz 0.5
      transitTime $tau
      density $dens
      rabi 3
      kappa 90
      invT2 0
      controlType contour_dpz0.5
      name pois_dt0.005_dZ0.5_dPz0.5_tau0${tau}_nBin20_nAtom${nAtom}_g3_k90_yWall1.0" > $iFile
    else
      printf "dt 0.005
      tmax 20
      nStore 1000
      nTrajectory 1000
      nBin 20
      yWall 1.0
      lambda 1.0
      deltaZ 0.5
      deltaPz 0.5
      transitTime $tau
      density $dens
      rabi 3
      kappa 90
      invT2 0
      controlType contour_dpz0.5
      name pois_dt0.005_dZ0.5_dPz0.5_tau${tau}_nBin20_nAtom${nAtom}_g3_k90_yWall1.0" > $iFile
    fi
    ./beamLaser -f $iFile
  done
  number=$((1+$i))
  echo "Run ${number} of" $nMax1
done
