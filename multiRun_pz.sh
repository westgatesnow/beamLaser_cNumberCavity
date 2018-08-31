#This bash file is to run multiple simulations with respect to 
#--different pz
#--different density
#for the beam laser
#using cNumberCavity theory.

iFile=input.txt
nMax1=10
nMax2=1
init_pz=1.0
interval_pz=1.0
init_dens=100
interval_dens=100

for ((i=0; i<nMax1; i+=1)) do 
  for ((j=0; j<nMax2; j+=1)) do 
    pz=$(echo "$init_pz + $interval_pz * $i" | bc -l)
    dens=$(echo "$init_dens + $interval_dens * $j" | bc -l)
    if (( $(bc <<< "$pz < 1") )) 
    then
      printf "dt 0.005
      tmax 20
      nStore 1000
      nTrajectory 1000
      nBin 30
      yWall 1.0
      lambda 1.0
      deltaZ 0.5
      deltaPz $pz
      transitTime 1.0
      density $dens
      rabi 3
      kappa 90
      invT2 0
      controlType dens100
      name dt0.005_dZ0.5_dPz0${pz}_tau1.0_nBin30_dens${dens}_g3_k90_yWall1.0" > $iFile
    else
      printf "dt 0.005
      tmax 20
      nStore 1000
      nTrajectory 1000
      nBin 30
      yWall 1.0
      lambda 1.0
      deltaZ 0.5
      deltaPz $pz
      transitTime 1.0
      density $dens
      rabi 3
      kappa 90
      invT2 0
      controlType dens100
      name dt0.005_dZ0.5_dPz${pz}_tau1.0_nBin30_dens${dens}_g3_k90_yWall1.0" > $iFile
    fi
    ./beamLaser -f $iFile
  done
  number=$((1+$i))
  echo "Run ${number} of" $nMax1
done
