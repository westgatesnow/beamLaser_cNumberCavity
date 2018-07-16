#This bash file is to run multiple simulations with respect to 
#--different pz
#--different density
#for the beam laser
#using cNumberCavity theory.

iFile=input.txt
nMax1=10
nMax2=1
init_pz=0.1
interval_pz=0.1
init_dens=100
interval_dens=100

for ((i=0; i<nMax1; i+=1)) do 
  for ((j=0; j<nMax2; j+=1)) do 
    pz=$(echo "$init_pz + $interval_pz * $i" | bc -l)
    dens=$(echo "$init_dens + $interval_dens * $j" | bc -l)
    if (( $(bc <<< "$pz < 1") )) 
    then
      printf "dt 0.01
      tmax 100
      nTrajectory 1000
      nstore 1000
      yWall 10
      sigmaXX 0
      sigmaXZ 0
      transitTime 1.0
      sigmaPX 0
      sigmaPY 0
      sigmaPZ $pz
      density $dens
      rabi 3
      kappa 90
      lambda 1.0
      invT2 0
      name latt_tau1.0_g3_k90_dens${dens}_l1.0_pz0${pz}_dz0.25" > $iFile
    else
      printf "dt 0.01
      tmax 100
      nTrajectory 1000
      nstore 1000
      yWall 10
      sigmaXX 0
      sigmaXZ 0
      transitTime 1.0
      sigmaPX 0
      sigmaPY 0
      sigmaPZ $pz
      density $dens
      rabi 3
      kappa 90
      lambda 1.0
      invT2 0
      name latt_tau1.0_g3_k90_dens${dens}_l1.0_pz${pz}_dz0.25" > $iFile
    fi
    ./beamLaser -f $iFile
  done
  number=$((1+$i))
  echo "Run ${number} of" $nMax1
done
