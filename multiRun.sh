#This bash file is to run multiple simulations with respect to 
#--different transit time tau
#--different average atom number nAtom
#for the beam laser
#using cNumberCavity theory.

iFile=input.txt
nMax1=1
nMax2=4
init_tau=0.65
interval_tau=0.5
init_nAtom=400
interval_nAtom=25

for ((i=0; i<nMax1; i+=1)) do 
  for ((j=0; j<nMax2; j+=1)) do 
    tau=$(echo "$init_tau + $interval_tau * $i" | bc -l)
    nAtom=$(echo "$init_nAtom + $interval_nAtom * $j" | bc -l)
    dens=$(echo "$nAtom / $tau" | bc -l)
    if (( $(bc <<< "$tau < 1") )) 
    then
      printf "dt 0.01
      tmax 100
      nTrajectory 1000
      nstore 1000
      yWall 10
      sigmaXX 0
      sigmaXZ 0
      transitTime $tau
      sigmaPX 0
      sigmaPY 0
      sigmaPZ 0
      density $dens
      rabi 2
      kappa 40
      name pois_tau0${tau}_g2_k40_N${nAtom}" > $iFile
    else
      printf "dt 0.01
      tmax 100
      nTrajectory 1000
      nstore 1000
      yWall 10
      sigmaXX 0
      sigmaXZ 0
      transitTime $tau
      sigmaPX 0
      sigmaPY 0
      sigmaPZ 0
      density $dens
      rabi 2
      kappa 40
      name pois_tau${tau}_g2_k40_N${nAtom}" > $iFile
    fi
    ./beamLaser -f $iFile
  done
  number=$((1+$i))
  echo "Run ${number} of" $nMax1
done
