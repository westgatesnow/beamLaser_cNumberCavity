#This bash file is to run multiple simulations with respect to 
#--different transit time tau
#--different average atom number nAtom
#for the beam laser
#using cNumberCavity theory.

iFile=input.txt
nMax1=10
nMax2=1
init_tau=1.0
interval_tau=0.5
init_nAtom=400
interval_nAtom=25

for ((i=0; i<nMax1; i+=1)) do 
  for ((j=0; j<nMax2; j+=1)) do 
    tau=$(echo "$init_tau + $interval_tau * $i" | bc -l)
    nAtom=$(echo "$init_nAtom + $interval_nAtom * $j" | bc -l)
    #dens=$(echo "$nAtom / $tau" | bc -l)
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
      density 100
      rabi 3
      kappa 90
      lambda 1.0
      invT2 0
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
      density 100
      rabi 3
      kappa 90
      lambda 1.0
      invT2 0
      name latt_tau${tau}_g3_k90_den100_l1.0_pz0_dz0_dt0.01" > $iFile
    fi
    ./beamLaser -f $iFile
  done
  number=$((1+$i))
  echo "Run ${number} of" $nMax1
done
