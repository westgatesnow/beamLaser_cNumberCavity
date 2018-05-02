#This bash file is to run multiple simulations with respect to transit time tau
#for the beam laser laser using cNumberCavity theory.

iFile=input.txt
nMax=100
init=0.1
interval=0.1
nAtom=100

for ((i=0; i<nMax; i+=1)) 
do

tau=$(echo "$init + $interval * $i" | bc -l)
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
name pois_tau0${tau}_g2_k40_N50" > $iFile

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
name pois_tau${tau}_g2_k40_N50" > $iFile

fi

./beamLaser -f $iFile

number=$((1+$i))

echo "Run ${number} of" $nMax

done
