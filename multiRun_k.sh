#This bash file is to run multiple simulations with respect to Kappa
#for the beam laser using cNumberCavity theory.

iFile=input.txt
nMax=40
init=-1.0
interval=0.1

for ((i=0; i<nMax; i+=1)) 
do

kappa=$(echo "e(($init + $interval * $i)*l(10))" | bc -l)
exp=$(echo "$init + $interval * $i" | bc -l | awk '{printf "%.1f", $0}')
#$(echo "10+$i*10" | bc -l);

printf "dt 0.001
tmax 20
nStore 1000
nTrajectory 1000
nBin 20
yWall 10
lambda 1.0
deltaZ 0
deltaPz 0
transitTime 1.0
density 100
rabi 3
kappa $kappa
invT2 0
controlType kappa
name dt0.001_dZ0_dPz0_tau1.0_nBin20_dens100_g3_kexp$exp" > $iFile

./beamLaser -f $iFile

number=$((1+$i))

echo "Run ${number} of" $nMax

done
