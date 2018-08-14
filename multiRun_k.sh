#This bash file is to run multiple simulations with respect to Kappa
#for the beam laser using cNumberCavity theory.

iFile=input.txt
nMax=8
init=21
interval=10

for ((i=0; i<nMax; i+=1)) 
do

kappa=$(($init+$i*$interval))

#$(echo "10+$i*10" | bc -l);

printf "dt 0.01
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
kappa 1.0
invT2 0
controlType k
name dt0.02_dZ0_dPz0_tau1.0__nBin20_dens100_g3_k${kappa}" > $iFile

./beamLaser -f $iFile

number=$((1+$i))

echo "Run ${number} of" $nMax

done
