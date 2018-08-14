#This bash file is to run multiple simulations with respect to nAtom
#for the beam laser using cNumberCavity theory.

iFile=input.txt
nMax=1
init=10
interval=10

for ((i=0; i<nMax; i+=1)) 
do

nAtom=$(($init+$i*$interval))

#$(echo "10+$i*10" | bc -l);

printf "dt 0.01
tmax 20
nTrajectory 1000
nStore 100
yWall 10
sigmaXX 0
sigmaXZ 0
transitTime 1
sigmaPX 0
sigmaPY 0
sigmaPZ 0
density $nAtom
rabi 1.6
kappa 25
name N${density}_transitTime1" > $iFile

./superradiantLaser -f $iFile

number=$((1+$i))

echo "Run ${number} of" $nMax

done
