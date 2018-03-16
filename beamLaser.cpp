//This program is used to simulate the beam laser 
//  --using the cNumber method 
//  --with the cavity variables
//  --using the individual variables.
#include "beamLaser.hpp"
#include "config.hpp"

//Changes required subject to the definition of Param 
void getParam(const char* filename, Param *param) 
{
  std::ifstream configInput(filename);
  std::string dummy;

  while (!configInput.eof()) {
    configInput >> dummy;
    if (configInput.eof()) break;
    if (!configInput.good()) {
      std::cout << "Bad read in input file" << std::endl;
      exit(-1);
    }
    if (dummy.compare("dt") == 0)
      configInput >> param->dt;
    else if (dummy.compare("tmax") == 0)
      configInput >> param->tmax;
    else if (dummy.compare("nstore") == 0)
      configInput >> param->nstore;
    else if (dummy.compare("nTrajectory") == 0)
      configInput >> param->nTrajectory;
    else if (dummy.compare("yWall") == 0)
      configInput >> param->yWall;
    else if (dummy.compare("sigmaXX") == 0)
      configInput >> param->sigmaX[0];
    else if (dummy.compare("sigmaXZ") == 0)
      configInput >> param->sigmaX[1];
    else if (dummy.compare("transitTime") == 0)
      configInput >> param->transitTime;
    else if (dummy.compare("sigmaPX") == 0)
      configInput >> param->sigmaP[0];
    else if (dummy.compare("sigmaPY") == 0)
      configInput >> param->sigmaP[1];
    else if (dummy.compare("sigmaPZ") == 0)
      configInput >> param->sigmaP[2];
    else if (dummy.compare("density") == 0)
      configInput >> param->density;
    else if (dummy.compare("rabi") == 0)
      configInput >> param->rabi;
    else if (dummy.compare("kappa") == 0)
      configInput >> param->kappa;
    else if (dummy.compare("name") == 0)
      configInput >> param->name;
    else {
      std::cout << "Error: invalid label " << dummy << " in "
          << filename << std::endl;
      exit(-1);
    }
  }
}

void generateInitialField(Ensemble& ensemble, const Param& param) {
  int nTimeStep = param.tmax/param.dt+0.5;
  ensemble.cavity.q.setZero(param.nTrajectory, nTimeStep+1);
  ensemble.cavity.p.setZero(param.nTrajectory, nTimeStep+1);
}

void generateExternalState(Atom& newAtom, const Param& param, const double meanP)
{
  //Initial values
  Vector3d X (0, -param.yWall, 0);
  Vector3d P (0, meanP, 0);
/*  Vector3d X (
   rng.get_gaussian_rn (param.sigmaX [0]),
   param.yWallCreate,
    rng.get_gaussian_rn (param.sigmaX [1]));
  //Only keep atoms moving to the positive y direction.
*/ 

/*  Vector3d P (
    rng.get_gaussian_rn (param.sigmaP [0]),
    rng.get_gaussian_rn (param.sigmaP [1]) + meanP,
    rng.get_gaussian_rn (param.sigmaP [2]));
  while (P[1] <= 0)
    P[1] = rng.get_gaussian_rn (param.sigmaP [1]) + meanP;
*/
  //Complete initiation
  External newExternal = {X,P};
  newAtom.external = newExternal;
}

void generateInternalState(Atom& newAtom, const Param& param)
{ 
  //For convenience
  const int nTrajectory = param.nTrajectory;

  //Set up empty internal state vectors
  VectorXd newSx, newSy, newSz;
  newSx = VectorXd::Zero(nTrajectory);
   for (int j = 0; j < nTrajectory; j++) 
       newSx[j] = double(rng.get_binomial_int(0.5, 1))*2-1; //50percent giving 1 or -1
  newSy = VectorXd::Zero(nTrajectory);
   for (int j = 0; j < nTrajectory; j++) 
       newSy[j] = double(rng.get_binomial_int(0.5, 1))*2-1; //50percent giving 1 or -1
  newSz = VectorXd::Zero(nTrajectory);

  //Initial values
  newSz.fill(1); 

  //Complete initiation
  Internal newInternal = {newSx, newSy, newSz};
  newAtom.internal = newInternal;
}

void addAtomsFromSource(Ensemble& ensemble, const Param& param, const double meanP, int& m)
{
  unsigned long int nAtom;
  const double dN = param.density*param.dt;

  if (dN >= 1) {
    nAtom = dN;//rng.get_poissonian_int(dN);      

    for (unsigned long int n = 0; n < nAtom; n++) {
      Atom newAtom; //Create a new atom
      generateExternalState(newAtom, param, meanP);    //For each atom, generate its own x and p;
      generateInternalState(newAtom, param);           //For each atom, generate its sx, sy, and sz vectors
      ensemble.atoms.push_back(newAtom);
    }
  }
  else if (dN > 0){
    int rep = 1/dN;
    if (m == 1) {
      Atom newAtom; //Create a new atom
      generateExternalState(newAtom, param, meanP);    //For each atom, generate its own x and p;
      generateInternalState(newAtom, param);           //For each atom, generate its sx, sy, and sz vectors
      ensemble.atoms.push_back(newAtom);
    }
    m += 1;
    if (m == rep) {
      m = 0;
    }
  }
  else {
    std::cout << "Bad dN in input file" << std::endl;
    exit(-1);
  }
}
  
 
void removeAtomsAtWalls(Ensemble& ensemble, const Param& param) 
{
  std::vector<Atom> newAtoms;
  for (std::vector<Atom>::iterator a = ensemble.atoms.begin(); a != ensemble.atoms.end(); a++)
    if (a->external.X[1] < param.yWall)   
      newAtoms.push_back(*a);
  ensemble.atoms = newAtoms;
}

void advanceExternalStateOneTimeStep(Ensemble& ensemble, const Param& param) 
{
  for (std::vector<Atom>::iterator a = ensemble.atoms.begin(); a != ensemble.atoms.end(); a++)
    a->external.X += param.dt * a->external.P;
}

void getDriftVector(const VectorXd& sAtoms, VectorXd& drift, const Param& param) 
{
  //For convenience
  const double rabi = param.rabi;  
  const double kappa = param.kappa;
  const int size = sAtoms.size();
  const int nAtom = (size-2)/NVAR;

  //Definition of Jx and Jy
  double jx = 0, jy = 0;
  for (int j = 0; j < nAtom; j++) {
    jx += sAtoms[NVAR*j];
    jy += sAtoms[NVAR*j+1];
  }

  //Drift vector terms. Dimension 3*nAtom, structure {D1x,D1y,D1z, D2x, D2y, D2z,...., q, p}
  drift = VectorXd::Zero(size);
  for (int j = 0; j < nAtom; j++) {
    drift[NVAR*j] = rabi/2*sAtoms[NVAR*nAtom+1]*sAtoms[NVAR*j+2];
    drift[NVAR*j+1] = -rabi/2*sAtoms[NVAR*nAtom]*sAtoms[NVAR*j+2];
    drift[NVAR*j+2] = rabi/2*(sAtoms[NVAR*nAtom]*sAtoms[NVAR*j+1]
                      -sAtoms[NVAR*nAtom+1]*sAtoms[NVAR*j]);
  }
  drift[NVAR*nAtom] = -rabi/2*jy-kappa/2*sAtoms[NVAR*nAtom];
  drift[NVAR*nAtom+1] = rabi/2*jx-kappa/2*sAtoms[NVAR*nAtom+1];
}

void advanceInternalStateOneTimeStep(Ensemble& ensemble, const Param& param, const int nStep)
{
  //For convenience
  const double dt = param.dt;
  const double rabi = param.rabi;
  const double kappa = param.kappa;
  const int nTrajectory = param.nTrajectory;
  const int nAtom = ensemble.atoms.size();
  const int size = NVAR*nAtom + 2; //3N+2

  //Loop over all trajectories. "n" stands for the number of the current trajectory.
  for (int n = 0; n < nTrajectory; n++) {

    //Y_n. Notations see "Beam laser/11.1" on ipad pro.
    VectorXd sAtoms = VectorXd::Zero(size);;//a vector of spins for all the atoms
    for (int i = 0; i < nAtom; i++) {
      sAtoms[NVAR*i] = ensemble.atoms[i].internal.sx[n];
      sAtoms[NVAR*i+1] = ensemble.atoms[i].internal.sy[n];
      sAtoms[NVAR*i+2] = ensemble.atoms[i].internal.sz[n];
    }
    sAtoms[NVAR*nAtom] = ensemble.cavity.q(n,nStep);
    sAtoms[NVAR*nAtom+1] = ensemble.cavity.p(n,nStep);
 
    //drift as a function of sAtoms
    VectorXd drift = VectorXd::Zero(size);
    getDriftVector(sAtoms, drift, param);

    //dW
    //double dw = rng.get_gaussian_rn(sqrt(dt));
    //VectorXd dW = VectorXd::Ones(size)*dw;
    VectorXd dW = VectorXd::Zero(size);
    dW[NVAR*nAtom] = sqrt(kappa)*rng.get_gaussian_rn(sqrt(dt));
    dW[NVAR*nAtom+1] = sqrt(kappa)*rng.get_gaussian_rn(sqrt(dt));
  
    //Y'. Notations see "Beam laser/11.1" on ipad pro.
    VectorXd sAtomsTemp; //Temporary value for sAtoms
    sAtomsTemp = sAtoms+drift*dt+dW;

    //driftTemp as a function of sAtomsTemp
    VectorXd driftTemp = VectorXd::Zero(size);
    getDriftVector(sAtomsTemp, driftTemp, param);

    //Y_{n+1}
    sAtoms += (driftTemp+drift)*dt/2+dW;

    //Put back
    for (int i = 0; i < nAtom; i++) {
      ensemble.atoms[i].internal.sx[n] = sAtoms[NVAR*i];
      ensemble.atoms[i].internal.sy[n] = sAtoms[NVAR*i+1];
      ensemble.atoms[i].internal.sz[n] = sAtoms[NVAR*i+2];
    }
    ensemble.cavity.q(n,nStep+1) = sAtoms[NVAR*nAtom];
    ensemble.cavity.p(n,nStep+1) = sAtoms[NVAR*nAtom+1];
  }
}

void advanceAtomsOneTimeStep(Ensemble& ensemble, const Param& param, const int nStep)
{
  advanceExternalStateOneTimeStep(ensemble, param);
  advanceInternalStateOneTimeStep(ensemble, param, nStep); //Including both atoms and cavity
}

void advanceInterval(Ensemble& ensemble, const Param& param, 
                  const double meanP, const int nStep, int& m)
{
  addAtomsFromSource(ensemble, param, meanP, m);
  removeAtomsAtWalls(ensemble, param);
  advanceAtomsOneTimeStep(ensemble, param, nStep);
}

void storeObservables(Observables& observables, int s, Ensemble& ensemble, 
    const Param& param, int nStep)
{
  //For convenience
  const double kappa = param.kappa;
  const int nTrajectory = param.nTrajectory;
  const int nAtom = ensemble.atoms.size();
  const int nTimeStep = param.tmax/param.dt+0.5;
  
  //nAtom
  observables.nAtom(s) = ensemble.atoms.size();
  
  //inversion
  double inversion = 0;
  for (int i = 0; i < nAtom; i++) {
    inversion += ensemble.atoms[i].internal.sz.sum();
  }
  observables.inversionAve(s) = inversion/nAtom/nTrajectory;

  //intensity
  observables.intensity(s) = kappa/4*(ensemble.cavity.q.col(nStep).array().square().sum()/nTrajectory                  
                                      +ensemble.cavity.p.col(nStep).array().square().sum()/nTrajectory
                                      -2);
}

void evolve(Ensemble& ensemble, const Param& param, Observables& observables)
{
  //meanP
  double meanP = param.yWall*2/param.transitTime; //vy = deltay/tau
  
  //evolve
  int nTimeStep = param.tmax/param.dt+0.5;
  double tstep = param.dt, t=0;
  int m = 1; // for atom generation
  
  //For "nTimeStep" number of data, keep "nstore" of them. 
  for (int n = 0, s = 0; n <= nTimeStep; n++, t += tstep) {
    if ((long)(n+1)*param.nstore/(nTimeStep+1) > s) {
      storeObservables(observables, s++, ensemble, param, n);
      //debug
      std::cout << "This is timestep " << n << "/" << nTimeStep << std::endl << std::endl;
      //debug
    }
    if (n != nTimeStep)
      advanceInterval(ensemble, param, meanP, n, m);
  }
}

void writeObservables(ObservableFiles& observableFiles, 
    Observables& observables, Ensemble& ensemble)
{
  observableFiles.nAtom << observables.nAtom << std::endl;
  observableFiles.intensity << observables.intensity << std::endl;
  observableFiles.inversionAve << observables.inversionAve << std::endl;
  observableFiles.qMatrix << ensemble.cavity.q << std::endl;
  observableFiles.pMatrix << ensemble.cavity.p << std::endl;
}

int main(int argc, char *argv[])
{
  //Count time
  clock_t t1,t2;
  t1=clock();
/////////////////////////////////////////////////////////////////////////////

  //Configuration. Calling functions from "config.hpp".
  CmdLineArgs config;
  getOptions(argc, argv, &config);

  //Set up parameters
  Param param;
  getParam (config.configFile, &param);
	
  //Set up initial conditions
  Ensemble ensemble;
  generateInitialField(ensemble, param);
  Observables observables(param.nstore);

  //Start simulation
  evolve(ensemble, param, observables);

  //Write Observables
  ObservableFiles observableFiles;
  writeObservables(observableFiles, observables, ensemble);
  
  //Move .dat files into the directory named "name". Calling from "config.hpp".
  mkdir(param);
  
///////////////////////////////////////////////////////////////////////////////
  //Count time
  t2=clock();
  float diff = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
  std::cout << "\nThis program takes " << diff << " seconds." << std::endl << std::endl;

  return 0;
}