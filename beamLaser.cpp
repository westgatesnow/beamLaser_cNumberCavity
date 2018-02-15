//This program is used to simulate the beam laser using the cNumber Langevin method without eliminating
//the cavity mode.
#include "beamLaser.hpp"

//Routine
void getOptions(int argc, char** argv, CmdLineArgs* cmdLineArgs)
{
  cmdLineArgs->configFile="sampleSimulation.txt";
  while (1) {
    int c;
    static struct option long_options[] = {
      {"help", no_argument, 0, 'h'},
      {"file", required_argument, 0, 'f'},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    c = getopt_long(argc, argv, "hf:", long_options, &option_index);
    if (c == -1) break;
    switch (c) {
      case 'h': std::cout << usageHeader << usageMessage;
        exit(0);
      case 'f': cmdLineArgs->configFile = optarg;
        break;
      default: exit(-1);
    }
  }
  if (optind < argc) {
    std::cout << "Error: non-option arguments: ";
    while (optind < argc) std::cout << argv[optind++] << " ";
    std::cout << std::endl;
    exit(-1);
  }
  std::cout << "Using parameters file " << cmdLineArgs->configFile << std::endl;
}

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

void addAtomsFromSource(Ensemble& ensemble, const Param& param, const double meanP)
{
  unsigned long int nAtom;

/////////////////////////////////////
// No Beam Noise. N0 changes.
//   double mean = param.density*param.dt;
//   nAtom = mean;

/////////////////////////////////////
// No Beam Noise. N0 is const.
 // nAtom = param.meanAtomGeneratingNumber;

/////////////////////////////////////
// Beam Noise. N0 is const.
  const double dN = param.density*param.dt;
  nAtom = dN;//rng.get_poissonian_int(dN);      

  for (unsigned long int n = 0; n < nAtom; n++) {
    Atom newAtom; //Create a new atom
    generateExternalState(newAtom, param, meanP);    //For each atom, generate its own x and p;
    generateInternalState(newAtom, param);           //For each atom, generate its sx, sy, and sz vectors
    ensemble.atoms.push_back(newAtom);
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

  //Drift vector terms. Dimension 3*nAtom, structure {D1+,D1-,D1z, D2+, D2-, D2z,...., q, p}
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

void advanceInternalStateOneTimeStep(Ensemble& ensemble, const Param& param)
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
    sAtoms[NVAR*nAtom] = ensemble.cavity.q[n];
    sAtoms[NVAR*nAtom+1] = ensemble.cavity.p[n];
    
    //debug
    //std::cout << "here" << std::endl << sAtoms << std::endl << std::endl;
    //debug

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
    ensemble.cavity.q[n] = sAtoms[NVAR*nAtom];
    ensemble.cavity.p[n] = sAtoms[NVAR*nAtom+1];
  }
}

void advanceAtomsOneTimeStep(Ensemble& ensemble, const Param& param)
{
  advanceExternalStateOneTimeStep(ensemble, param);
  advanceInternalStateOneTimeStep(ensemble, param); //Including both atoms and cavity
}

void advanceInterval(Ensemble& ensemble, const Param& param, 
                  const double meanP)
{
  addAtomsFromSource(ensemble, param, meanP);
  removeAtomsAtWalls(ensemble, param);
  advanceAtomsOneTimeStep(ensemble, param);
}

void storeObservables(Observables& observables, int s, Ensemble& ensemble, 
    const Param& param)
{
  //For convenience
  const double kappa = param.kappa;
  const int nTrajectory = param.nTrajectory;
  const int nAtom = ensemble.atoms.size();
  
  //nAtom
  observables.nAtom(s) = ensemble.atoms.size();
  
  //inversion
  double inversion = 0;
  for (int i = 0; i < nAtom; i++) {
    inversion += ensemble.atoms[i].internal.sz.sum();
  }
  observables.inversion(s) = inversion/nAtom/nTrajectory;

  //intensity
  observables.intensity(s) = kappa/4*(ensemble.cavity.q.array().square().sum()/nTrajectory                  
                                      +ensemble.cavity.p.array().square().sum()/nTrajectory
                                      -2);
}

void evolve(Ensemble& ensemble, const Param& param, Observables& observables)
{
  //meanP
  double meanP = param.yWall*2/param.transitTime; //vy = deltay/tau

  //Integration conditions
  double dt = param.dt;//dt = dN/N0*tau
  double tmax = param.tmax;
  
  //evolve
  int nTimeStep = tmax/dt+0.5;
  double tstep = dt, t=0;

  for (int n = 0, s = 0; n <= nTimeStep; n++, t += tstep) {
    if ((long)(n+1)*param.nstore/(nTimeStep+1) > s) {
      storeObservables(observables, s++, ensemble, param);
      //debug
      std::cout << "This is timestep " << n << "/" << nTimeStep << std::endl << std::endl;
      //debug
    }
    if (n != nTimeStep)
      advanceInterval(ensemble, param, meanP);
  }
}

void writeObservables(ObservableFiles& observableFiles, 
    Observables& observables)
{
  observableFiles.nAtom << observables.nAtom << std::endl;
  observableFiles.intensity << observables.intensity << std::endl;
  observableFiles.inversion << observables.inversion << std::endl;
}

void mkdir(Param& param) {
  std::string mkdir = "mkdir "+param.name; //make a new directory to store data
  system(mkdir.c_str());
  std::string cpInput = "cp input.txt "+param.name;
  system(cpInput.c_str());  
  std::string moveparam = "mv *.dat "+param.name;
  system(moveparam.c_str());
};


int main(int argc, char *argv[])
{
  //Count time
  clock_t t1,t2;
  t1=clock();

  //Configuration
  CmdLineArgs config;
  getOptions(argc, argv, &config);

  //Set up parameters
  Param param;
  getParam (config.configFile, &param);
	
  //Set up initial conditions
  Ensemble ensemble;
  ensemble.cavity.q = VectorXd::Zero(param.nTrajectory);
  ensemble.cavity.p = VectorXd::Zero(param.nTrajectory);
  Observables observables(param.nstore);

  //Start simulation
  evolve(ensemble, param, observables);

  //Write Observables
  ObservableFiles observableFiles;
  writeObservables(observableFiles, observables);
  
  //Move .dat files into the directory named "name"
  mkdir(param);
  
  //Count time
  t2=clock();
  float diff = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
  std::cout << "\nThis program takes " << diff << " seconds." << std::endl << std::endl;

  return 0;
}


//debug
//std::cout << ensemble.cov << std::endl << std::endl;
//exit(-1);
//debug