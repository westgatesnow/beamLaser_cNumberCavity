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
    else if (dummy.compare("gammac") == 0)
      configInput >> param->gammac;
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

void addAtomsFromSource(Ensemble& ensemble, const Param& param, 
                      const double meanP)
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
  nAtom = rng.get_poissonian_int(dN);      

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


void getDiffusionMatrix(const Ensemble& ensemble, MatrixXd& diffusion, const Param& param)
{ 
  //For convenience
  const double gc = param.gammac;
  const int nAtom = ensemble.atoms.size();
  const int nTrajectory = param.nTrajectory;

  //Diagonal block matrices
  for (int i = 0; i < nAtom; i++) {
    diffusion(NVAR*i, NVAR*i) = gc;
    diffusion(NVAR*i, NVAR*i+1) = 0;
    diffusion(NVAR*i, NVAR*i+2) = gc*ensemble.atoms[i].internal.sx.sum()/nTrajectory;
    diffusion(NVAR*i+1, NVAR*i+1) = gc;
    diffusion(NVAR*i+1, NVAR*i+2) = gc*ensemble.atoms[i].internal.sy.sum()/nTrajectory;
    diffusion(NVAR*i+2, NVAR*i+2) = 2*gc*(ensemble.atoms[i].internal.sz.sum()/nTrajectory+1);
  }

  //Off-diagonal block matrices
  for (int i = 0; i < nAtom; i++) {
    for (int j = i+1; j < nAtom; j++) {
      diffusion(NVAR*i, NVAR*j) = gc/nTrajectory
        *(ensemble.atoms[i].internal.sz.cwiseProduct(ensemble.atoms[j].internal.sz)).sum();
      diffusion(NVAR*i, NVAR*j+1) = 0;
      diffusion(NVAR*i, NVAR*j+2) = -gc/nTrajectory
        *(ensemble.atoms[i].internal.sz.cwiseProduct(ensemble.atoms[j].internal.sx)).sum();
      diffusion(NVAR*i+1, NVAR*j) = 0;
      diffusion(NVAR*i+1, NVAR*j+1) = gc/nTrajectory
        *(ensemble.atoms[i].internal.sz.cwiseProduct(ensemble.atoms[j].internal.sz)).sum();
      diffusion(NVAR*i+1, NVAR*j+2) = -gc/nTrajectory
        *(ensemble.atoms[i].internal.sz.cwiseProduct(ensemble.atoms[j].internal.sy)).sum();
      diffusion(NVAR*i+2, NVAR*j) = -gc/nTrajectory
        *(ensemble.atoms[j].internal.sz.cwiseProduct(ensemble.atoms[i].internal.sx)).sum();
      diffusion(NVAR*i+2, NVAR*j+1) = -gc/nTrajectory
        *(ensemble.atoms[j].internal.sz.cwiseProduct(ensemble.atoms[i].internal.sy)).sum();
      diffusion(NVAR*i+2, NVAR*j+2) = gc/nTrajectory
        *((ensemble.atoms[i].internal.sx.cwiseProduct(ensemble.atoms[j].internal.sx)).sum()
         +(ensemble.atoms[i].internal.sy.cwiseProduct(ensemble.atoms[j].internal.sy)).sum());
    }
  }

  //The lower half of the symmetric matrix
  for (int i = 1; i < NVAR*nAtom; i++) {
    for (int j = 0; j < i; j++) {
      diffusion(i,j) = diffusion(j,i);
    }
  }
}

void getSqrtMatrix(const MatrixXd& diffusion, MatrixXd& bMatrix)
{
  if (diffusion.size() == 0) {
    return;
  }
  else {
    SelfAdjointEigenSolver<MatrixXd> es(diffusion);
    //MatrixXd bMatrix = es.operatorSqrt();
    MatrixXd v = es.eigenvectors();
    MatrixXd vt = v.transpose();
    VectorXd eigen = es.eigenvalues();
    for (int i=0; i<eigen.size(); i++)
      if(eigen[i] < 1.0E-10) {
        eigen[i] = 0;
      //          std::cout << "set zero \n" << std::endl;
      }

    MatrixXd l = eigen.asDiagonal();
    bMatrix = v*l.cwiseSqrt()*vt;
    //Quit if NaN
    //debug
    //std::cout << "mMatrix is " << std::endl << diffusion << std::endl << std::endl;
    std::cout << "Eigenvalues of mMatrix is " << std::endl << es.eigenvalues() << std::endl << std::endl;
    //std::cout << "bMatrix is " << std::endl << bMatrix << std::endl << std::endl;
    //debug
    if (bMatrix(1,1) != bMatrix(1,1)) {
      std::cout << "Not a number\n" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  
}

void getDriftVector(const VectorXd& sAtoms, VectorXd& drift, const Param& param) 
{
  //For convenience
  const double gc = param.gammac;
  const int size = sAtoms.size();
  const int nAtom = size/NVAR;
  
  //Definition of Jx and Jy
  double jx = 0, jy = 0;
  for (int j = 0; j < nAtom; j++) {
    jx += sAtoms[NVAR*j];
    jy += sAtoms[NVAR*j+1];
  }

  //Drift vector terms. Dimension 3*nAtom, structure {D1+,D1-,D1z, D2+, D2-, D2z,....}
  drift = VectorXd::Zero(size);
  for (int j = 0; j < nAtom; j++) {
    drift[NVAR*j] = gc/2*(jx*sAtoms[NVAR*j+2]-sAtoms[NVAR*j]);
    drift[NVAR*j+1] = gc/2*(jy*sAtoms[NVAR*j+2]-sAtoms[NVAR*j+1]);
    drift[NVAR*j+2] = -gc/2*(jx*sAtoms[NVAR*j]+jy*sAtoms[NVAR*j+1]-pow(sAtoms[NVAR*j],2)-pow(sAtoms[NVAR*j+1],2)
                     +2+2*sAtoms[NVAR*j+2]);
  }
  
}

void advanceInternalStateOneTimeStep(Ensemble& ensemble, const Param& param)
{
  //For convenience
  const double dt = param.dt;
  const double gc = param.gammac;
  const int nTrajectory = param.nTrajectory;
  const int nAtom = ensemble.atoms.size();
  const int size = NVAR*nAtom;

  //Diffusion matrix 2*m_ij. The data structure is {1+,1-,1z,2+,2-,2z,...}^2;
  MatrixXd diffusion = MatrixXd::Zero(size, size);
  getDiffusionMatrix(ensemble, diffusion, param);
  //Get the sqrt matrix.
  MatrixXd bMatrix = MatrixXd::Zero(size, size);
  getSqrtMatrix(diffusion, bMatrix);
  //Loop over all trajectories. "n" stands for the number of the current trajectory.
  for (int n = 0; n < nTrajectory; n++) {

    //Y_n. Notations see "Beam laser/11.1" on ipad pro.
    VectorXd sAtoms = VectorXd::Zero(size);;//a vector of spins for all the atoms
    for (int i = 0; i < nAtom; i++) {
      sAtoms[NVAR*i] = ensemble.atoms[i].internal.sx[n];
      sAtoms[NVAR*i+1] = ensemble.atoms[i].internal.sy[n];
      sAtoms[NVAR*i+2] = ensemble.atoms[i].internal.sz[n];
    }

    //drift as a function of sAtoms
    VectorXd drift = VectorXd::Zero(size);
    getDriftVector(sAtoms, drift, param);

    //dW
    //double dw = rng.get_gaussian_rn(sqrt(dt));
    //VectorXd dW = VectorXd::Ones(size)*dw;

    VectorXd dW = VectorXd::Zero(size);
    for (int i = 0; i < size; i++) {
      dW[i] = rng.get_gaussian_rn(sqrt(dt));
    }
    
    //Y'. Notations see "Beam laser/11.1" on ipad pro.
    VectorXd sAtomsTemp; //Temporary value for sAtoms
    sAtomsTemp = sAtoms+drift*dt+bMatrix*dW;

    //driftTemp as a function of sAtomsTemp
    VectorXd driftTemp = VectorXd::Zero(size);
    getDriftVector(sAtomsTemp, driftTemp, param);

    //Y_{n+1}
    sAtoms += (driftTemp+drift)*dt/2+bMatrix*dW;

    //Put back
    for (int i = 0; i < nAtom; i++) {
      ensemble.atoms[i].internal.sx[n] = sAtoms[NVAR*i];
      ensemble.atoms[i].internal.sy[n] = sAtoms[NVAR*i+1];
      ensemble.atoms[i].internal.sz[n] = sAtoms[NVAR*i+2];
    }
  }
}

void advanceAtomsOneTimeStep(Ensemble& ensemble, const Param& param)
{
  advanceExternalStateOneTimeStep(ensemble, param);
  advanceInternalStateOneTimeStep(ensemble, param);
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
  const double gc = param.gammac;
  const int nTrajectory = param.nTrajectory;
  const int nAtom = ensemble.atoms.size();
  
  //The order of the definitions of observables matters. 

  //nAtom
  observables.nAtom(s) = ensemble.atoms.size();
  
  //inversion
  double inversion = 0;
  for (int i = 0; i < nAtom; i++) {
    inversion += ensemble.atoms[i].internal.sz.sum();
  }
  observables.inversion(s) = inversion/nAtom/nTrajectory;

  //intensityUnCor
  double intensityUnCor = 0;
  observables.intensityUnCor(s) = gc/2*(nAtom+inversion/nTrajectory);

  //spinSpinCor
  double spinSpinCor =0;
  for (int n = 0; n < nTrajectory; n++) {
    for (int i = 0; i < nAtom; i++) {
      for (int j = 0; j < nAtom; j++) {
        spinSpinCor += (ensemble.atoms[i].internal.sx[n]*ensemble.atoms[j].internal.sx[n]
                       +ensemble.atoms[i].internal.sy[n]*ensemble.atoms[j].internal.sy[n])/4;
      }
    }
  }
  observables.spinSpinCor(s) = spinSpinCor/nAtom/(nAtom-1)/nTrajectory;

  //intensity
  observables.intensity(s) = gc*spinSpinCor/nTrajectory+observables.intensityUnCor(s);
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
  observableFiles.intensityUnCor << observables.intensityUnCor << std::endl;
  observableFiles.inversion << observables.inversion << std::endl;
  observableFiles.spinSpinCor << observables.spinSpinCor << std::endl;
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
//debug