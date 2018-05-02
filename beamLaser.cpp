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

void generateInitialField(Ensemble& ensemble, const Param& param) 
{
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
  // newSx.fill(0);
  // newSy.fill(0);
  newSz.fill(1); 

  //Complete initiation
  Internal newInternal = {newSx, newSy, newSz};
  newAtom.internal = newInternal;
}

void addAtomsFromSource(Ensemble& ensemble, const Param& param, const double meanP, int& m)
{
  unsigned long int nAtom;
  const double dN = param.density*param.dt;

  ///////////////////////////////////////////////////////////////////////////
  // Uniform atom generation
  ///////////////////////////////////////////////////////////////////////////
  if (dN >= 1) {
    nAtom = dN;     

    for (unsigned long int n = 0; n < nAtom; n++) {
      Atom newAtom; //Create a new atom
      generateExternalState(newAtom, param, meanP);    //For each atom, generate its own x and p;
      generateInternalState(newAtom, param);           //For each atom, generate its sx, sy, and sz vectors
      ensemble.atoms.push_back(newAtom);
    }
  }

  else if (dN > 0){
    int rep = 1/dN;
    if (m == rep) {
      Atom newAtom; //Create a new atom
      generateExternalState(newAtom, param, meanP);    //For each atom, generate its own x and p;
      generateInternalState(newAtom, param);           //For each atom, generate its sx, sy, and sz vectors
      ensemble.atoms.push_back(newAtom);
      m = 0;
    }
    m++;
  }
  else {
    std::cout << "Bad dN in input file" << std::endl;
    exit(-1);
  }

  ///////////////////////////////////////////////////////////////////////////
  //Poissonian atom generation
  ///////////////////////////////////////////////////////////////////////////
  // nAtom = rng.get_poissonian_int(dN);
  // for (unsigned long int n = 0; n < nAtom; n++) {
  //     Atom newAtom; //Create a new atom
  //     generateExternalState(newAtom, param, meanP);    //For each atom, generate its own x and p;
  //     generateInternalState(newAtom, param);           //For each atom, generate its sx, sy, and sz vectors
  //     ensemble.atoms.push_back(newAtom);
  // }
}
  
 
void removeAtomsAtWalls(Ensemble& ensemble, const Param& param, Vector3d& spinVar) 
{
  std::vector<Atom> newAtoms;
  int nLeaving = 0;
  for (std::vector<Atom>::iterator a = ensemble.atoms.begin(); a != ensemble.atoms.end(); a++) {
    if (a->external.X[1] < param.yWall) {
      newAtoms.push_back(*a);
    } 
    else {
      nLeaving += 1;
      spinVar(0) += a->internal.sx.sum();
      spinVar(1) += a->internal.sy.sum();
      spinVar(2) += a->internal.sz.sum();
    }
  }
  ensemble.atoms = newAtoms;
  spinVar = spinVar/param.nTrajectory/nLeaving;
}

void advanceExternalStateOneTimeStep(Ensemble& ensemble, const Param& param) 
{
  for (std::vector<Atom>::iterator a = ensemble.atoms.begin(); a != ensemble.atoms.end(); a++)
    a->external.X += param.dt * a->external.P;
}

void getDiffusionVector(VectorXd& dW, const Param& param) 
{
  //For convenience
  const double dt = param.dt;
  const double kappa = param.kappa;
  const int nAtom = (dW.size()-2)/NVAR;
  //Diffusion for the field
  dW[NVAR*nAtom] = sqrt(kappa)*rng.get_gaussian_rn(sqrt(dt));
  dW[NVAR*nAtom+1] = sqrt(kappa)*rng.get_gaussian_rn(sqrt(dt));
}

void getDriftVector(const VectorXd& sVar, VectorXd& drift, const Param& param) 
{
  //For convenience
  const double rabi = param.rabi;  
  const double kappa = param.kappa;
  const int size = sVar.size();
  const int nAtom = (size-2)/NVAR;

  //Definition of Jx and Jy
  double jx = 0, jy = 0;
  for (int j = 0; j < nAtom; j++) {
    jx += sVar[NVAR*j];
    jy += sVar[NVAR*j+1];
  }

  //Drift vector terms. Dimension 3*nAtom, structure {D1x,D1y,D1z, D2x, D2y, D2z,...., q, p}
  for (int j = 0; j < nAtom; j++) {
    drift[NVAR*j] = rabi/2*sVar[NVAR*nAtom+1]*sVar[NVAR*j+2];
    drift[NVAR*j+1] = -rabi/2*sVar[NVAR*nAtom]*sVar[NVAR*j+2];
    drift[NVAR*j+2] = rabi/2*(sVar[NVAR*nAtom]*sVar[NVAR*j+1]
                      -sVar[NVAR*nAtom+1]*sVar[NVAR*j]);
  }
  drift[NVAR*nAtom] = -rabi/2*jy-kappa/2*sVar[NVAR*nAtom];
  drift[NVAR*nAtom+1] = rabi/2*jx-kappa/2*sVar[NVAR*nAtom+1];
}

void stochasticIntegration(VectorXd& sVar, const VectorXd& drift, const VectorXd& dW, const Param& param) {
    //For convenience
    double dt = param.dt;
    int dim = drift.size();
    //Y'. Notations see "Beam laser/11.1" on ipad pro.
    VectorXd sVarTemp = sVar+drift*dt+dW;//Temporary value for sVar
    //driftTemp as a function of sVarTemp
    VectorXd driftTemp = VectorXd::Zero(dim);
    getDriftVector(sVarTemp, driftTemp, param);
    //Y_{n+1}
    sVar += (driftTemp+drift)*dt/2+dW;
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
    VectorXd sVar = VectorXd::Zero(size);;//a vector of spins for all the atoms
    for (int i = 0; i < nAtom; i++) {
      sVar[NVAR*i] = ensemble.atoms[i].internal.sx[n];
      sVar[NVAR*i+1] = ensemble.atoms[i].internal.sy[n];
      sVar[NVAR*i+2] = ensemble.atoms[i].internal.sz[n];
    }
    sVar[NVAR*nAtom] = ensemble.cavity.q(n, nStep);
    sVar[NVAR*nAtom+1] = ensemble.cavity.p(n, nStep);
 
    //drift
    VectorXd drift = VectorXd::Zero(size);
    getDriftVector(sVar, drift, param);
    //diffusion
    VectorXd dW = VectorXd::Zero(size);
    getDiffusionVector(dW, param);
    //integration
    stochasticIntegration(sVar, drift, dW, param);
    //Put back
    for (int i = 0; i < nAtom; i++) {
      ensemble.atoms[i].internal.sx[n] = sVar[NVAR*i];
      ensemble.atoms[i].internal.sy[n] = sVar[NVAR*i+1];
      ensemble.atoms[i].internal.sz[n] = sVar[NVAR*i+2];
    }
    ensemble.cavity.q(n,nStep+1) = sVar[NVAR*nAtom];
    ensemble.cavity.p(n,nStep+1) = sVar[NVAR*nAtom+1];
  }
}

void advanceAtomsOneTimeStep(Ensemble& ensemble, const Param& param, const int nStep)
{
  advanceExternalStateOneTimeStep(ensemble, param);
  advanceInternalStateOneTimeStep(ensemble, param, nStep); //Including both atoms and cavity
}

void advanceInterval(Ensemble& ensemble, const Param& param, 
                  const double meanP, const int nStep, int& m, Vector3d& spinVar)
{
  //Newly added atoms are in the tail of the "atoms" vector, so the first atoms 
  //in the "atoms" vector will be the first to be removed.
  addAtomsFromSource(ensemble, param, meanP, m);
  removeAtomsAtWalls(ensemble, param, spinVar);
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
  
  //intensity
  observables.intensity(s) = kappa/4*(ensemble.cavity.q.col(nStep).array().square().sum()/nTrajectory                  
                                      +ensemble.cavity.p.col(nStep).array().square().sum()/nTrajectory
                                      -2);
  //inversionAve
  double inversionAve = 0;
  for (int i = 0; i < nAtom; i++)
    inversionAve += ensemble.atoms[i].internal.sz.sum();
  observables.inversionAve(s) = inversionAve/nAtom/nTrajectory;
  
  //qMatrix and pMatrix
  for (int i = 0; i < nTrajectory; i++) {
    observables.qMatrix(i,s) = ensemble.cavity.q(i,nStep);
    observables.pMatrix(i,s) = ensemble.cavity.p(i,nStep);
  }

  //spinSpinCor
  double binSize = param.yWall*2/NBIN;
  VectorXd xx = VectorXd::Zero(NBIN);
  VectorXd xy = VectorXd::Zero(NBIN);
  VectorXd yx = VectorXd::Zero(NBIN);
  VectorXd yy = VectorXd::Zero(NBIN);
  VectorXi m_y = VectorXi::Zero(NBIN);
  //Find the atoms in the first bin
  VectorXd jx_0 = VectorXd::Zero(nTrajectory);
  VectorXd jy_0 = VectorXd::Zero(nTrajectory);
  for (std::vector<Atom>::iterator a = ensemble.atoms.begin(); a != ensemble.atoms.end(); a++) {
    if (a->external.X[1] < -param.yWall+binSize) {
      jx_0 += a->internal.sx;
      jy_0 += a->internal.sy;
      m_y(0) += 1;
    }
  }               //Can also start from the end of the array and then jump out of the loop 
                  //when no new atoms appear for a while.

  if(m_y(0) == 0) {
   std::cout << "\nBin size too small." << std::endl << std::endl; 
  }
  //Get spinSpinCor
  for (std::vector<Atom>::iterator a = ensemble.atoms.begin(); a != ensemble.atoms.end(); a++) {
    int binNumber = (a->external.X[1]+param.yWall)/binSize;
    if (binNumber > NBIN-1)
      binNumber = NBIN-1;
    if (binNumber != 0) {
      m_y(binNumber) += 1;
      xx(binNumber) += (jx_0.cwiseProduct(a->internal.sx)).sum()/nTrajectory;
      xy(binNumber) += (jx_0.cwiseProduct(a->internal.sy)).sum()/nTrajectory;
      yx(binNumber) += (jy_0.cwiseProduct(a->internal.sx)).sum()/nTrajectory;
      yy(binNumber) += (jy_0.cwiseProduct(a->internal.sy)).sum()/nTrajectory;
    }
  } 
  //Store spinSpinCor
  for (int i = 0; i < NBIN; i++) {
    observables.spinSpinCor_re(i,s) = 1.0/4*(xx(i)+yy(i))/m_y(0)/m_y(i);
    observables.spinSpinCor_im(i,s) = 1.0/4*(yx(i)-xy(i))/m_y(0)/m_y(i);
  }
}

void storeSpinVariables(SpinVariables& spinVariables, const Ensemble& ensemble, 
  const Param& param, int n, Vector3d& spinVar)
{
  spinVariables.sxFinal(n) = spinVar(0);
  spinVariables.syFinal(n) = spinVar(1);
  spinVariables.szFinal(n) = spinVar(2);
}


void evolve(Ensemble& ensemble, const Param& param, Observables& observables, SpinVariables& spinVariables)
{
  //meanP
  double meanP = param.yWall*2/param.transitTime; //vy = deltay/tau
  //evolve
  int nTimeStep = param.tmax/param.dt+0.5;
  double t = 0;
  //for atom generation
  int m = 0; 
  //for spinVariables
  Vector3d spinVar;    

  //For "nTimeStep" number of data, keep "nstore" of them. 
  for (int n = 0, s = 0; n <= nTimeStep; n++, t += param.dt) {
    if ((long)(n+1)*param.nstore/(nTimeStep+1) > s) {
      storeObservables(observables, s++, ensemble, param, n);
      std::cout << "Data " << s << "/" << param.nstore << " stored." << std::endl << std::endl;
    }
    if (n != nTimeStep) {
      spinVar << 0, 0, 0;
      advanceInterval(ensemble, param, meanP, n, m, spinVar);
      //store all the spinVariables; since when dN < 1, there are many NaN's
      storeSpinVariables(spinVariables, ensemble, param, n, spinVar);
    }
  }
}

void writeObservables(ObservableFiles& observableFiles, 
    Observables& observables, SpinVariables& spinVariables, Ensemble& ensemble)
{
  std::cout << "Writing data... (This may take several minutes.)" << std::endl << std::endl;
  observableFiles.nAtom << observables.nAtom << std::endl;
  observableFiles.intensity << observables.intensity << std::endl;
  observableFiles.inversionAve << observables.inversionAve << std::endl;
  observableFiles.qMatrix << observables.qMatrix << std::endl;
  observableFiles.pMatrix << observables.pMatrix << std::endl;
  observableFiles.spinSpinCor_re << observables.spinSpinCor_re << std::endl;
  observableFiles.spinSpinCor_im << observables.spinSpinCor_im << std::endl;
  observableFiles.sxFinal << spinVariables.sxFinal << std::endl;
  observableFiles.syFinal << spinVariables.syFinal << std::endl;
  observableFiles.szFinal << spinVariables.szFinal << std::endl;  
}

void mkdir(Param& param) 
{
  std::string mkdir = "mkdir "+param.name; //make a new directory to store data
  system(mkdir.c_str());
  std::string cpInput = "cp input.txt "+param.name;
  system(cpInput.c_str());  
  std::string moveparam = "mv *.dat "+param.name;
  system(moveparam.c_str());
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
  int nTimeStep = param.tmax/param.dt+0.5;

  //Set up initial conditions
  Ensemble ensemble;
  generateInitialField(ensemble, param);
  Observables observables(param.nstore, param.nTrajectory, NBIN);
  SpinVariables spinVariables(nTimeStep);

  //Start simulation
  evolve(ensemble, param, observables, spinVariables);

  //Write Observables
  ObservableFiles observableFiles;
  writeObservables(observableFiles, observables, spinVariables, ensemble);
  
  //Move .dat files into the directory named "name". 
  mkdir(param);
  
///////////////////////////////////////////////////////////////////////////////
  //Count time
  t2=clock();
  float diff = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
  std::cout << "\nThis program takes " << diff << " seconds." << std::endl << std::endl;

  return 0;
}