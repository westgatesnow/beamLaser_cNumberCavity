#ifndef __BEAMLASER__HPP__
#define __BEAMLASER__HPP__
//This program is used to simulate the beam laser 
//  --using the cNumber method 
//  --with the cavity variables
//  --using the individual variables.

//Include Eigen package
//Work in Eigen namespace
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/Eigenvalues>
using namespace Eigen;

//Include standard packages
#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>
#include <time.h>
#include <stdlib.h> 

//Include and define RNG
#include "RNG.hpp"
RNG rng(time(NULL));

//Some global constants
#define NVAR 3 // 3 variables for each atom for this code;
#define NBIN 10 //seperate the cavity length into 100 bins 

//Define data structure

//Atom external states
typedef struct 
{
  Vector3d X;     //position
  Vector3d P;     //momentum. We suppose mass is one, so momentum is velocity.
} External;

//Atom internal states
typedef struct 
{
  VectorXd sx;       //sigma_x. Dim: nTrajectory
  VectorXd sy;       //sigma_y. Dim: nTrajectory
  VectorXd sz;       //sigma_z. Dim: nTrajectory
} Internal;

//Atom total states
typedef struct 
{
  External external;  //The position and velocity of an atom.
  Internal internal;  //The internal states of an atom. We keep track of sx, sy, and sz
                        //of a single atom at all times for all trajectories.
} Atom;

//Cavity field
typedef struct 
{
  MatrixXd q;        //q = a^dagger + a. Dim: nTrajectory*(nTimeStep+1)
  MatrixXd p;        //p = - I (a^dagger - a). Dim: nTrajectory*(nTimeStep+1)
} Cavity;

//Ensemble of the system
typedef struct 
{
  std::vector<Atom> atoms;
  Cavity cavity;
} Ensemble;

//Simulation parameters
typedef struct Param 
{
  //parameter definitions
  double dt; 
  double tmax;
  int nstore; // number of times to store observables
  int nTrajectory; //number of trajectories
  //beam parameters
  double yWall; //position of the wall where atoms are destroyed
                //The coordinated are chosen s.t. atoms are created at -yWall.
                //The walls are assumed to be in xz plane.
  Vector2d sigmaX;  //standard deviation of position in xz. y deviation is 
                    //taken care of by the Poisson distribution.
  double transitTime;     //the transit time tau1>0, unit 1/gammaC.
  Vector3d sigmaP;  //standard deviation of momentum
  double density;   //the mean number of atoms per unit time;
  double rabi;    //single-atom rabi frequency
  double kappa;   //caivty decay rate. Condition of bad cavity: dt>>1/kappa.
  //Other parameters
  std::string name; //name of the directory to store results

  //Constructor; initial values
  Param() : dt(0.01), tmax(10), nstore(10), nTrajectory(1), yWall(5.0e0), 
    sigmaX(0.0e0,0.0e0), transitTime(1.0e0), sigmaP(0.0e0,0.0e0,0.0e0), 
    density(1.0), rabi(10), kappa(1000), name("abracadabra")
  {}

} Param;

std::ostream& operator<< (std::ostream& o, const Param& s)
{
  o << s.dt << std::endl;
  o << s.tmax << std::endl;
  o << s.nstore << std::endl;
  o << s.nTrajectory << std::endl;
  o << s.yWall << std::endl;
  o << s.sigmaX << std::endl;
  o << s.transitTime << std::endl;
  o << s.sigmaP << std::endl;
  o << s.density << std::endl;
  o << s.rabi << std::endl;
  o << s.kappa << std::endl; 
  return o;
}

//Observables; n is the nTimeStep
typedef struct Observables 
{
  VectorXi nAtom; 
  VectorXd intensity;
  VectorXd inversionAve;
  VectorXd spinSpinCorAve_re;
  VectorXd spinSpinCorAve_im;
  MatrixXd qMatrix;
  MatrixXd pMatrix;
  MatrixXd spinSpinCor_re;
  MatrixXd spinSpinCor_im;

  Observables(const int n, const int m, const int l) 
  {
    nAtom = VectorXi(n); 
    intensity = VectorXd(n);
    inversionAve = VectorXd(n);
    spinSpinCorAve_re = VectorXd(n);
    spinSpinCorAve_im = VectorXd(n);
    qMatrix = MatrixXd(m,n);
    pMatrix = MatrixXd(m,n);
    spinSpinCor_re = MatrixXd(l,n);
    spinSpinCor_im = MatrixXd(l,n);
  } 

} Observables;

typedef struct SpinVariables 
{
  //Definition
  VectorXd sxFinal;
  VectorXd syFinal;
  VectorXd szFinal;

  //Constructor
  SpinVariables(const int m)                                      
  { 
    sxFinal = VectorXd(m);
    syFinal = VectorXd(m);
    szFinal = VectorXd(m);
  }

} SpinVariables;

typedef struct ObservableFiles 
{
  //Definition
  std::ofstream nAtom, intensity, inversionAve, spinSpinCorAve_re, spinSpinCorAve_im,
    sxFinal, syFinal, szFinal, qMatrix, pMatrix, spinSpinCor_re, spinSpinCor_im;

  //Constructor              
  ObservableFiles() : nAtom("nAtom.dat"), 
                      intensity("intensity.dat"), 
                      inversionAve("inversionAve.dat"),
                      spinSpinCorAve_re("spinSpinCorAve_re.dat"),
                      spinSpinCorAve_im("spinSpinCorAve_im.dat"),
                      sxFinal("sxFinal.dat"),
                      syFinal("syFinal.dat"),
                      szFinal("szFinal.dat"),
                      qMatrix("qMatrix.dat"),
                      pMatrix("pMatrix.dat"),
                      spinSpinCor_re("spinSpinCor_re.dat"),
                      spinSpinCor_im("spinSpinCor_im.dat")
  {}
  
  //Deconstructor
  ~ObservableFiles() 
  {
    nAtom.close();
    intensity.close();
    inversionAve.close();
    spinSpinCorAve_re.close();
    spinSpinCorAve_im.close();
    sxFinal.close();
    syFinal.close();
    szFinal.close();
    qMatrix.close();
    pMatrix.close();
    spinSpinCor_re.close();
    spinSpinCor_im.close();
  }
  
} ObservableFiles;

#endif