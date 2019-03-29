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
#include <getopt.h>
#include <time.h>
#include <stdlib.h> 
//Define pi
#define _USE_MATH_DEFINES
#include <cmath>

//Include and define RNG
#include "RNG.hpp"
RNG rng(time(NULL));

//Some global constants
#define NVAR 3 // 3 variables for each atom for this code;
//#define NBIN 10 //seperate the cavity length into 100 bins 

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
  // MatrixXd Jx;      //Jx = \sum \cos{kz_j}s_j^x Dim: nTrajectory*(nTimeStep+1)
  MatrixXd Jy;      //Jy = \sum \cos{kz_j}s_j^y Dim: nTrajectory*(nTimeStep+1)
  MatrixXd Jz;      //Jz = \sum \cos{kz_j}s_j^z Dim: nTrajectory*(nTimeStep+1)
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
  int nStore; // number of times to store observables
  int nTrajectory; //number of trajectories
  int nBin; //number of bins along the y direction
  //beam parameters
  double yWall; //position of the wall where atoms are destroyed
                //The coordinated are chosen s.t. atoms are created at -yWall.
                //The walls are assumed to be in xz plane.
  double lambda;  //The wavelength of the laser light
  double deltaZ;  //z direction standard deviation. 
                    //y direction is taken care of by the Poisson distribution. 
                    //x direction is irrelevant, but we still keep it here.
  double deltaPz;  //standard deviation of z velocity
  double transitTime;     //the transit time tau1>0, unit 1/gammaC.
  double density;   //the mean number of atoms per unit time;
  double rabi;    //single-atom rabi frequency
  double kappa;   //caivty decay rate. Condition of bad cavity: dt>>1/kappa.
  //Other parameters
  double invT2;  //The T2 dephasing time inverse
  std::string controlType; //name of the parent directory for a certain controlle variable
  std::string name; //name of the directory to store results
  int pois; //If the Poissonian noise is included; 0 = no, 1 = yes

  //Constructor; initial values
  Param() : dt(0.01), tmax(10), 
            nStore(10), nTrajectory(1), nBin(10), 
            yWall(10.0), lambda(1.0),
            deltaZ(0.0), deltaPz(0.0), 
            transitTime(1.0), density(1.0), 
            rabi(10.0), kappa(1000.0), invT2(0),
            controlType("test"), name("aProgramHasNoName"), pois(0),
  {}

} Param;

std::ostream& operator<< (std::ostream& o, const Param& s)
{
  o << s.dt << std::endl;
  o << s.tmax << std::endl;
  o << s.nStore << std::endl;
  o << s.nTrajectory << std::endl;
  o << s.nBin << std::endl;
  o << s.yWall << std::endl;
  o << s.lambda << std::endl;
  o << s.deltaZ << std::endl;
  o << s.deltaPz << std::endl;
  o << s.transitTime << std::endl;
  o << s.density << std::endl;
  o << s.rabi << std::endl;
  o << s.kappa << std::endl;
  o << s.invT2 << std::endl;    
  return o;
}

//Observables; n is the nTimeStep
typedef struct Observables 
{
  VectorXi nAtom; 
  VectorXd intensity;
  VectorXd inversionAve;
  MatrixXd qMatrix;
  MatrixXd pMatrix;
  // MatrixXd JxMatrix;
  MatrixXd JyMatrix;
  MatrixXd JzMatrix;
  // MatrixXd sxMatrix;
  // MatrixXd syMatrix;
  MatrixXd szMatrix;
  VectorXd spinSpinCorAve_re;
  // VectorXd spinSpinCorAve_im;
  MatrixXd spinSpinCor_re;
  // MatrixXd spinSpinCor_im;

  Observables(const int nStore, const int nTrajectory, const int nBin)
  {
    nAtom = VectorXi(nStore); 
    intensity = VectorXd(nStore);
    inversionAve = VectorXd(nStore);
    qMatrix = MatrixXd(nTrajectory, nStore);
    pMatrix = MatrixXd(nTrajectory, nStore);
    // JxMatrix = MatrixXd(nTrajectory, nStore);
    JyMatrix = MatrixXd(nTrajectory, nStore);
    JzMatrix = MatrixXd(nTrajectory, nStore);
    // sxMatrix = MatrixXd(nBin, nStore);
    // syMatrix = MatrixXd(nBin, nStore);
    szMatrix = MatrixXd(nBin, nStore);
    spinSpinCorAve_re = VectorXd(nStore);
    // spinSpinCorAve_im = VectorXd(nStore);
    int nBinSquare = nBin * nBin;
    spinSpinCor_re = MatrixXd(nBinSquare, nStore);
    // spinSpinCor_im = MatrixXd(nBinSquare, nStore);
  } 

} Observables;

typedef struct SpinVariables 
{
  //Definition
  //Definition for the average values of sx, sy, and sz as each atom leaves the cavity
  // VectorXd sxFinal;
  // VectorXd syFinal;
  VectorXd szFinal;

  //Constructor
  SpinVariables(const int m)                                      
  { 
    // sxFinal = VectorXd(m);
    // syFinal = VectorXd(m);
    szFinal = VectorXd(m);
  }

} SpinVariables;

typedef struct ObservableFiles 
{
  //Definition
  std::ofstream nAtom, 
                intensity, 
                inversionAve, 
                qMatrix, 
                pMatrix,
                // JxMatrix, 
                JyMatrix,
                JzMatrix,
                spinSpinCorAve_re, 
                // spinSpinCorAve_im,
                spinSpinCor_re, 
                // spinSpinCor_im,
                // sxMatrix,
                // syMatrix,
                szMatrix,
                // sxFinal, 
                // syFinal, 
                szFinal;

  //Constructor              
  ObservableFiles() : nAtom("nAtom.dat"), 
                      intensity("intensity.dat"), 
                      inversionAve("inversionAve.dat"),
                      qMatrix("qMatrix.dat"),
                      pMatrix("pMatrix.dat"),
                      // JxMatrix("JxMatrix.dat"),
                      JyMatrix("JyMatrix.dat"),
                      JzMatrix("JzMatrix.dat"),
                      spinSpinCorAve_re("spinSpinCorAve_re.dat"),
                      // spinSpinCorAve_im("spinSpinCorAve_im.dat"),
                      spinSpinCor_re("spinSpinCor_re.dat"),
                      // spinSpinCor_im("spinSpinCor_im.dat"),
                      // sxMatrix("sxMatrix.dat"),
                      // syMatrix("syMatrix.dat"),
                      szMatrix("szMatrix.dat"),
                      // sxFinal("sxFinal.dat"),
                      // syFinal("syFinal.dat"),
                      szFinal("szFinal.dat")
  {}
  
  //Deconstructor
  ~ObservableFiles() 
  {
    nAtom.close();
    intensity.close();
    inversionAve.close();
    qMatrix.close();
    pMatrix.close();
    // JxMatrix.close();
    JyMatrix.close();
    JzMatrix.close();
    spinSpinCorAve_re.close();
    // spinSpinCorAve_im.close();
    spinSpinCor_re.close();
    // spinSpinCor_im.close();
    // sxMatrix.close();
    // syMatrix.close();
    szMatrix.close();
    // sxFinal.close();
    // syFinal.close();
    szFinal.close();
  }
  
} ObservableFiles;

#endif