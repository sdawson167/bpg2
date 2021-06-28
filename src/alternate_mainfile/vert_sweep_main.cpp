#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <utility>
#include <cstring>
#include <fstream>
#include <cmath>
#include <chrono>

#include "utils.h"
#include "fftw3.h"

#include "BpgMinimizer.h"
#include "FieldProvider.h"
#include "GenericProvider.h"

#include "PwFunctionalCalculator.h"

int main(int argc, char** argv)
{

  /*
   * ===============================================================================
   *		             set lam1, lam2 values
   * ===============================================================================
   */

  // LB: 0,0
  // OK: 1,1

  const double lam1 = 0.0;
  const double lam2 = 0.0;


  /*
   * ===============================================================================
   *                      get gamma value and phaseID(s) 
   * ===============================================================================
   */ 
  
  double gamma;
  std::vector<int> phaseIDList;

  // make sure there are enough input args to extract tau and phaseID vals:
  if (argc == 1) {
    std::cout << "error: must choose gamma value and one or more phase IDs (1 - 7) to run program." << std::endl;

    return 1;
  
  } else if (argc == 2)  {
    std::cout << "error: must choose one or more phase IDs (1 - 7) to run program." << std::endl; 

    return 1;

  } else {
    gamma = std::atof(argv[1]);

    for (int i = 2; i < argc; i++) {
      int p = std::atoi(argv[i]);
      phaseIDList.push_back(p);
    }
  }

  /*
   * ==============================================================================
   *               Start optimization for all points listed in input
   * ==============================================================================
   */

  // minimization parameters (max iterations and error tolerance)
  const double errorTol = 1e-9;
  const int    maxIter  = 15;
  
  // vertical resolution: 
  const double dTau = 0.01;

  // ensure we don't print more digits of precion than we know (based on errorTol):
  std::cout << std::fixed << std::setprecision(8);

  // create minimizer (object that handles implementation of minimization algorithm) 
  BpgMinimizer<PwFunctionalCalculator> minimizer(errorTol, maxIter);

  // create vector of vectors used to store free-energies:
  typedef std::vector<double>    doubleVec;
  typedef std::vector<doubleVec> doubleVecVec;

  doubleVecVec freeEnergiesList;

  // loop through all phases
  for (std::vector<int>::iterator it = phaseIDList.begin(); 
       it != phaseIDList.end();
       it++) 
    {
      int phaseID = *it;

      // start timer
      // auto start = std::chrono::steady_clock::now();
    
      // get starting value of tau for this phase: 
      const double tau0 = -0.05;

      // make list of tau values starting at tau0:
      doubleVec tauVec;
      
      tauVec.push_back(tau0);

      double tau = tau0 + dTau;
      while (tau < 0.51) {
	tauVec.push_back(tau);

        tau += dTau;
      }

      // create field-provider object (represents the order parameter in real and fourier space) 
      // initialize field-provider in chosen phase
      GenericPhaseProvider provider(phaseID); 
      FieldProvider *field = provider.generateInitialCondition();

      // create vector to store free-energy for this phase
      doubleVec freeEnergies;

      // loop through phase points (tau, gamma, ... vals) - compute free-energy at each point
      for (std::vector<double>::const_iterator it = tauVec.begin(); it != tauVec.end(); it++)
      {
        // get gamma value
        double tau = *it; 
  
        // create 'functional calculator' object - handles free-energy calculations
        PwFunctionalCalculator calculator(tau, gamma, lam1, lam2);
  
        // minimize field, print result (if success)
        bool movingResetFlag = false;
        bool firstTry;
        try {
          firstTry = true;
          minimizer.minimize(*field, calculator);
        
          // if convergence is achieved: 
      
          // calculate final free energy:
            double freeEnergy = calculator.f(*field);
        
          // check to see if we've decayed to the disordered state:
          if (std::abs(freeEnergy) < 1e-8)
            movingResetFlag = true;
  
          // print free-energy:
	  freeEnergies.push_back(freeEnergy);
      
        // if we fail:	
        } catch (std::string errorMessage) {
        
          // try resetting field and minimizing again:
          if (firstTry) {
            provider.resetCondition(*field);
            firstTry = false;
          
            try {
              minimizer.minimize(*field, calculator);
  	      freeEnergies.push_back(calculator.f(*field));
  
            } catch (std::string errorMessage) { 
	      freeEnergies.push_back(100);
              movingResetFlag = true;
            } 
  
          } // end second try
        }
  
        // reset initial condition
        if (movingResetFlag) 
          provider.resetCondition(*field);
  
      } // end loop over gamma values
  
      // save free-energy:
    
      // reverse order of points if needed:
      freeEnergiesList.push_back(freeEnergies);

      // clear FieldProvider object for this phase
      delete(field);

      // end timer, print elapsed time
      /*
      auto end = std::chrono::steady_clock::now();
      std::cout << "phase " << phaseID 
                << " elapsed time: " 
		<< std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
		<< " s" << std::endl;
      */
  
  } // end loop over phaseIDs

  // print free-energies to file

  // create filename:
  std::stringstream gammaVal;
  gammaVal << std::fixed << std::setprecision(2) << gamma;
  std::string filename = "gamma_" + gammaVal.str() + "_phases";
  for (std::vector<int>::iterator it = phaseIDList.begin();
       it != phaseIDList.end();
       it++)
    filename += "_" + std::to_string(*it);

  // create output file
  std::ofstream stream(filename);
  stream << std::fixed << std::setprecision(8); 
  
  // loop through gamma values, print to file:
  double tau = -0.05;
  int index = 0;
  while (tau < 0.51) {
    stream << tau << ", " << gamma;

    for (doubleVecVec::iterator it = freeEnergiesList.begin();
         it != freeEnergiesList.end();
         it++)
    {
      doubleVec fP = *it;
      stream << ", " << fP[index];
    }
    
    stream << std::endl;

    tau += dTau;
    index++;
  }

  stream.close();

  return 0;
} // end main function
