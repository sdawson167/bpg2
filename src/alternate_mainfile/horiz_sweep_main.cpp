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

  const double lam1 = 1.0;
  const double lam2 = 1.0;


  /*
   * ===============================================================================
   *                      get tau value and phaseID(s) 
   * ===============================================================================
   */ 
  
  double tau;
  std::vector<int> phaseIDList;

  // make sure there are enough input args to extract tau and phaseID vals:
  if (argc == 1) {
    std::cout << "error: must choose tau value and one or more phase IDs (1 - 7) to run program." << std::endl;

    return 1;
  
  } else if (argc == 2)  {
    std::cout << "error: must choose one or more phase IDs (1 - 7) to run program." << std::endl; 

    return 1;

  } else {
    tau = std::atof(argv[1]);

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
  
  // horizontal resolution: 
  const double dGamma = 0.1;

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
    
      // get starting value of gamma for this phase: 
      const double gamma0 = getGamma0(phaseID);

      // make list of gamma values starting at gamma0:
      doubleVec gammaVec;
      
      gammaVec.push_back(gamma0);

      // first try to move fwd (right) in phase space:
      if (gamma0 < 2.0) {
        double gamma = gamma0 + dGamma;

        while (gamma < 2.01) {
	  gammaVec.push_back(gamma);

	  gamma += dGamma;
        }
      }

      // then try to move bwd (left) in phase space:
      if (gamma0 > 0.0) {
        double gamma = gamma0 - dGamma;

        while (gamma > -0.01) {
	  gammaVec.push_back(gamma);

	  gamma -= dGamma;
        }
      }

      // create field-provider object (represents the order parameter in real and fourier space) 
      // initialize field-provider in chosen phase
      GenericPhaseProvider provider(phaseID); 
      FieldProvider *field = provider.generateInitialCondition();

      // create vector to store free-energy for this phase
      doubleVec freeEnergies;

      // loop through phase points (tau, gamma, ... vals) - compute free-energy at each point
      for (std::vector<double>::const_iterator it = gammaVec.begin(); it != gammaVec.end(); it++)
      {
        // get gamma value
        double gamma = *it; 
  
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
      if (gamma0 > 0.0) {
        doubleVec flippedFreeEnergies;
        for (doubleVec::reverse_iterator it = freeEnergies.rbegin(); 
             it != freeEnergies.rend();
	     it++) 
          flippedFreeEnergies.push_back(*it);

        freeEnergiesList.push_back(flippedFreeEnergies);
      } else 
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
  std::stringstream tauVal;
  tauVal << std::fixed << std::setprecision(2) << tau;
  std::string filename = "tau_" + tauVal.str() + "_phases";
  for (std::vector<int>::iterator it = phaseIDList.begin();
       it != phaseIDList.end();
       it++)
    filename += "_" + std::to_string(*it);

  // create output file
  std::ofstream stream(filename);
  stream << std::fixed << std::setprecision(8); 
  
  // loop through gamma values, print to file:
  double gamma = 0.0;
  int index = 0;
  while (gamma < 2.1) {
    stream << tau << ", " << gamma;

    for (doubleVecVec::iterator it = freeEnergiesList.begin();
         it != freeEnergiesList.end();
         it++)
    {
      doubleVec fP = *it;
      stream << ", " << fP[index];
    }
    
    stream << std::endl;

    gamma += dGamma;
    index++;
  }

  stream.close();

  return 0;
} // end main function
