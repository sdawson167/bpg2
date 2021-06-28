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
   *             Parse input args - determine which calculation to do
   * ===============================================================================
   */
  
  paramList phasePoints;
  std::vector<int> phaseIDList;
  bool resetFlag = 1;

  // attempt to extract list of calculations to perform 
  try {
    parseInputArgs(argc, argv, phasePoints, phaseIDList, resetFlag);

  } catch (std::string errorMessage){
    std::cout << errorMessage << std::endl;
    return 1;
  }

  /*
   * ==============================================================================
   *               Start optimization for all points listed in input
   * ==============================================================================
   */

  // minimization parameters (max iterations and error tolerance)
  const double errorTol = 1e-9;
  const int    maxIter  = 15;
  
  // ensure we don't print more digits of precion than we know (based on errorTol):
  std::cout << std::fixed << std::setprecision(8);

  // create minimizer (object that handles implementation of minimization algorithm) 
  BpgMinimizer<PwFunctionalCalculator> minimizer(errorTol, maxIter);

  typedef std::vector<double> doubleVec;
  std::vector<doubleVec> freeEnergies;
  // loop through phases
  
  for (std::vector<int>::const_iterator idIter = phaseIDList.begin(); idIter != phaseIDList.end(); idIter++) { 
   
    // get phaseID - integer that indicates which phase we are optimizing
    int phaseID = *idIter;

    // create field-provider object (represents the order parameter in real and fourier space) 
    // initialize field-provider in chosen phase
    GenericPhaseProvider provider(phaseID); 
    FieldProvider *field = provider.generateInitialCondition();

    // vector of free-energies:
    doubleVec fP;

    // loop through phase points (tau, gamma, ... vals) - compute free-energy at each point
    
    for (paramList::const_iterator it = phasePoints.begin(); it != phasePoints.end(); it++)
    {
      // get phasePoint
      std::vector<double> otherParams = *it;

      // get tau/gamma and lam1/lam2 values
      double tau   = otherParams[0];
      double gamma = otherParams[1];
      double lam1  = otherParams[2];
      double lam2  = otherParams[3];

      // create 'functional calculator' object - handles free-energy calculations
      PwFunctionalCalculator calculator(tau, gamma, lam1, lam2);

      // minimize field, print result (if success) or '100' (if failure)
      bool movingResetFlag = resetFlag;
      bool firstTry;
      try {
	firstTry = true;
        minimizer.minimize(*field, calculator);
	
        // if we succeed - save free-energy:
	fP.push_back(calculator.f(*field));
      
      // if we fail:	
      } catch (std::string errorMessage) {
 	
	// try resetting field and minimizing again:
        if (firstTry && !resetFlag) {
	  provider.resetCondition(*field);
	  firstTry = false;
	  
	  try {
	    minimizer.minimize(*field, calculator);
	    fP.push_back(calculator.f(*field));

	  } catch (std::string errorMessage) { 
	    fP.push_back(100);
	    movingResetFlag = true;
	  } 

        } // end second try
      }

      // reset initial condition
      if (movingResetFlag) 
        provider.resetCondition(*field);

    } // end loop over phase points

    // add list of free-energies to vector
    freeEnergies.push_back(fP);

    // clear FieldProvider object for this phase
    delete(field);

  } // end loop over phases
  
  // print free-energies to command line:
  for (size_t i = 0; i < phasePoints.size(); i++) {
    std::vector<double> x = phasePoints[i];
    std::cout << x[0] << ", " << x[1];

    for (size_t j = 0; j < phaseIDList.size(); j++) { 
      doubleVec fP = freeEnergies[j];
      std::cout << ", " << fP[i];
    }
    std::cout << std::endl;
  }

  return 0;
}
