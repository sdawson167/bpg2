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
   *            		get tau value and phaseID 
   * ===============================================================================
   */ 
  
  double tau = 0.0;
  if (argc == 1) {
    // allow user to enter tau value from command line:
    std::cout << "choose tau value" << std::endl;

    std::string input;
    std::getline(std::cin, input);

    // make sure tau value is valid:
    std::stringstream stream(input);
    if (!(stream >> tau)) {
      std::cout << "invalid choice for tau parameter" << std::endl;
      return 1;
    }

  } else {
    // use tau value provided as input arg to program
    tau = std::atof(argv[1]);
  }

  int phaseID = 1;
  if (argc < 3) {
    // allow user to enter phaseID from command line:
    std::cout << "choose phase ID" << std::end;;
    
    std::string input;
    std::getline(std::cin, input);

    // make sure phaseID is valid:
    std::stringstream stream(input);
    if (!(stream >> phaseID)) {
      std::cout << "invalud choice for phaseID parameter" << std::endl;
      return 1;
    }

  } else {
    phaseID = std::atof(argv[2]);
  }


  /*
   * =============================================================================
   *    	      get gamma0 and make list of phase points
   * =============================================================================
   */
  
  // LB calculation:
  const double lam1 = 0.0;
  const double lam2 = 0.0;

  // get gamma0 val
  const double gamma0 = getGamma0(phaseID);
  const double dGamma = 0.1;

  // make list of phase points starting at tau, gamma0
  paramList phasePoints;

  std::vector<double> point0{ tau, gamma0, lam1, lam2 };
  phasePoints.push_back(point0);

  // first move backwards in phase space:
  if (gamma0 > 0.0) {
    double gamma = gamma0 - dGamma;

    while (gamma > -0.1) {
      std::vector<double> point{ tau, gamma, lam1, lam2 };
      phasePoints.push_back(point);

      gamma -= dGamma;
    }

    // redo gamma0 calc:
    phasePoints.push_back(point0);
  }

  // now move forward in phase space:
  if (gamma0 < 2.5) {
    double gamma = gamma0 + dGamma;

    while (gamma < 2.6) {
      std::vector<double> point{ tau, gamma, lam1, lam2 };
      phasePoints.push_back(point);

      gamma += dGamma;
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
  
  // ensure we don't print more digits of precion than we know (based on errorTol):
  std::cout << std::fixed << std::setprecision(8);

  // create minimizer (object that handles implementation of minimization algorithm) 
  // BpgMinimizer<PwFunctionalCalculator> minimizer(errorTol, maxIter);

  // typedef std::vector<double> doubleVec;
  // std::vector<doubleVec> freeEnergies;
   
  // create field-provider object (represents the order parameter in real and fourier space) 
  // initialize field-provider in chosen phase
  GenericPhaseProvider provider(phaseID); 
  // FieldProvider *field = provider.generateInitialCondition();

  // vector of free-energies:
  // doubleVec fP;

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

    // print phase point
    std::cout << tau << ", " << gamma << std::endl;

    // create 'functional calculator' object - handles free-energy calculations
    // PwFunctionalCalculator calculator(tau, gamma, lam1, lam2);

    // minimize field, print result (if success) or '100' (if failure)
    /*
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
    */

  } // end loop over phase points

  // add list of free-energies to vector
  //freeEnergies.push_back(fP);

  // clear FieldProvider object for this phase
  // delete(field);

  // print free-energies to command line:
  /*
  for (size_t i = 0; i < phasePoints.size(); i++) {
    std::vector<double> x = phasePoints[i];
    std::cout << x[0] << ", " << x[1];

    for (size_t j = 0; j < phaseIDList.size(); j++) { 
      doubleVec fP = freeEnergies[j];
      std::cout << ", " << fP[i];
    }
    std::cout << std::endl;
  }
  */
  return 0;
}
