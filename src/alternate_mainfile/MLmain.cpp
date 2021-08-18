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

#include "MlFunctionalCalculator.h"

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
    parseInputArgsMl(argc, argv, phasePoints, phaseIDList, resetFlag);

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
  const double errorTol = 1e-7;
  const int    maxIter  = 40;
  
  // ensure we don't print more digits of precion than we know (based on errorTol):
  // std::cout << std::fixed << std::setprecision(10);

  // create minimizer (object that handles implementation of minimization algorithm) 
  BpgMinimizer<MlFunctionalCalculator> minimizer(errorTol, maxIter);

  // loop through phases
  for (std::vector<int>::const_iterator idIter = phaseIDList.begin(); idIter != phaseIDList.end(); idIter++) { 
   
    // get phaseID - integer that indicates which phase we are optimizing
    int phaseID = *idIter;

    // create field-provider object (represents the order parameter in real and fourier space) 
    // initialize field-provider in chosen phase
    GenericPhaseProvider provider(phaseID); 
    FieldProvider *field = provider.generateInitialCondition();

    // flag used to track if phase was freshly initialized: 
    bool reset = true;

    // create file to store results 
    std::stringstream stream;
    stream << "ML_phase" << phaseID << ".csv";
    std::string fileName = stream.str();
    
    std::ofstream outFileStream;
    outFileStream.open(fileName, std::ios_base::app);
    outFileStream << std::fixed << std::setprecision(8);

    // loop through phase points (tau, gamma, ... vals) - compute free-energy at each point
    for (paramList::const_iterator it = phasePoints.begin(); it != phasePoints.end(); it++)
    {
      // get phasePoint
      std::vector<double> otherParams = *it;

      // get tau/gamma and lam1/lam2 values
      double tau   = otherParams[0];
      double gamma = otherParams[1];
      double eps  = otherParams[2];
      double f  = otherParams[3];
      bool movingResetFlag = resetFlag;
      
      if ((phaseID == 2 && (tau > 0.0 || gamma > 0.7)) || (phaseID > 6 && gamma < 0.5))
        outFileStream << 100 << std::endl;
	//std::cout << 100 << std::endl;
      else {
        // create 'functional calculator' object - handles free-energy calculations
        MlFunctionalCalculator calculator(tau, gamma, eps, f);

        // flag use to track if we should reset the phase between points:
        bool movingResetFlag = resetFlag;

        // minimize field, print result (if success) or error message (if failure)
        //auto start = std::chrono::steady_clock::now();
        try {
          minimizer.minimize(*field, calculator);
	
          // if we succeed - print free-energy:
          //std::cout << calculator.f(*field) << std::endl;
	  outFileStream << calculator.f(*field) << std::endl;

	  // we can use this phase for the next loop
	  reset = false;
      
        // if we fail:	
        } catch (std::string errorMessage) {
 	
	  // if we didn't start with a freshly initialized phase, try resetting and minimizing again:
          if (!reset) {
	    provider.resetCondition(*field);
	    reset = true;
	  
	    try {
	      minimizer.minimize(*field, calculator);

	      // if we succeed - print free-energy:
	      //std::cout << calculator.f(*field) << std::endl;
	      outFileStream << calculator.f(*field) << std::endl;

	      // we can use this phase for the next loop
	      reset = false;

	    } catch (std::string errorMessage) { 
	      // if we fail the second time, give up on this point
	      //std::cout << 100 << std::endl;
	      outFileStream << 100 << std::endl;
	    
	      // need to reset before the next loop
	      movingResetFlag = true;
	    }
          } // end second chance

	  // if we initialized the phase at the beginning of the loop there's nothing to do 
	  // except give up on this point :(
	  else {  
            //std::cout << 100 << std::endl;
	    outFileStream << 100 << std::endl;
	  }  
        }


        //auto end = std::chrono::steady_clock::now();
        //std::cout << "elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << std::endl;

        // reset initial condition
        if (movingResetFlag) { 
          provider.resetCondition(*field);
	  reset = true;
        }

      } // end loop over phase points
    
    }
    // clear FieldProvider object for this phase
    delete(field);

  } // end loop over phases
  
  return 0;
}
