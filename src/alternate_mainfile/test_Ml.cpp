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

//#include "PwFunctionalCalculator.h"
#include "MlFunctionalCalculator.h"

int main( )
{
  // minimization parameters (max iterations and error tolerance)
  const double errorTol = 1e-7;
  const int    maxIter  = 40;
  
  // ensure we don't print more digits of precion than we know (based on errorTol):
  // std::cout << std::fixed << std::setprecision(10);

  // create minimizer (object that handles implementation of minimization algorithm) 
  BpgMinimizer<MlFunctionalCalculator> minimizer(errorTol, maxIter);

  // get phaseID - integer that indicates which phase we are optimizing
  int phaseID = 7;
 
  // create field-provider object (represents the order parameter in real and fourier space) 
  // initialize field-provider in chosen phase
  GenericPhaseProvider provider(phaseID); 
  FieldProvider *field = provider.generateInitialCondition();

  std::cout << "minimizing phase " << phaseID << std::endl;

  // get phasePoint
  std::vector<double> otherParams = *it;

  // get tau/gamma and eps/f values
  double tau   = 0.0; 
  double gamma = 1.0;
  double eps   = 1.0;
  double f     = 0.5;

  // create 'functional calculator' object - handles free-energy calculations
  MlFunctionalCalculator calculator(tau, gamma, eps, f);

  std::cout << "initial free energy: " << calculator.f(*field) << std::endl;

  // minimize field, print result (if success) or error message (if failure)
  //auto start = std::chrono::steady_clock::now();
  try {
    minimizer.minimize(*field, calculator);
  
    // if we succeed - print free-energy:
    std::cout << tau << ", " << gamma << ", " << calculator.f(*field) << std::endl;

  // if we fail:	
  } catch (std::string errorMessage) {
    std::cout << tau << ", " << gamma << ", " << 100 << std::endl;
  }
  
  //auto end = std::chrono::steady_clock::now();
  //std::cout << "elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << std::endl;
  
  // clear FieldProvider object for this phase
  delete(field);

  return 0;
}
