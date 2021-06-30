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

int main( )
{
  GenericPhaseProvider initializer("c15");

  FieldProvider *field = initializer.generateInitialCondition();
  //field.saveForPlotting("csv_files/c140.csv");

  // parameters
  double tau = -0.2;
  double gamma = 1.2;
  double lam1 = 0.0;
  double lam2 = 0.0;

  // functional calculator
  PwFunctionalCalculator calculator(tau, gamma, lam1, lam2);

  std::cout << "initial free energy " << calculator.f(*field) << std::endl << std::endl;
  
  // minimization parameters (max iterations and error tolerance)
  const double errorTol = 1e-7;
  const int    maxIter  = 40;
  
  // create minimizer (object that handles implementation of minimization algorithm) 
  BpgMinimizer<PwFunctionalCalculator> minimizer(errorTol, maxIter);

  // minimize field
  minimizer.minimize(*field, calculator);
  std::cout << "final free energy " << calculator.f(*field) << std::endl << std::endl; 

  // save final field to plot
  //field.saveForPlotting("csv_files/final.csv");

  // reset field
  (*field).reset();
  std::cout << "free energy after reset " << calculator.f(*field) << std::endl;

  delete(field);

  return 0;
}
