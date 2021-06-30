#ifndef _BPG_MINIMIZER_GUARD
#define _BPG_MINIMIZER_GUARD

#include <cmath>
#include <complex.h>
#include <cstring>
#include <iostream>
#include <stdlib.h>

#include "fftw3.h"
#include "FieldProvider.h"

template<typename TFunctionalCalculator>
class BpgMinimizer {
  
  /*
   * ===========================================================================
   * 				 CLASS ATTRIBUTES
   * ===========================================================================
   */

private:
  double m_errorTolerance;	// used for convergence test
  int	 m_maxIterations;	// used to limit total no. of iterations
  int    m_maxFieldIterations  = 500; // max iterations of field opt. algorithm
  int    m_maxPeriodIterations = 1000;   // max iterations of period opt. algorithm
  int 	 m_iterator;		// used to track total iterations (alternating field/period optimization)
  int 	 m_fieldIterator;	// used to track no. of iterations of field minimization alg.
  int	 m_periodIterator;	// used to track no. of iterations of period optimization alg.

public:
  // constructor
  BpgMinimizer(double errorTolerance, int maxIterations)
  : m_errorTolerance(errorTolerance),
    m_maxIterations(maxIterations)
  {
    m_iterator       = 0;
    m_fieldIterator  = 0;
    m_periodIterator = 0;
  }

  // getter functions for iterators
  int getIterator()       { return m_iterator; }
  int getFieldIterator()  { return m_fieldIterator; }
  int getPeriodIterator() { return m_periodIterator; }

  /*
   * =============================================================================
   *                      FULL OPTIMIZATION (FIELD + PERIOD)
   * =============================================================================
   */
  void minimize( FieldProvider &field, TFunctionalCalculator &calculator )
  {
    // initialize iterator
    m_iterator = 0;

    /*
     * ================= iteration variables =====================
     */
    
    double fNew = calculator.f(field);
    double fOld = 0.0;
    double currentError = 0.0;
    bool stopCriterion = false;

    /*
     * ======================== loop =============================
     */
    while (!stopCriterion)
    {
      // increment iterator
      m_iterator++;
      std::cout << std::endl << "iteration no. " << m_iterator << std::endl;

      // optimize field (fixed period)
      optimizeField(field, calculator);
      std::cout << "free energy after field opt. " << calculator.f(field) << std::endl;
      
      // optimize period (fixed densities)
      optimizePeriods(field, calculator);
      std::cout << "free energy after period opt. " << calculator.f(field) << std::endl;

      // recompute free-energy
      fOld = fNew;
      fNew = calculator.f(field);

      // update error
      currentError = fNew - fOld;

      // update stop criterion
      if (std::abs(currentError) < m_errorTolerance || m_iterator == m_maxIterations)
        stopCriterion = true;
    }
    
    /*
     * ==================== finish up ============================
     */

    if (m_iterator == m_maxIterations) {
      std::string message = "maximum iterations reached";
      throw message;
    }
      
  } // end minimize method


  /*
   * =============================================================================
   * 			        FIELD OPTIMIZATION
   * =============================================================================
   */
  void optimizeField( FieldProvider &field, TFunctionalCalculator &calculator )
  {
    // start field optimization - set no. of iterations to zero
    m_fieldIterator = 0;


    /* 
     * =================== Initial field data =======================
     */
    fftw_complex* cplxFieldData = field.getCplxDataPointer();
    fftw_complex* realFieldData = field.getRealDataPointer();
    const int     d             = field.getDimension();
    const int*    gridSizes     = field.getGridSizes();
    const int     N             = field.getNumFieldElements();
    double*       dx            = field.getDx();
    double*       dq            = field.getDq(); 
    const int     phaseID       = field.getPhaseID();

    // allocate memory for Laplacian and initialize
    double* laplacian = (double*) malloc(N * sizeof(double));
    field.laplacian(laplacian);


    /*
     * ============== Nesterov field initialization ==================
     */ 
    
    // first need to make new grid sizes and spacings arrays - these are copied from initial field 
    int*    nestGridSizes = (int*)    malloc(d * sizeof(int));
    double* nestDq	      = (double*) malloc(d * sizeof(double));
    memcpy(nestGridSizes, gridSizes, d * sizeof(int));	// memcpy(dest, src, size)
    memcpy(nestDq,	      dq,	       d * sizeof(double));

    // create field
    FieldProvider nestField(
      cplxFieldData,    // field data, note: data works different from other array initializations
      d,                // dimension
      nestGridSizes,    // grid sizes
      nestDq,           // grid spacing
      false,            // initialized using cplx data
      phaseID);	        

    // get pointers to Nesterov field
    fftw_complex* realNestData = nestField.getRealDataPointer();
    fftw_complex* cplxNestData = nestField.getCplxDataPointer();


    /*
     * ============= NL derivative field initialization ============== 
     */

    // initialize data
    fftw_complex* initNlDerivData = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));

    // compute NL deriv from real field data
    calculator.nlDeriv(realNestData, initNlDerivData, N);

    // need to make new grid sizes/spacing arrays
    int*    nlDerivGridSizes = (int*)    malloc(d * sizeof(int));
    double* nlDerivDx        = (double*) malloc(d * sizeof(double));
    memcpy(nlDerivGridSizes, gridSizes, d * sizeof(int));
    memcpy(nlDerivDx,        dx,        d * sizeof(double));

    // create field 
    FieldProvider nlDerivField(
      initNlDerivData,
      d,
      nlDerivGridSizes,
      nlDerivDx,
      true,      // initialized using real data
      phaseID);      

    // free pointer to initializer data
    fftw_free(initNlDerivData);

    // get pointers to Nl deriv field
    fftw_complex* realNlDerivData = nlDerivField.getRealDataPointer();
    fftw_complex* cplxNlDerivData = nlDerivField.getCplxDataPointer();


    /*
     * ================= Quadratic coefficients ===================
     */

    // make and store array of quadratic coefficients - Gamma(q)
    double* Gamma = (double*) malloc(N * sizeof(double));
    calculator.makeGammaArray(Gamma, laplacian, N);


    /* 
     * ============== initialize iteration variables ==============
     */

    // free-energy
    double fNew = calculator.f(cplxFieldData, realFieldData, laplacian, N);
    double fOld = 0.0;

    // time step
    double tStep = 0.1;
    
    // Nesterov step variables
    double theta = 1.0;
    double beta  = 0.0;

    // error
    double currentError = 0.0;

    // boolean that determines whether to exit iterative loop
    bool stopCriterion = false;

    /* 
     * ====================== while loop ==========================  
     */
    
    while (!stopCriterion)
    {
      // increment iterator
      m_fieldIterator++;

      // update field
      for (int i = 0; i < N; i++) {
        // matrix element for update step
        double A = 1.0 - tStep * laplacian[i] * Gamma[i];

        // store old field vals
        double cplxFieldOldReVal = cplxFieldData[i][0];
        double cplxFieldOldImVal = cplxFieldData[i][1];

        // update real and imaginary parts of field
        cplxFieldData[i][0] = (cplxNestData[i][0] + tStep * laplacian[i] * cplxNlDerivData[i][0]) / A;
	cplxFieldData[i][1] = (cplxNestData[i][1] + tStep * laplacian[i] * cplxNlDerivData[i][1]) / A * 0;

        // update Nesterov field
        cplxNestData[i][0] = (1.0 + beta) * cplxFieldData[i][0] - beta * cplxFieldOldReVal;
        cplxNestData[i][1] = (1.0 + beta) * cplxFieldData[i][1] - beta * cplxFieldOldImVal;	
      }

      // transform field and Nesterov field to real space
      field.transformC2R();
      nestField.transformC2R();

      // update NL deriv field in real space and then transform to cplx space
      calculator.nlDeriv(realNestData, realNlDerivData, N);
      nlDerivField.transformR2C();

      // compute change in free-energy - used as error 
      fOld = fNew;
      fNew = calculator.f(cplxFieldData, realFieldData, laplacian, N);
      currentError = fOld - fNew;

      // update Nesterov coefficient - use test to damp oscillations
      if (currentError < 0) theta = 1.0;
      double theta2   = theta * theta;
      double newTheta = -0.5 * theta2 + std::sqrt(0.25 * theta2 * theta2 + theta2);
      beta = theta * (1 - theta) / (theta2 + newTheta);
      theta = newTheta;  

      // update stop criterion
      if (m_fieldIterator == m_maxFieldIterations || std::abs(currentError) < m_errorTolerance)
        stopCriterion = true;

    } // end of while loop

    /*
     * ====================== wrap it up! =============================
     */

    // if we've reached the max no. of iterations, throw an error
    if (m_fieldIterator == m_maxFieldIterations) {
      std::string message = "Maximum iterations reached in field optimization";
      throw message;
    }

    // free laplacian and quad coefficent arrays
    free(laplacian);
    free(Gamma);  

  } // end of field optimization method


  /*
   * =========================================================================================
   * 			              PERIOD OPTIMIZATION
   * =========================================================================================
   */

  void optimizePeriods( FieldProvider &field, TFunctionalCalculator &calculator) 
  {
    // start iterator at zero
    m_periodIterator = 0;

    /*
     * ==================== unpack field provider =======================
     */

    fftw_complex* cplxFieldData = field.getCplxDataPointer();
    const int     d             = field.getDimension();
    const int     N             = field.getNumFieldElements();
    const int*    gridSizes     = field.getGridSizes();
    double*       dQ            = field.getDq();

    /*
     * ======================= initializations ===========================
     */

    double* b    = (double*) calloc(3, sizeof(double)); // new reciprocal vec
    double* rOld = (double*) calloc(3, sizeof(double)); // old gradient vec
    double* rNew = (double*) calloc(3, sizeof(double)); // new gradient vec
    double* s    = (double*) calloc(3, sizeof(double)); // step direction vec

    /* 
     * ================== initialize reciprocal vector ==================
     */
    
    b[0] = dQ[d - 1];
    if (d > 1) {
      b[1] = dQ[d - 2];
      if (d > 2)
        b[2] = dQ[d - 3];
    }
    
    /*
     * =========== compute initial residual (gradient vector) ===========
     */
    
    computeGradient(rNew, b, calculator, cplxFieldData, d, N, gridSizes);
    double gradMagnitudeSquared = (rNew[0] * rNew[0]) + (rNew[1] * rNew[1]) + (rNew[2] * rNew[2]);

    /*
     * ====================== iteration variables =======================
     */
    double alpha = 1.0; // step size
    double beta  = 0.0; // conjugate coefficient

    double fNew = calculator.fQuad(field); // initial free-energy
    double fOld = 0.0;
    double currentError = 0.0;
    bool   posFlag = false; // becomes true if free-energy increases during loop
    double alphaMin = 0.1;  // minimum step-size (can be decreased if oscillations occur)

    // stop criterion 
    // note: if grad magnitude is small then we may not need to do anything 
    bool   stopCriterion = gradMagnitudeSquared < m_errorTolerance;   

    /*
     * ========================= while loop =============================
     */

    while (!stopCriterion) {
      m_periodIterator++;

      // compute line search direction s:
      s[0] = rNew[0] + beta * rOld[0];
      s[1] = rNew[1] + beta * rOld[1];
      s[2] = rNew[2] + beta * rOld[2];

      // compute optimal step-size (alpha) using backtracking line search
      // and update lattice vector b using optimal alpha
      // posFlag used here to shrink minimum allowed step size (crude)
      findAlpha(alpha, b, s, gradMagnitudeSquared, calculator, field, posFlag, alphaMin);

      // store residuals
      double gradMagnitudeOld = gradMagnitudeSquared;
      rOld[0] = rNew[0];
      rOld[1] = rNew[1];
      rOld[2] = rNew[2];

      // compute new residuals
      computeGradient(rNew, b, calculator, cplxFieldData, d, N, gridSizes);
      gradMagnitudeSquared = (rNew[0] * rNew[0]) + (rNew[1] * rNew[1]) + (rNew[2] * rNew[2]); 

      // compute conjugate coefficient for next time (Fletcher-Reeves formula)
      // use restart to minimize oscillations
      if (gradMagnitudeSquared > gradMagnitudeOld) beta = 0.0;
      else beta = gradMagnitudeSquared / gradMagnitudeOld;

      // compute change in free-energy - used as error
      fOld = fNew;
      fNew = calculator.fQuadB(field, b);
      currentError = fNew - fOld;

      if (currentError > 0) posFlag = true;
      else posFlag = false;

      // update stop criterion
      if (std::abs(currentError) < m_errorTolerance || m_periodIterator == m_maxPeriodIterations)
        stopCriterion = true;
    } // end while loop

    // if we've reached maximum iterations throw an error
    if (m_periodIterator == m_maxPeriodIterations) {
      std::string message = "Maximum iterations reached in period optimization";
      throw message;
    }

    // repackage reciprocal lattice vector b into format that field-provider likes
    double* newDq = (double*) malloc(d * sizeof(double));
    newDq[d-1] = b[0];
    if (d > 1) {
      newDq[d-2] = b[1];
      if (d > 2)
        newDq[d-3] = b[2];
    }

    // update field provider lattice spacings:
    field.setDq(newDq);	// note - do it this way because it also updates real spacing

    // free all arrays
    // free(bOld);
    free(b);
    free(rOld);
    free(rNew);
    free(s);   
    
  } // end of period optimization method

  /* 
   * ===================================================================
   * 		   period optimization - helper functions
   * ===================================================================
   */
private:

  /*
   * ===================================================================
   *  			COMPUTE GRADIENT VECTOR
   * ===================================================================
   */
  void computeGradient(
    double* r,
    double* b,
    TFunctionalCalculator &calculator,
    fftw_complex* cplxFieldData, 
    const int d,
    const int N, 
    const int* gridSizes) 
  {
    // initialize residuals 
    r[0] = 0.0;
    r[1] = 0.0;
    r[2] = 0.0;

    // loop over reciprocal lattice points
    for (int index = 0; index < N; index++) {
      
      // determine x, y, z indices of point:
      // need to also shift indicies to get negative values
      int i = 0, j = 0, k = 0;
      int iMod = 0, jMod = 0, kMod = 0;

      int Nx = gridSizes[d-1];
      i = index % Nx;
      iMod = i < Nx / 2 ? i : Nx - i;

      if (d > 1) {
        int Ny = gridSizes[d-2];
	j = (index - i) / Nx % Ny;
	jMod = j < Ny / 2 ? j : Ny - j;

	if (d > 2) {
	  int Nz = gridSizes[d-3];
	  k = ((index - i) / Nx - j) / Ny;
	  kMod = k < Nz / 2 ? k : Nz - k;
	}
      }

      // calculate squared reciprocal lattice vector (q2) at this point:
      double q2 = (iMod * iMod) * (b[0] * b[0]) +
	          (jMod * jMod) * (b[1] * b[1]) +
		  (kMod * kMod) * (b[2] * b[2]);

      // get value of derivative Gamma'(q2)
      double GammaPrime = calculator.quadraticCoeffDeriv(q2);

      // squared amplitude of fourier peak:
      double phi2 = cplxFieldData[index][0] * cplxFieldData[index][0] +
	            cplxFieldData[index][1] * cplxFieldData[index][1];

      // add value onto residual
      r[0] += -GammaPrime * phi2 * (iMod * iMod) * b[0];
      r[1] += -GammaPrime * phi2 * (jMod * jMod) * b[1];
      r[2] += -GammaPrime * phi2 * (kMod * kMod) * b[2];

    } // end loop over reciprocal lattice points
  } // end compute gradient method

  /*
   * ============================================================================
   *                              COMPUTE STEP SIZE
   * ============================================================================
   */

  void findAlpha(
    double& alpha,
    double* b,
    double* s,
    double gradMagnitude,
    TFunctionalCalculator& calculator,
    FieldProvider &field,
    bool posFlag,
    double &alphaMin)
  {
    // initialize alpha to max val.
    const double alphaMax = 1.0;
    alpha = alphaMax;

    // initial free-energy
    double f0 = calculator.fQuadB(field, b);

    // algorithm variables
    const double c        = 0.1;
    const double rho      = 0.9;

    // minimum step size:
    // decrease if free-energy increases OR if we've taken too many steps (may indicate oscillation)
    if (posFlag) alphaMin *= 0.1;
    if (m_periodIterator % 50 == 0) alphaMin *= 0.1;

    double* bNew = (double*) calloc(3, sizeof(double));
    bool stopCriterion = false;

    while (!stopCriterion) {
      // update b-vectors
      bNew[0] = b[0] + alpha * s[0];
      bNew[1] = b[1] + alpha * s[1];
      bNew[2] = b[2] + alpha * s[2];

      // calculate left and right-hand sides of inequality
      double lhs = f0 - calculator.fQuadB(field, bNew);
      double rhs = alpha * c * gradMagnitude;

      // test inequality:
      // if inequality holds or alpha is already too small - accept this value of alpha
      if (lhs >= rhs || alpha < alphaMin) stopCriterion = true;
      // otherwize shrink alpha and try again!
      else alpha *= rho;
    }

    // update step vector
    b[0] = bNew[0];
    b[1] = bNew[1];
    b[2] = bNew[2];

    free(bNew);
  } // end findAlpha method

}; // end class definition

#endif
