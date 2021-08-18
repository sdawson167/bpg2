#include <cmath>

#include "MlFunctionalCalculator.h"

/*
 *  ============================================================
 *                     Quadratic free-energy
 *  ============================================================
 */

// function to compute the value of q^2 at which the unscaled gamma function is
// minimized
// this minimization will be done using Newton's method
double MlFunctionalCalculator::findQ02()
{
  const int maxIterator = 1000;
  const double errorTol = 1e-7;

  double x0 = 1.0;

  bool stopCriterion = false;
  int  iterator = 0;
  while (!stopCriterion) {
    iterator++;

    double A =  dUnscaledGamma(x0);
    double B = d2UnscaledGamma(x0);

    double x = x0 - A/B;
    
    double error = std::abs(x - x0);
    
    if (error < errorTol || iterator == maxIterator)
      stopCriterion = true;  

    x0 = x;
  }

  if (iterator == maxIterator)
    throw "maximum iterations reached when attempting to find minimum of ML quadratic coefficient";

  return x0;
}

// compute quadratic free-energy from field and laplacian (requires laplacian initialization)
double MlFunctionalCalculator::fQuad(fftw_complex* cplxFieldData, double* laplacian, int numFieldElements)
{
  double fQuad = 0.0;
  for (int index = 0; index < numFieldElements; index++) {
    auto fieldVal = cplxFieldData[index];
    double fieldVal2 = (fieldVal[0] * fieldVal[0]) + (fieldVal[1] * fieldVal[1]);

    double coeff = quadraticCoeff(-laplacian[index]);

    fQuad += 0.5 * coeff * fieldVal2;
  }

  return fQuad;
}

void MlFunctionalCalculator::makeGammaArray(double* Gamma, double* laplacian, int N)
{
  for (int index = 0; index < N; index++)
    Gamma[index] = quadraticCoeff(-laplacian[index]);
}

// compute quadratic free-energy from field provider (initializes laplacian)
double MlFunctionalCalculator::fQuad(FieldProvider &field)
{
  fftw_complex* cplxField = field.getCplxDataPointer();

  const int numFieldElements = field.getNumFieldElements();

  // initialize laplacian
  double* laplacian = (double*) malloc(sizeof(double) * numFieldElements);
  field.laplacian(laplacian);

  double f = fQuad(cplxField, laplacian, numFieldElements);

  free(laplacian);

  return f;
}

/*
 *  ============================================================
 *                    non-linear free-energy
 *  ============================================================
 */

// compute non-linear free-energy from real field pointer
double MlFunctionalCalculator::fNL(fftw_complex* realFieldData, int numFieldElements)
{
  double fNL = 0.0;
  for (int index = 0; index < numFieldElements; index++) {
    double fieldVal  = realFieldData[index][0];
    double fieldVal2 = fieldVal * fieldVal;
    double fieldVal3 = fieldVal * fieldVal2;
    double fieldVal4 = fieldVal * fieldVal3;

    fNL += (m_tau / 2) * fieldVal2 - (m_gamma / 6) * fieldVal3 + (1.0 / 24) * fieldVal4;
  }

  fNL *= 1.0 / numFieldElements;
  return fNL;
}

// compute non-linear free-energy from field provider
double MlFunctionalCalculator::fNL(FieldProvider &field)
{
  fftw_complex* realField = field.getRealDataPointer();
  int numFieldElements = field.getNumFieldElements();
  return fNL(realField, numFieldElements);
}

/*
 *  ============================================================
 *                  derivative of NL free-energy
 *  ============================================================
 */
void MlFunctionalCalculator::nlDeriv(
  fftw_complex* realFieldData,
  fftw_complex* realNLFieldData,
  int numFieldElements)
{
  for (int index = 0; index < numFieldElements; index++) {
    double fieldVal  = realFieldData[index][0];
    double fieldVal2 = fieldVal * fieldVal;
    double fieldVal3 = fieldVal * fieldVal2;

    realNLFieldData[index][0] = (m_tau) * fieldVal - (m_gamma / 2.0) * fieldVal2 + (1.0 / 6.0) * fieldVal3;
    realNLFieldData[index][1] = 0.0;
  }
}

/*
 *  ===========================================================
 *        free-energy as a function of lattice spacing
 *  ===========================================================
 */
double MlFunctionalCalculator::fQuadB(FieldProvider &field, double* b)
{
  // total number of lattice points:
  int numFieldElements = field.getNumFieldElements();

  // dimension of field:
  int d = field.getDimension();

  // vector of lattice sizes in z, y, x-directions:
  int* N = field.getGridSizes();

  // squared lattice spacings and grid sizes:
  int Nx{0}, Ny{0}, Nz{0};
  double bx2{0.0}, by2{0.0}, bz2{0.0};
  Nx  = N[d - 1];
  bx2 = b[0] * b[0];
  if (d > 1) {
    Ny  = N[d - 2];
    by2 = b[1] * b[1];
    if (d > 2) {
      Nz  = N[d - 3];
      bz2 = b[2] * b[2];
    }
  }

  // field data  (fourier peaks \hat{phi}_n)
  fftw_complex* cplxFieldData = field.getCplxDataPointer();

  // laplacian
  double* laplacian = (double*) fftw_malloc(sizeof(double) * numFieldElements);
  field.laplacian(laplacian);

  // compute quadratic free-energy by summing over reciprocal lattice points:
  double fb{0.0};
  for (int index = 0; index < numFieldElements; index++) {

    // get phi2 val:
    double phi2 = (cplxFieldData[index][0] * cplxFieldData[index][0]) + (cplxFieldData[index][1] * cplxFieldData[index][1]);

    // determine x, y, z indices of lattice point:
    int i{0}, j{0}, k{0};
    int iMod{0}, jMod{0}, kMod{0};

    i = index % Nx;
    iMod = i < Nx / 2 ? i : Nx - i;

    if (d > 1) {
      j = (index - i) / Nx % Ny;
      jMod = j < Ny / 2 ? j : Ny - j;

      if (d > 2) {
        k = ((index - i) / Nx - j) / Ny;
        kMod = k < Nz / 2 ? k : Nz - k;
      }
    }

    // calculate q2:
    double q2 = (iMod * iMod) * bx2 + (jMod * jMod) * by2 + (kMod * kMod) * bz2;

    // calculate free-energy contribution from this mode:
    fb += 0.5 * quadraticCoeff(q2) * phi2;

  } // end loop over modes

  free(laplacian);

  return fb;
} // end fQuadB function

