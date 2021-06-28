#ifndef _PW_CALCULATOR_GUARD
#define _PW_CALCULATOR_GUARD

#include "FieldProvider.h"
#include "fftw3.h"

/*
 * Class represents the piece-wise free-energy functional
 * this model computes the (q<1) and (q>1) parts of the quadratic coefficient 
 * independently, with both halves interpolating between the LB and OK models
 *
 *   tau   -- quadratic coefficient for nl part of free-energy
 *   gamma -- cubic coefficient for nl part of free-energy
 *   lam1  -- interpolation parameter for q < 1 (lam1 = 0.0 : LB, lam1 = 1.0 : OK)
 *   lam2  -- interpolation parameter for q > 1 (lam2 = 0.0 : LB, lam2 = 1.0 : OK)
 */

class PwFunctionalCalculator {

private:
    double  m_tau;
    double  m_gamma;
    double  m_lam1;
    double  m_lam2;

    const int m_paramNum = 4;

public:
    // constructor
    PwFunctionalCalculator(double tau, double gamma, double lam1, double lam2)
    : m_tau{tau},
      m_gamma{gamma},
      m_lam1{lam1},
      m_lam2{lam2}
    { }

    // getter and setter functions
    double getTau() { return m_tau; }
    void setTau(double tau) { m_tau = tau; }

    double getGamma() { return m_gamma; }
    void setGamma(double gamma) { m_gamma = gamma; }

    double getLam1() { return m_lam1; }
    void setLam1(double lam1) { m_lam1 = lam1; }

    double getLam2() { return m_lam2; }
    void setLam2(double lam2) { m_lam2 = lam2; }

    int getParamNum() { return m_paramNum; }

    // find quadratic potential coefficients, D(q)
    double lbCoeff(double q2) { 
      return (q2 - 1) * (q2 - 1); 
    }
    double okCoeff(double q2) { 
      if (q2 == 0) return 0.0;
      else return q2 + 1.0 / q2 - 2; 
    }
    double quadraticCoeff(double q2) {
      if (q2 <= 1.0)
	return m_lam1 * okCoeff(q2) + (1.0 - m_lam1) * lbCoeff(q2);
      else
	return m_lam2 * okCoeff(q2) + (1.0 - m_lam2) * lbCoeff(q2);
    }
    void   makeGammaArray(double* Gamma, double* laplacian, int N);

    // compute quadratic free-energy
    double fQuad(FieldProvider &field);
    double fQuad(fftw_complex* cplxFieldData, double* laplacian, int numFieldElements);

    // compute non-linear free-energy
    double fNL(FieldProvider &field);
    double fNL(fftw_complex* realFieldData, int numFieldElements);

    // compute full free-energy
    double f(FieldProvider &field) {return fNL(field) + fQuad(field);}
    double f(fftw_complex* cplxFieldData, fftw_complex* realFieldData, double* laplacian, int numFieldElements)
    {
      return fQuad(cplxFieldData, laplacian, numFieldElements) + fNL(realFieldData, numFieldElements);
    }

    // compute derivative of NL free-energy and store in field
    void nlDeriv(fftw_complex* realFieldData, fftw_complex* realNLFieldData, int numFieldElements);

    // period optimization method
    double lbDeriv(double q2) {
      return 2 * (q2 - 1);
    }
    double okDeriv(double q2) {
      if (q2 == 0) return 0;
      else return 1.0 - 1.0 / (q2 * q2);
    }
    double quadraticCoeffDeriv(double q2) { 
      if (q2 <= 1.0)
        return m_lam1 * okDeriv(q2) + (1.0 - m_lam1) * lbDeriv(q2);
      else
	return m_lam2 * okDeriv(q2) + (1.0 - m_lam2) * lbDeriv(q2);
    }

    // compute quadratic free-energy for given lattice spacings
    double fQuadB(FieldProvider &field, double* b);
};
#endif
