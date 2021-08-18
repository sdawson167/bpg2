#ifndef _ML_CALCULATOR_GUARD
#define _ML_CALCULATOR_GUARD


#include "FieldProvider.h"
#include "fftw3.h"

/*
 * Class represents the modified Liebler free-energy functional, which describes
 * an AB-type diblock copolymer with segment asymmetry
 *
 *   tau   -- quadratic coefficient for nl part of free-energy
 *   gamma -- cubic coefficient for nl part of free-energy
 *   eps   -- squared ratio of Kuhn lengths - represents asymmetry of chains
 *   f     -- fraction of diblock chain composed of A-type polymer segments
 */

class MlFunctionalCalculator {

private:
    double  m_tau;
    double  m_gamma;
    double  m_eps;
    double  m_f;

    double m_q02;  // wavevector for which unscaled quadratic coefficient is minimized
    double m_scaleFactor; // rescaling of quadratic coefficient

    const int m_paramNum = 4;

public:
    // constructor
    MlFunctionalCalculator(double tau, double gamma, double eps, double f)
    : m_tau{tau},
      m_gamma{gamma},
      m_eps{eps},
      m_f{f}
    { 
      m_q02         = findQ02();
      m_scaleFactor = findScaleFactor();
    }

    // getter and setter functions
    double getTau() { return m_tau; }
    void setTau(double tau) { m_tau = tau; }

    double getGamma() { return m_gamma; }
    void setGamma(double gamma) { m_gamma = gamma; }

    double getEps() { return m_eps; }
    void setEps(double eps) { m_eps = eps; }

    double getF() { return m_f; }
    void setF(double f) { m_f = f; }

    double getQ02() { return m_q02; }

    double getScaleFactor() { return m_scaleFactor; }

    int getParamNum() { return m_paramNum; }

    /*
     * compute quadratic coefficient - Gamma(q)
     */

    // helper functions used to compute correlation functions and their
    // derivatives:
    
    double   g(double alpha, double y) { return 2 * ((alpha * y) + exp(-alpha * y) - 1) / (y * y); }
    //derivative of g wrt y:
    double  dg(double alpha, double y) { return 4 * (1 - (alpha * y / 2) - (1 + alpha * y / 2) * exp(-alpha * y)) / (y * y * y); }
    // second derivative of g wrt y:
    double d2g(double alpha, double y) { return -12 * (1 - (alpha * y / 3) - (1 + (2 * alpha * y / 3) + (alpha * alpha * y * y / 6)) * exp(-alpha * y)) / (y * y * y * y); }

    double   h(double alpha, double y) { return (1 - exp(-alpha * y)) / y; }
    //derivative of h wrt y:
    double  dh(double alpha, double y) { return ((1 + alpha * y) * exp(-alpha * y) - 1) / (y * y);}
    // second derivative of h wrt y:
    double d2h(double alpha, double y) { return 2 * (1 - (1 + (alpha * y) + (alpha * alpha * y * y / 2)) * exp(-alpha * y)) / (y * y * y); }

    // correlation functions and their derivatives

    // AA-correlation function
    double   SAA(double q2) { return g(m_f, m_eps * q2); }
    // derivative of AA-correlation function wrt q2
    double  dSAA(double q2) { return m_eps * dg(m_f, m_eps * q2); }
    // second derivative of AA-correlation function wrt q2
    double d2SAA(double q2) { return m_eps * m_eps * d2g(m_f, m_eps * q2); }

    // BB-correlation function
    double   SBB(double q2) { return g(1 - m_f, q2); }
    // derivative of BB-correlation function wrt q2
    double  dSBB(double q2) { return dg(1 - m_f, q2); }
    // second derivative of BB-correlation function wrt q2
    double d2SBB(double q2) { return d2g(1 - m_f, q2); }

    // AB-correlation function
    double   SAB(double q2) { return h(m_f, m_eps * q2) * h(1 - m_f, q2); }
    // derivative of AB-correlation function wrt q2
    double  dSAB(double q2) { return m_eps * dh(m_f, m_eps * q2) * h(1 - m_f, q2) + h(m_f, m_eps * q2) * dh(1 - m_f, q2); }
    // second derivative of AB-correlation function wrt q2
    double d2SAB(double q2) { 
      return m_eps * m_eps * d2h(m_f, m_eps * q2) * h(1 - m_f, q2) + 
             2 * m_eps * dh(m_f, m_eps * q2) * dh(1 - m_f, q2) + 
             h(m_f, m_eps * q2) * d2h(1 - m_f, q2); 
    }

    // combinations of correlation functions (and their derivatives) found by
    // RPA approximation to quadratic coefficient

    double   S(double q2) { return SAA(q2) + 2 * SAB(q2) + SBB(q2);}
    // derivative of S wrt q2
    double  dS(double q2) { return dSAA(q2) + 2 * dSAB(q2) + dSBB(q2); }
    // second derivative of S wrt q2
    double d2S(double q2) { return d2SAA(q2) + 2 * d2SAB(q2) + d2SBB(q2); }

    double   W(double q2) { return SAA(q2) * SBB(q2) - SAB(q2) * SAB(q2); }
    // derivative of W wrt q2
    double  dW(double q2) { return dSAA(q2) * SBB(q2) + SAA(q2) * dSBB(q2) - 2 * SAB(q2) * dSAB(q2); }
    // second derivative of W wrt q2
    double d2W(double q2) { return d2SAA(q2) * SBB(q2) + 2 * dSAA(q2) * dSBB(q2) + SAA(q2) * d2SBB(q2) - 2 * (dSAB(q2) * dSAB(q2) + SAB(q2) * d2SAB(q2)); }

    // unscaled quadratic coefficient
    double   unscaledGamma(double q2) { return S(q2) / W(q2); }
    // derivative of unscaled quadratic coefficient wrt q2
    double  dUnscaledGamma(double q2) {
      if (q2 == 0) return 0.0;
      else return (dS(q2) * W(q2) - S(q2) * dW(q2)) / (W(q2) * W(q2));
    }
    // second derivative of unscalled coeff wrt q2
    double d2UnscaledGamma(double q2) {
      if (q2 == 0) return 0.0; 
      else return (W(q2) * (d2S(q2) * W(q2) - S(q2) * d2W(q2)) - 2 * dW(q2) * (dS(q2) * W(q2) - S(q2) * dW(q2)) ) / (W(q2) * W(q2) * W(q2));
    }

    // function to find minimum value of q2 - used to rescale lengths in theory
    double findQ02();

    // function used to calculate scale factor of quadratic coefficient
    double findScaleFactor() {
      return 2 / ((m_q02 * m_q02) * d2UnscaledGamma(m_q02)); 
    }

    // return fully scaled quadratic coefficient
    double quadraticCoeff(double q2) {
      if (q2 == 0)
        return 0;
      else {
        double scaledQ = m_q02 * q2;
        return m_scaleFactor * (unscaledGamma(scaledQ) - unscaledGamma(m_q02));
      }
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

    // compute derivative of quadratic coefficient - Gamma'(q2) 
    // needed for period optimization code in BPG algorithm
    double quadraticCoeffDeriv(double q2) { 
      double scaledQ = m_q02 * q2;
      return m_scaleFactor * m_q02 * dUnscaledGamma(scaledQ);
    }

    // compute quadratic free-energy for given lattice spacings
    double fQuadB(FieldProvider &field, double* b);
};
#endif
