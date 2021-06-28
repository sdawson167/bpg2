#ifndef _ML_CALCULATOR_GUARD
#define _ML_CALCULATOR_GUARD

#include <cmath>
#include "FieldProvider.h"
#include "fftw3.h"

/*
 * Class contains information about the Modified Liebler free-energy functional
 */

class MlFunctionalCalculator {

private:
    double m_tau;
    double m_gamma;
    double m_f;
    double m_epsilon;

    double m_q02;
    double m_scaleFactor;

    const int m_paramNum = 4;

public:
    // constructor
    MlFunctionalCalculator(double tau, double gamma, double f, double epsilon)
    : m_tau{tau},
      m_gamma{gamma},
      m_f{f},
      m_epsilon{epsilon}
    {
      m_q02 = computeQ0();
      m_scaleFactor = computeScaleFactor(); 
    }

    // getter and setter functions
    double getTau() { return m_tau; }
    void setTau(double tau) { m_tau = tau; }

    double getGamma() { return m_gamma; }
    void setGamma(double gamma) { m_gamma = gamma; }

    int getParamNum() { return m_paramNum; }

    double  g(double y, double alpha) { return 2 * (alpha * y + std::exp(-alpha * y) - 1) / (y * y); }
    double dg(double y, double alpha) { return 4 * (1 - (alpha / 2) * y - (1 + (alpha / 2) * y) * std::exp(-alpha * y) / (y * y * y); }  

    double  h(double y, double alpha) { return (1 - std::exp(-alpha * y)) / y; }
    double dh(double y, double alpha) { return ((1 + alpha * y) * std::exp(-alpha * y) - 1) / (y * y); }

    double  S11(double q2) { return              g(m_epsilon * q2, m_f); }
    double dS11(double q2) { return m_epsilon * dg(m_epsilon * q2, m_f); }
    
    double  S22(double q2) { return  g(q2, 1 - m_f); }
    double dS22(double q2) { return dg(q2, 1 - m_f); }

    double  S12(double q2) { return h(m_epsilon * q2, m_f) * h(q2, 1 - m_f); } 
    double dS12(double q2) { return m_epsilon * dh(m_epsilon * q2, m_f) * h(q2, 1 - m_f) + h(m_epsilon * q2, m_f) * dh(q2, 1 - m_f); }

    double  S(double q2) { return  S11(q2) + 2 *  S12(q2) +  S22(q2); }
    double dS(double q2) { return dS11(q2) + 2 * dS12(q2) + dS22(q2); }

    double  W(double q2) { return S11(q2) * S22(q2) - S12(q2) * S12(q2); }
    double dW(double q2) { return dS11(q2) * S22(q2) + S11(q2) * dS22(q2) - 2 * S12(q2) * dS12(q2); }

    double computeQ0();

    double computeScaleFactor();

    // find quadratic potential coefficients, D(q)
    double quadraticCoeff(double q2) { return S(q2) / W(q2); }
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
    double quadraticCoeffDeriv(double q2) { return 2 * (q2 - 1); }

    // compute quadratic free-energy for given lattice spacings
    double fQuadB(FieldProvider &field, double* b);

};
#endif
