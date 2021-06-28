#ifndef _HEX_GUARD
#define _HEX_GUARD

#include "Phase.h"
#include "FieldProvider.h"

/*
 * class object facilitates creation of field provider in cylindrical hexagonal phase
 */

class CylindricalHexagonalPhaseProvider {
private:
    double m_period;
    double m_avDensity;
    double m_amplitude;
    const int m_dimension = 2;
    const int m_phaseID = 3;

    // range of tau values for which phase will "definitely" converge
    const double m_tauMin = -0.4;
    const double m_tauMax = -0.01;

    // range of gamma values for which phase will "definitely" converge
    const double m_gammaMin = 0.0;
    const double m_gammaMax = 1.8;

    typedef std::vector<int>    intPoint;                                      	
    typedef std::tuple<intPoint, double> point;
    point makePoint(intPoint coords, double amp) { return point(coords, amp); }

public:
    CylindricalHexagonalPhaseProvider(double period, double avDensity, double amplitude)
    : m_period{ period },
      m_avDensity{ avDensity },
      m_amplitude{ amplitude }
    {};

    /*
     * ======================================
     *          getters and setters
     * ======================================
     */
    double getPeriod(){ return m_period; }

    double getAverageDensity() { return m_avDensity; }

    double getAmplitude() { return m_amplitude; }

    Phase getPhase()
    {
      return Phase::hex;
    }

    double getTauMin() { return m_tauMin; }       
    double getTauMax() { return m_tauMax; }
    
    double getGammaMin() { return m_gammaMin; }
    double getGammaMax() { return m_gammaMax; }

    /*
     * ======================================
     *      initialize field provider
     * ======================================
     */
    FieldProvider generateInitialCondition(int gridSize);

    /*
     * ======================================
     *	       reset field provider
     * ======================================
     */
    void resetCondition(FieldProvider &field);

    /*
     * =====================================
     *       initialize array values 
     * =====================================
     */
    void populateDataArray(double* dqVec, fftw_complex* data, int numFieldElements, int* gridSizes);
};
#endif
