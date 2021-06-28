#ifndef _GYROID_GUARD
#define _GYROID_GUARD

#include "Phase.h"
#include "FieldProvider.h"

/*
 * class object facilitates creation of field provider in gyroid phase
 */

class GyroidPhaseProvider {
private:
    double m_period;
    double m_avDensity;
    double m_amplitude;
    const int m_dimension = 3;
    const int m_phaseID = 2;

    // range of tau values for which newly initialized phase will "definitely" converge
    const double m_tauMin = -0.4;
    const double m_tauMax = -0.15;

    // range of gamma values for which newly initialized phase will "definitely" converge
    const double m_gammaMin = 0.1;
    const double m_gammaMax = 1.0;

public:
    GyroidPhaseProvider(double period, double avDensity, double amplitude)
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
      return Phase::gyr;
    }

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
     * ======================================
     *      initialize array elements
     * ======================================
     */
    void populateDataArray(double* dxVec, fftw_complex* data, int numFieldElements, int* gridSizes);
};
#endif
