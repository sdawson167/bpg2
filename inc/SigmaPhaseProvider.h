#ifndef _SIGMA_GUARD
#define _SIGMA_GUARD

#include "Phase.h"
#include "FieldProvider.h"

/*
 * class object facilitates creation of field provider in Frank-Kasper A15 phase
 */

class SigmaPhaseProvider {
private:
    double m_periodX;
    double m_periodZ;
    double m_avDensity;
    double m_amplitude;
    const int m_dimension = 3;
    const int m_phaseID = 7;

    typedef std::vector<int>    intPoint;
    typedef std::tuple<intPoint, double> point;
    point makePoint(intPoint coords, double amp) { return point(coords, amp); }

public:
    SigmaPhaseProvider(double periodX, double periodZ, double avDensity, double amplitude)
      : m_periodX{ periodX },
        m_periodZ{ periodZ },
        m_avDensity{ avDensity },
        m_amplitude{ amplitude }
    { };

    /*
     * ======================================
     *          getters and setters
     * ======================================
     */
    double getXPeriod(){ return m_periodX; }

    double getZPeriod(){ return m_periodZ; }

    double getAverageDensity() { return m_avDensity; }

    double getAmplitude() { return m_amplitude; }

    Phase getPhase()
    {
      return Phase::sig;
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
     * =====================================
     *        initialize array values
     * =====================================
     */
    void populateDataArray(double* dqVec, fftw_complex* data, int numFieldElements, int* gridSizes);

};
#endif
