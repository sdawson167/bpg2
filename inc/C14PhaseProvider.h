#ifndef _C14_GUARD
#define _C14_GUARD

#include "Phase.h"
#include "FieldProvider.h"

/*
 * class object facilitates creation of field provider in Frank-Kasper A15 phase
 */

class C14PhaseProvider {
private:
    double m_periodX;
    double m_periodY;
    double m_periodZ;
    double m_avDensity;
    double m_amplitude;
    const int m_dimension = 3;
    const int m_phaseID = 8;

    /*
    double m_tauMin;
    double m_tauMax;
    double m_tau0;
    double m_gammaMin;
    double m_gammaMax;
    double m_gamma0;
    */

    typedef std::vector<int>    intPoint;
    typedef std::tuple<intPoint, double> point;
    point makePoint(intPoint coords, double amp) { return point(coords, amp); }

public:
    C14PhaseProvider(double periodX, double periodY, double periodZ, double avDensity, double amplitude)
      : m_periodX{ periodX },
        m_periodY{ periodY },
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

    double getYPeriod(){ return m_periodY; }

    double getZPeriod(){ return m_periodZ; }

    double getAverageDensity() { return m_avDensity; }

    double getAmplitude() { return m_amplitude; }

    Phase getPhase()
    {
      return Phase::c14;
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
    void populateDataArray(double* dxVec, fftw_complex* data, int numFieldElements, int* gridSizes);

};
#endif
