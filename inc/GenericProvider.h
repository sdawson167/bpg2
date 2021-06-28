#ifndef _GENERIC_GUARD
#define _GENERIC_GUARD

#include <stdexcept>
#include <string>
#include <math.h>
#include "Phase.h"
#include "FieldProvider.h"

// generic phase provider object - allows easy initialization of any phase!

class GenericPhaseProvider
{
private:
    const int m_phaseID;

    const double m_amplitude  = 1.0;
    const double m_avgDensity = 0.0;

    const double m_lamPeriod  = 2 * M_PI;
    const double m_hexPeriod  = 4 * M_PI / sqrt(3);
    const double m_gyrPeriod  = 15.0;
    const double m_bccPeriod  = 2 * sqrt(2) * M_PI;
    const double m_fccPeriod  = 2 * sqrt(3) * M_PI;
    const double m_a15Period  = 2 * sqrt(5) * M_PI;
    const double m_sigPeriodX = 8.88 * M_PI;
    const double m_sigPeriodZ = 4.54 * M_PI;
    const double m_c14PeriodX = 5 * M_PI;
    const double m_c14PeriodY = sqrt(3) * m_c14PeriodX;
    const double m_c14PeriodZ = 1.6 * m_c14PeriodX;
    const double m_c15Period  = 21;

    int stringToPhaseID(std::string phaseID) {
      if      (phaseID == "lam") {
        return 1;
      }
      else if (phaseID == "gyr") {
        return 2;
      }
      else if (phaseID == "hex") {
        return 3;
      }
      else if (phaseID == "bcc") {
        return 4;
      }
      else if (phaseID == "fcc") {
        return 5;
      }
      else if (phaseID == "a15") {
        return 6;
      }
      else if (phaseID == "sig") {
        return 7;
      }
      else if (phaseID == "c14") {
        return 8;
      }
      else if (phaseID == "c15") {
        return 9;
      }
      else {
        return 0;
      }
    }

public:
    GenericPhaseProvider(int phaseID)
      : m_phaseID{ phaseID }
    { }

    GenericPhaseProvider(std::string phaseID)
      : m_phaseID{ stringToPhaseID(phaseID) }
    { }

    Phase getPhase();

    int getPhaseID() { return m_phaseID; }

    /*
     * =======================================
     *     Method to initialize field
     * =======================================
     */
    FieldProvider* generateInitialCondition();

    /*
     * ======================================
     *  Method to reinitialize field values
     * ======================================
     */
    void resetCondition(FieldProvider &field);

};  // end class definition
#endif
