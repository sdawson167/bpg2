#include "GenericProvider.h"

#include "LamellarPhaseProvider.h"
#include "GyroidPhaseProvider.h"
#include "CylindricalHexagonalPhaseProvider.h"
#include "BccPhaseProvider.h"
#include "FccPhaseProvider.h"
#include "A15PhaseProvider.h"
#include "SigmaPhaseProvider.h"
#include "C14PhaseProvider.h"
#include "C15PhaseProvider.h"


FieldProvider* GenericPhaseProvider::generateInitialCondition()
{
  switch(m_phaseID) {
    case 1: {   // lam
      LamellarPhaseProvider lamProvider{ m_lamPeriod, m_avgDensity, m_amplitude };
      FieldProvider* initialFieldProvider = new FieldProvider(lamProvider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 2: {   // gyr
      GyroidPhaseProvider gyrProvider{ m_gyrPeriod, m_avgDensity, m_amplitude };
      FieldProvider* initialFieldProvider = new FieldProvider(gyrProvider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 3: {   // hex
      CylindricalHexagonalPhaseProvider hexProvider{ m_hexPeriod, m_avgDensity, m_amplitude };
      FieldProvider* initialFieldProvider = new FieldProvider(hexProvider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 4: {   // bcc
      BccPhaseProvider bccProvider{ m_bccPeriod, m_avgDensity, m_amplitude};
      FieldProvider* initialFieldProvider = new FieldProvider(bccProvider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 5: {   // fcc
      FccPhaseProvider fccProvider{ m_fccPeriod, m_avgDensity, m_amplitude};
      FieldProvider* initialFieldProvider = new FieldProvider(fccProvider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 6: {   // a15
      A15PhaseProvider a15Provider{ m_a15Period, m_avgDensity, m_amplitude};
      FieldProvider* initialFieldProvider = new FieldProvider(a15Provider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 7: {
      SigmaPhaseProvider sigProvider{ m_sigPeriodX, m_sigPeriodZ, m_avgDensity, m_amplitude};
      FieldProvider* initialFieldProvider = new FieldProvider(sigProvider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 8: {
      C14PhaseProvider c14Provider{ m_c14PeriodX, m_c14PeriodY, m_c14PeriodZ, m_avgDensity, m_amplitude};
      FieldProvider* initialFieldProvider = new FieldProvider(c14Provider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    case 9: {
      C15PhaseProvider c15Provider{ m_c15Period, m_avgDensity, m_amplitude};
      FieldProvider* initialFieldProvider = new FieldProvider(c15Provider.generateInitialCondition(64));
      return initialFieldProvider;
    }
    default: {  // dis
      int* gridSize = (int*) malloc(sizeof(int));
      gridSize[0] = 1;      

      double* dr = (double*) malloc(sizeof(double));
      dr[0] = 0.0;

      fftw_complex* data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex));
      FieldProvider* initialFieldProvider = 
      new FieldProvider{
        data, 
        1,
        gridSize,
        dr,
        false,
        0};
      fftw_free(data);
      
      return initialFieldProvider;
    }
  } // end switch
} // end generateInitialCondition

void GenericPhaseProvider::resetCondition(FieldProvider &field)
{
  switch(m_phaseID) {
    case 1: {
      LamellarPhaseProvider lamProvider{ m_lamPeriod, m_avgDensity, m_amplitude };
      lamProvider.resetCondition(field);
      break;
    }
    case 2: {
      GyroidPhaseProvider gyrProvider{ m_gyrPeriod, m_avgDensity, m_amplitude };
      gyrProvider.resetCondition(field);
      break;
    }
    case 3: {
      CylindricalHexagonalPhaseProvider hexProvider{ m_hexPeriod, m_avgDensity, m_amplitude };
      hexProvider.resetCondition(field);
      break;
    }
    case 4: {
      BccPhaseProvider bccProvider{ m_bccPeriod, m_avgDensity, m_amplitude };
      bccProvider.resetCondition(field);
      break;
    }
    case 5: {
      FccPhaseProvider fccProvider{ m_fccPeriod, m_avgDensity, m_amplitude };
      fccProvider.resetCondition(field);
      break;
    }
    case 6: {
      A15PhaseProvider a15Provider{ m_a15Period, m_avgDensity, m_amplitude };
      a15Provider.resetCondition(field);
      break;
    }
    case 7: {
      SigmaPhaseProvider sigProvider{ m_sigPeriodX, m_sigPeriodZ, m_avgDensity, m_amplitude };
      sigProvider.resetCondition(field);
      break;
    }
    case 8: {
      C14PhaseProvider c14Provider{ m_c14PeriodX, m_c14PeriodY, m_c14PeriodZ, m_avgDensity, m_amplitude };
      c14Provider.resetCondition(field);
      break;
    }
    case 9: {
      C15PhaseProvider c15Provider{ m_c15Period, m_avgDensity, m_amplitude };
      c15Provider.resetCondition(field);
      break;
    }
    default: {
      int phaseID = field.getPhaseID();
      if (phaseID != m_phaseID)
	throw std::runtime_error("Unable to reset DIS phase - wrong phase ID");
      
      int N = field.getNumFieldElements();
      fftw_complex* realFieldData = field.getRealDataPointer();
      fftw_complex* cplxFieldData = field.getCplxDataPointer();

      for (int i = 0; i < N; i++) {
	realFieldData[i][0] = 0.0;
	realFieldData[i][1] = 0.0;

	cplxFieldData[i][0] = 0.0;
	cplxFieldData[i][1] = 0.0;
      }
    }
  } // end switch
} // end resetCondition
