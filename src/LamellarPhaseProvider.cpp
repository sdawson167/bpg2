#include "LamellarPhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

FieldProvider LamellarPhaseProvider::generateInitialCondition(int gridSize) {

  // initialize field provider variables:
  int* gridSizes = (int*) malloc(sizeof(int));
  gridSizes[0] = gridSize;

  // grid spacing based on period (box size)
  double* dq = (double*) malloc(sizeof(double));

  // initialize field in complex space:
  const bool real = false;

  // initialize data
  fftw_complex* data;
  data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * gridSize);
  populateDataArray(dq, data, gridSize, gridSizes);

  // create field provider object
  FieldProvider initialCondition{
    data,
    m_dimension,
    gridSizes,
    dq,
    real,
    m_phaseID};

  fftw_free(data);

  return initialCondition;
}

void LamellarPhaseProvider::resetCondition(FieldProvider &field)
{
  // verify phaseID
  const int phaseID = field.getPhaseID();
  if (phaseID != m_phaseID)
    throw std::runtime_error("Could not reset LAM phase - incorrect phase ID");

  // unpack field provider
  int* gridSizes         = field.getGridSizes();
  int  numFieldElements  = field.getNumFieldElements();
  double* dqVec		 = field.getDq();
  fftw_complex* cplxData = field.getCplxDataPointer();
  
  // set values of cplxData
  populateDataArray(dqVec, cplxData, numFieldElements, gridSizes);
  field.updateDx();

  // update cplx data
  field.transformC2R();
}

void LamellarPhaseProvider::populateDataArray(double* dqVec, fftw_complex* data, int numFieldElements, int* gridSizes)
{
  // update grid spacing
  const double gridSpacing = 2 * M_PI / m_period;
  dqVec[0] = gridSpacing;

  // set all array values to zero
  for (int index = 0; index < numFieldElements; index++) {
    data[index][0] = 0.0;
    data[index][1] = 0.0;
  }

  // initialize lamellar phase - 1D cosine function
  data[0][0]                = 0.0;
  data[1][0]                = 0.5;
  data[numFieldElements-1][0] = 0.5; 
}
