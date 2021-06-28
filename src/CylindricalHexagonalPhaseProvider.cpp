#include "CylindricalHexagonalPhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

FieldProvider CylindricalHexagonalPhaseProvider::generateInitialCondition(int gridSize) {

  // initialize field provider variables:

  // Ny = 2 * Nx for hexagonal phase
  const int N = gridSize;

  int* gridSizes = (int*) malloc(2 * sizeof(int));
  gridSizes[0] = 2 * N; gridSizes[1] = N;

  const int numFieldElements = 2 * N * N;
  
  double* dq = (double*) malloc(2 * sizeof(double));

  // initialize field in complex space:
  const bool real = false;

  // initialize data -
  fftw_complex* data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);
  populateDataArray(dq, data, numFieldElements, gridSizes);

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
} // end of generateInitialCondition method

void CylindricalHexagonalPhaseProvider::resetCondition(FieldProvider &field)
{
  // verify phaseID
  const int phaseID = field.getPhaseID();
  if (phaseID != m_phaseID)
    throw std::runtime_error("Could not reset HEX phase - incorrect phase ID");

  // unpack field provider
  int* gridSizes         = field.getGridSizes();
  int  numFieldElements  = field.getNumFieldElements();
  double* dqVec		 = field.getDq();
  fftw_complex* cplxData = field.getCplxDataPointer();
  
  // set values of cplxData
  populateDataArray(dqVec, cplxData, numFieldElements, gridSizes);
  field.updateDx();

  // update real data
  field.transformC2R();
} // end of resetCondition method

void CylindricalHexagonalPhaseProvider::populateDataArray(double* dqVec, fftw_complex* data, int numFieldElements, int* gridSizes)
{
  // grid spacing based on period (box size)
  const double dqx = 2 * M_PI / m_period;
  const double dqy = 2 * M_PI / (2 * sqrt(3) * m_period);
  dqVec[0] = dqy; dqVec[1] = dqx;

  // set all array values to zero
  for (int index = 0; index < numFieldElements; index++) {
    data[index][0] = 0.0;
    data[index][1] = 0.0;
  }
  
  // set of fourier peaks
  std::vector<point> initVals;
  initVals.push_back( makePoint({ 1,  2}, 1.0/6 ));
  initVals.push_back( makePoint({-1,  2}, 1.0/6 ));
  initVals.push_back( makePoint({ 1, -2}, 1.0/6 ));
  initVals.push_back( makePoint({-1, -2}, 1.0/6 ));
  initVals.push_back( makePoint({ 0,  4}, 1.0/6 ));
  initVals.push_back( makePoint({ 0, -4}, 1.0/6 ));

  // initialize non-zero array values using peaks
  for (std::vector<point>::iterator it = initVals.begin(); it != initVals.end(); it++) {
    std::vector<int> k = std::get<0>(*it);	  
    double	     u = std::get<1>(*it);

    int ky = k[1] < 0 ? k[1] + gridSizes[0] : k[1];
    int kx = k[0] < 0 ? k[0] + gridSizes[1] : k[0];

    int index = kx + (gridSizes[1] * ky);

    data[index][0] = u;
  }

} // end populateDataArray method
