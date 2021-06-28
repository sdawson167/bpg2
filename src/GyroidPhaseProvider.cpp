#include "GyroidPhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

typedef std::vector<int>    intPoint;
typedef std::tuple<intPoint, double> point;
point makePoint(intPoint coords, double amp) { return point(coords, amp); }

FieldProvider GyroidPhaseProvider::generateInitialCondition(int gridSize) {

  // initialize field provider variables:
  int* gridSizes = (int*) malloc(3 * sizeof(int));
  gridSizes[0] = gridSize;  gridSizes[1] = gridSize;  gridSizes[2] = gridSize;

  const int numFieldElements = gridSize * gridSize * gridSize;

  // grid spacing based on period (box size)
  double* dxVec = (double*) malloc(3 * sizeof(double));

  // initialize field in real space:
  const bool real = true;

  // initialize data -
  fftw_complex* data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);
  populateDataArray(dxVec, data, numFieldElements, gridSizes);
  
  // create field provider object
  FieldProvider initialCondition{
    data,
    m_dimension,
    gridSizes,
    dxVec,
    real,
    m_phaseID};

  fftw_free(data);

  return initialCondition;
}

void GyroidPhaseProvider::resetCondition(FieldProvider &field)
{
  // verify phaseID
  const int phaseID = field.getPhaseID();
  if (phaseID != m_phaseID)
    throw std::runtime_error("Could not reset GYR phase - incorrect phase ID");

  // unpack field provider
  int* gridSizes         = field.getGridSizes();
  int  numFieldElements  = field.getNumFieldElements();
  double* dxVec		 = field.getDx();
  fftw_complex* realData = field.getRealDataPointer();
  
  // set values of cplxData
  populateDataArray(dxVec, realData, numFieldElements, gridSizes);
  field.updateDq();

  // update cplx data
  field.transformR2C();

} // end resetCondition method

void GyroidPhaseProvider::populateDataArray(double* dxVec, fftw_complex* data, int numFieldElements, int* gridSizes)
{
  const double dx = m_period / gridSizes[0];
  dxVec[0] = dx; dxVec[1] = dx; dxVec[2] = dx;

  // set all array values to zero
  for (int index = 0; index < numFieldElements; index++) {
    data[index][0] = 0.0;
    data[index][1] = 0.0;
  }

  // compute real gyroid array values
  double w = 2 * M_PI / m_period;
  double sum = 0.0;
  for (int k = 0, index = 0; k < gridSizes[0]; k++) {
    double z = k * dx;

    for (int j = 0; j < gridSizes[1]; j++) {
      double y = j * dx;

      for (int i = 0; i < gridSizes[2]; i++, index++) {
        double x = i * dx;

        double val = sin(w * x) * cos(w * y) +
                     sin(w * y) * cos(w * z) +
                     sin(w * z) * cos(w * x);

        if (val > 1 || val < -1)
          data[index][0] = m_amplitude;
        else
          data[index][0] = 0.0;

        sum += data[index][0];
      }
    }
  }

  // shift everything by avg. density
  double avDensity = sum / numFieldElements;
  for(int index = 0; index < numFieldElements; index++)
    data[index][0] += (m_avDensity - avDensity);

} // end populate array method
