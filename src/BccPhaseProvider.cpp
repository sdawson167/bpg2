#include "BccPhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

FieldProvider BccPhaseProvider::generateInitialCondition(int gridSize) {

  // initialize field provider variables:

  const int N = gridSize;

  int* gridSizes = (int*) malloc(3 * sizeof(int));
  gridSizes[0] = N; gridSizes[1] = N; gridSizes[2] = N;

  const int numFieldElements = N * N * N;

  // grid spacing based on period (box size)
  double* dqVec = (double*) malloc(3 * sizeof(double));

  // initialize field in real space:
  const bool real = true;
  fftw_complex* data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numFieldElements);
  populateDataArray(dqVec, data, numFieldElements, gridSizes);
  
  // create field provider object
  FieldProvider initialCondition{
    data,
    m_dimension,
    gridSizes,
    dqVec,
    real,
    m_phaseID};

  fftw_free(data);

  return initialCondition;
}

// method to reset field to initial BCC condition
void BccPhaseProvider::resetCondition(FieldProvider &field)
{
  // verify phaseID
  const int phaseID = field.getPhaseID();
  if (phaseID != m_phaseID)
    throw std::runtime_error("Could not reset BCC phase - incorrect phase ID");

  // unpack field provider
  int* gridSizes         = field.getGridSizes();
  int  numFieldElements  = field.getNumFieldElements();
  double* dxVec    = field.getDx();
  fftw_complex* realData = field.getRealDataPointer();
  
  // set values of realData
  populateDataArray(dxVec, realData, numFieldElements, gridSizes);
  
  // update complex grid spacing
  field.updateDq();

  // update complex data
  field.transformR2C();
}

// initialize BCC phase in real space:
void BccPhaseProvider::populateDataArray(double* dxVec, fftw_complex* realData, int numFieldElements, int* gridSizes)
{
  int Nx = gridSizes[2];
  int Ny = gridSizes[1];
  int Nz = gridSizes[0];

  // reset grid spacing
  const double dx = m_period / Nx;
  dxVec[0] = dx;  dxVec[1] = dx;  dxVec[2] = dx;

  // initialize values to zero
  for (int index = 0; index < numFieldElements; index++) {
    realData[index][0] = 0.0;
    realData[index][1] = 0.0;
  }

  // locations of peaks
  double coords[9][3] = {{0.5, 0.5, 0.5},
                         {0, 0, 0},
                         {1, 0, 0},
                         {0, 1, 0},
                         {0, 0, 1},
                         {1, 1, 0},
                         {1, 0, 1},
                         {0, 1, 1},
                         {1, 1, 1}};

  double r0sqrd = m_period;

  double sum = 0.0;
  double fMin = 0.0 - 1e-8;
  double fMax = 0.0 + 1e-8;

  for (int k = 0; k < Nz; k++) 
  for (int j = 0; j < Ny; j++) 
  for (int i = 0; i < Nx; i++) {

    int index = i + Nx * (j + Ny * k);

    double rx = i * dx;
    double ry = j * dx;
    double rz = k * dx;

    // loop over list of peaks:
    for (int m = 0; m < 9; m++) {
      double rx0 = coords[m][0] * m_period;
      double ry0 = coords[m][1] * m_period;
      double rz0 = coords[m][2] * m_period;

      double rSqrd = (rx - rx0) * (rx - rx0) + (ry - ry0) * (ry - ry0) + (rz - rz0) * (rz - rz0);

      if (rSqrd < r0sqrd) {
        double x = rSqrd / r0sqrd;
        if (x < sqrt(1 + 0.01)) 
          realData[index][0] = exp(-1.0 / (1.0 - x * x) + 1.0);
      }
    }

    sum += realData[index][0];
    if (realData[index][0] < fMin)
      fMin = realData[index][0];
    if (realData[index][0] > fMax)
      fMax = realData[index][0];
  }

  double avg = sum / numFieldElements;
  double hlfRng = (fMax - fMin) / 2.0;

  for (int index = 0; index < numFieldElements; index++) {
    realData[index][0] -= avg;
    realData[index][0] /= hlfRng;
  }

} // end of populateDataArray method

