#include "C15PhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

FieldProvider C15PhaseProvider::generateInitialCondition(int gridSize) {

  // initialize field provider variables:

  const int N = gridSize;

  int* gridSizes = (int*) malloc(3 * sizeof(int));
  gridSizes[0] = N; gridSizes[1] = N; gridSizes[2] = N;

  const int numFieldElements = N * N * N;

  // grid spacing based on period (box size)
  double* dxVec = (double*) malloc(3 * sizeof(double));

  // initialize field in real space:
  const bool real = true;
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

// method to reset field to initial C15 condition
void C15PhaseProvider::resetCondition(FieldProvider &field)
{
  // verify phaseID
  const int phaseID = field.getPhaseID();
  if (phaseID != m_phaseID)
    throw std::runtime_error("Could not reset C15 phase - incorrect phase ID");

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

// initialize C15 phase in real space:
void C15PhaseProvider::populateDataArray(double* dxVec, fftw_complex* realData, int numFieldElements, int* gridSizes)
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
  double coords[34][3] = {
    {0, 0, 0},
    {1.00000, 1.00000, 0.00000},
		{1.00000, 0.50000, 0.50000},
		{1.00000, 0.00000, 1.00000},
		{1.00000, 1.00000, 1.00000},
		{0.75000, 0.75000, 0.25000},
		{0.75000, 0.25000, 0.75000},
		{0.50000, 0.50000, 0.00000},
		{0.50000, 0.00000, 0.50000},
		{0.50000, 1.00000, 0.50000},
		{0.50000, 0.50000, 1.00000},
		{0.25000, 0.25000, 0.25000},
		{0.25000, 0.75000, 0.75000},
		{0.00000, 0.00000, 0.00000},
		{0.00000, 1.00000, 0.00000},
		{0.00000, 0.50000, 0.50000},
		{0.00000, 0.00000, 1.00000},
		{0.00000, 1.00000, 1.00000},
		{0.87500, 0.37500, 0.12500},
		{0.87500, 0.12500, 0.37500},
		{0.87500, 0.87500, 0.62500},
		{0.87500, 0.62500, 0.87500},
		{0.62500, 0.87500, 0.87500},
		{0.62500, 0.62500, 0.62500},
		{0.62500, 0.37500, 0.37500},
		{0.62500, 0.12500, 0.12500},
		{0.37500, 0.87500, 0.12500},
		{0.37500, 0.62500, 0.37500},
		{0.37500, 0.37500, 0.62500},
		{0.37500, 0.12500, 0.87500},
		{0.12500, 0.37500, 0.87500},
		{0.12500, 0.12500, 0.62500},
		{0.12500, 0.87500, 0.37500},
		{0.12500, 0.62500, 0.12500}
  };

  double r0sqrd = m_period / 8.0;

  double sum = 0.0;

  for (int k = 0; k < Nz; k++) 
  for (int j = 0; j < Ny; j++) 
  for (int i = 0; i < Nx; i++) {

    int index = i + Nx * (j + Ny * k);

    double rx = i * dx;
    double ry = j * dx;
    double rz = k * dx;

    // loop over list of peaks:
    for (int m = 0; m < 34; m++) {
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
  }

  double avg = sum / numFieldElements;

  for (int index = 0; index < numFieldElements; index++) {
    realData[index][0] = realData[index][0] - avg;
  }

} // end of populateDataArray method
