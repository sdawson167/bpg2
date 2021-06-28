#include "C14PhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

FieldProvider C14PhaseProvider::generateInitialCondition(int gridSize) {

  // initialize field provider variables:

  const int N = gridSize;

  int* gridSizes = (int*) malloc(3 * sizeof(int));
  gridSizes[0] = 2 * N; gridSizes[1] = 2 * N; gridSizes[2] = N;

  const int numFieldElements = 4 * N * N * N;

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

// method to reset field to initial C14 condition
void C14PhaseProvider::resetCondition(FieldProvider &field)
{
  // verify phaseID
  const int phaseID = field.getPhaseID();
  if (phaseID != m_phaseID)
    throw std::runtime_error("Could not reset C14 phase - incorrect phase ID");

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

// initialize C14 phase in real space:
void C14PhaseProvider::populateDataArray(double* dxVec, fftw_complex* realData, int numFieldElements, int* gridSizes)
{
  int Nx = gridSizes[2];
  int Ny = gridSizes[1];
  int Nz = gridSizes[0];

  // reset grid spacing
  const double dx = m_periodX / Nx;
  const double dy = m_periodY / Ny;
  const double dz = m_periodZ / Nz;
  dxVec[0] = dz;  dxVec[1] = dy;  dxVec[2] = dx;

  // initialize values to zero
  for (int index = 0; index < numFieldElements; index++) {
    realData[index][0] = 0.0;
    realData[index][1] = 0.0;
  }

  // locations of peaks
  double coords[41][3] = {
    {0, 0, 0},
    {0, 0, 1},
		{0, 1, 0},
		{0, 1, 1},
		{1, 1, 0},
		{1, 1, 1},
		{1, 0, 0},
		{1, 0, 1},
		{0.5, 0.5, 0.5},
		{0.5, 0.5, 0},
		{0.5, 0.5, 1},
		{1, 1, 0.5},
		{1, 0, 0.5},
		{0, 0, 0.5},
		{0, 1, 0.5},
		{0.74575, 0.91525, 0.75},
		{0.25425, 0.91525, 0.75},
		{1, 0.8305, 0.25},
		{0, 0.8305, 0.25},
		{0.5, 0.6695, 0.75},
		{0.75425, 0.58475, 0.25},
		{0.24575, 0.58475, 0.25},
		{0.75425, 0.41525, 0.75},
		{0.24575, 0.41525, 0.75},
		{0.5, 0.3305, 0.25},
		{1, 0.1695, 0.75},
		{0.74575, 0.08475, 0.25},
		{0.25425, 0.08475, 0.25},
		{0, 0.1695, 0.75},
		{0.5, 0.83334, 0.438},
		{0.5, 0.83334, 0.062},
		{1, 0.66666, 0.562},
		{1, 0.66666, 0.938},
		{0, 0.66666, 0.562},
		{0, 0.66666, 0.938},
		{1, 0.33334, 0.062},
		{1, 0.33334, 0.438},
		{0.5, 0.16666, 0.938},
		{0.5, 0.16666, 0.562},
		{0, 0.33334, 0.438},
		{0, 0.33334, 0.062}
  };

  double r0sqrd = m_periodX / 6.0;

  double sum = 0.0;

  for (int k = 0; k < Nz; k++) 
  for (int j = 0; j < Ny; j++) 
  for (int i = 0; i < Nx; i++) {

    int index = i + Nx * (j + Ny * k);

    double rx = i * dx;
    double ry = j * dy;
    double rz = k * dz;

    // loop over list of peaks:
    for (int m = 0; m < 41; m++) {
      double rx0 = coords[m][0] * m_periodX;
      double ry0 = coords[m][1] * m_periodY;
      double rz0 = coords[m][2] * m_periodZ;

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

