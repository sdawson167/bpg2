#include "SigmaPhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

FieldProvider SigmaPhaseProvider::generateInitialCondition(int gridSize) {

  // initialize field provider variables:

  const int N = gridSize;

  int* gridSizes = (int*) malloc(3 * sizeof(int));
  gridSizes[0] = N; gridSizes[1] = 2 * N; gridSizes[2] = 2 * N;

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

// method to reset field to initial sigma condition
void SigmaPhaseProvider::resetCondition(FieldProvider &field)
{
  // verify phaseID
  const int phaseID = field.getPhaseID();
  if (phaseID != m_phaseID)
    throw std::runtime_error("Could not reset sigma phase - incorrect phase ID");

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

// initialize A15 phase in real space:
void SigmaPhaseProvider::populateDataArray(double* dxVec, fftw_complex* realData, int numFieldElements, int* gridSizes)
{
  int Nx = gridSizes[2];
  int Ny = gridSizes[1];
  int Nz = gridSizes[0];

  // reset grid spacing
  const double dx = m_periodX / Nx;
  const double dz = m_periodZ / Nz;
  dxVec[0] = dz;  dxVec[1] = dx;  dxVec[2] = dx;

  // initialize values to zero
  for (int index = 0; index < numFieldElements; index++) {
    realData[index][0] = 0.0;
    realData[index][1] = 0.0;
  }

  // locations of peaks
  double coords[47][3] = {
    {0, 0, 0},
    {1, 0, 1},
    {0, 1, 0},
    {0, 0, 1},
    {1, 1, 0},
    {1, 0, 0},
    {0, 1, 1},
    {1, 1, 1},
    {0.5, 0.5, 0.5},
    {0.36840, 0.96320, 0.50000},
    {0.53680, 0.86840, 0.00000},
    {0.53680, 0.86840, 1.00000},
    {0.86840, 0.53680, 0.00000},
    {0.96320, 0.36840, 0.50000},
    {0.86840, 0.53680, 1.00000},
    {0.03680, 0.63160, 0.50000},
  	{0.13160, 0.46320, 1.00000},
  	{0.13160, 0.46320, 0.00000},
  	{0.46320, 0.13160, 0.00000},
  	{0.63160, 0.03680, 0.50000},
  	{0.46320, 0.13160, 1.00000},
  	{0.31770, 0.68230, 0.75240},
  	{0.31770, 0.68230, 0.24760},
  	{0.81770, 0.81770, 0.25240},
  	{0.81770, 0.81770, 0.74760},
  	{0.68230, 0.31770, 0.75240},
  	{0.68230, 0.31770, 0.24760},
  	{0.18230, 0.18230, 0.74760},
  	{0.18230, 0.18230, 0.25240},
  	{0.10190, 0.89810, 0.50000},
  	{0.60190, 0.60190, 1.00000},
  	{0.60190, 0.60190, 0.00000},
  	{0.39810, 0.39810, 1.00000},
  	{0.39810, 0.39810, 0.00000},
  	{0.89810, 0.10190, 0.50000},
  	{0.06530, 0.73760, 0.00000},
  	{0.06530, 0.73760, 1.00000},
  	{0.26240, 0.93470, 0.00000},
  	{0.26240, 0.93470, 1.00000},
  	{0.56530, 0.76240, 0.50000},
  	{0.23760, 0.43470, 0.50000},
  	{0.76240, 0.56530, 0.50000},
  	{0.43470, 0.23760, 0.50000},
  	{0.73760, 0.06530, 0.00000},
  	{0.73760, 0.06530, 1.00000},
  	{0.93470, 0.26240, 0.00000},
  	{0.93470, 0.26240, 1.00000}
  };

  double r0sqrd = m_periodX / 4.0;

  double sum = 0.0;
  double fMin = 0.0 - 1e-8;
  double fMax = 0.0 + 1e-8;

  for (int k = 0; k < Nz; k++) 
  for (int j = 0; j < Ny; j++) 
  for (int i = 0; i < Nx; i++) {

    int index = i + Nx * (j + Ny * k);

    double rx = i * dx;
    double ry = j * dx;
    double rz = k * dz;

    // loop over list of peaks:
    for (int m = 0; m < 47; m++) {
      double rx0 = coords[m][0] * m_periodX;
      double ry0 = coords[m][1] * m_periodX;
      double rz0 = coords[m][2] * m_periodZ;

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
