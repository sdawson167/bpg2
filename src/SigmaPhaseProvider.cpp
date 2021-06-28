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

  double r0sqrd = m_periodX / 8.0;

  double sum = 0.0;

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
  }

  double avg = sum / numFieldElements;

  for (int index = 0; index < numFieldElements; index++) {
    realData[index][0] = realData[index][0] - avg;
  }

} // end of populateDataArray method

/*
void SigmaPhaseProvider::populateDataArray(double* dqVec, fftw_complex* data, int numFieldElements, int* gridSizes)
{
  // update grid spacing
  const double dqx = 2 * M_PI/m_periodX;
  const double dqz = 2 * M_PI/m_periodZ;
  dqVec[0] = dqz; dqVec[1] = dqx; dqVec[2] = dqx;

  // set all array values to zero
  for (int index = 0; index < numFieldElements; index++) {
    data[index][0] = 0.0;
    data[index][1] = 0.0;
  }

  // set of fourier peaks 
  typedef std::vector<int>    intPoint;
  typedef std::tuple<intPoint, double> point;
  std::vector<point> initVals;
  initVals.push_back(makePoint( {   0,   0,  0 },  m_avDensity));
  initVals.push_back(makePoint( {   4,   0,  0 }, -0.0167031));
  initVals.push_back(makePoint( { 124,   0,  0 }, -0.0167031));
  initVals.push_back(makePoint( {   3,   1,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   4,   1,  0 },  0.1859610));
  initVals.push_back(makePoint( {   5,   1,  0 },  0.0155903));
  initVals.push_back(makePoint( { 123,   1,  0 },  0.0155903));
  initVals.push_back(makePoint( { 124,   1,  0 }, -0.1859610));
  initVals.push_back(makePoint( { 125,   1,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   3,   2,  0 },  0.0133784));
  initVals.push_back(makePoint( { 125,   2,  0 }, -0.0133784));
  initVals.push_back(makePoint( {   1,   3,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   2,   3,  0 },  0.0133784));
  initVals.push_back(makePoint( {   3,   3,  0 },  0.1906830));
  initVals.push_back(makePoint( {   5,   3,  0 },  0.0105081));
  initVals.push_back(makePoint( { 123,   3,  0 },  0.0105081));
  initVals.push_back(makePoint( { 125,   3,  0 },  0.1906830));
  initVals.push_back(makePoint( { 126,   3,  0 }, -0.0133784));
  initVals.push_back(makePoint( { 127,   3,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   0,   4,  0 }, -0.0167031));
  initVals.push_back(makePoint( {   1,   4,  0 },  0.1859610));
  initVals.push_back(makePoint( { 127,   4,  0 }, -0.1859610));
  initVals.push_back(makePoint( {   1,   5,  0 },  0.0155903));
  initVals.push_back(makePoint( {   3,   5,  0 },  0.0105081));
  initVals.push_back(makePoint( {   5,   5,  0 },  0.0107489));
  initVals.push_back(makePoint( { 123,   5,  0 },  0.0107489));
  initVals.push_back(makePoint( { 125,   5,  0 },  0.0105081));
  initVals.push_back(makePoint( { 127,   5,  0 },  0.0155903));
  initVals.push_back(makePoint( {   1, 123,  0 },  0.0155903));
  initVals.push_back(makePoint( {   3, 123,  0 },  0.0105081));
  initVals.push_back(makePoint( {   5, 123,  0 },  0.0107489));
  initVals.push_back(makePoint( { 123, 123,  0 },  0.0107489));
  initVals.push_back(makePoint( { 125, 123,  0 },  0.0105081));
  initVals.push_back(makePoint( { 127, 123,  0 },  0.0155903));
  initVals.push_back(makePoint( {   0, 124,  0 }, -0.0167031));
  initVals.push_back(makePoint( {   1, 124,  0 }, -0.1859610));
  initVals.push_back(makePoint( { 127, 124,  0 },  0.1859610));
  initVals.push_back(makePoint( {   1, 125,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   2, 125,  0 }, -0.0133784));
  initVals.push_back(makePoint( {   3, 125,  0 },  0.1906830));
  initVals.push_back(makePoint( {   5, 125,  0 },  0.0105081));
  initVals.push_back(makePoint( { 123, 125,  0 },  0.0105081));
  initVals.push_back(makePoint( { 125, 125,  0 },  0.1906830));
  initVals.push_back(makePoint( { 126, 125,  0 },  0.0133784));
  initVals.push_back(makePoint( { 127, 125,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   3, 126,  0 }, -0.0133784));
  initVals.push_back(makePoint( { 125, 126,  0 },  0.0133784));
  initVals.push_back(makePoint( {   3, 127,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   4, 127,  0 }, -0.1859610));
  initVals.push_back(makePoint( {   5, 127,  0 },  0.0155903));
  initVals.push_back(makePoint( { 123, 127,  0 },  0.0155903));
  initVals.push_back(makePoint( { 124, 127,  0 },  0.1859610));
  initVals.push_back(makePoint( { 125, 127,  0 }, -0.0152148));
  initVals.push_back(makePoint( {   3,   0,  1 }, -0.0127882));
  initVals.push_back(makePoint( { 125,   0,  1 }, -0.0127882));
  initVals.push_back(makePoint( {   3,   1,  1 }, -0.0340823));
  initVals.push_back(makePoint( {   4,   1,  1 },  0.1662350));
  initVals.push_back(makePoint( {   5,   1,  1 },  0.0224927));
  initVals.push_back(makePoint( { 123,   1,  1 }, -0.0224927));
  initVals.push_back(makePoint( { 124,   1,  1 },  0.1662350));
  initVals.push_back(makePoint( { 125,   1,  1 },  0.0340823));
  initVals.push_back(makePoint( {   2,   2,  1 }, -0.0142726));
  initVals.push_back(makePoint( {   3,   2,  1 }, -0.0145262));
  initVals.push_back(makePoint( {   4,   2,  1 },  0.0103089));
  initVals.push_back(makePoint( {   5,   2,  1 },  0.0102342));
  initVals.push_back(makePoint( { 123,   2,  1 },  0.0102342));
  initVals.push_back(makePoint( { 124,   2,  1 }, -0.0103089));
  initVals.push_back(makePoint( { 125,   2,  1 }, -0.0145262));
  initVals.push_back(makePoint( { 126,   2,  1 },  0.0142726));
  initVals.push_back(makePoint( {   0,   3,  1 }, -0.0127882));
  initVals.push_back(makePoint( {   1,   3,  1 }, -0.0340823));
  initVals.push_back(makePoint( {   2,   3,  1 }, -0.0145262));
  initVals.push_back(makePoint( {   3,   3,  1 }, -0.1606770));
  initVals.push_back(makePoint( {   4,   3,  1 },  0.0191507));
  initVals.push_back(makePoint( { 124,   3,  1 },  0.0191507));
  initVals.push_back(makePoint( { 125,   3,  1 },  0.1606770));
  initVals.push_back(makePoint( { 126,   3,  1 }, -0.0145262));
  initVals.push_back(makePoint( { 127,   3,  1 },  0.0340823));
  initVals.push_back(makePoint( {   1,   4,  1 },  0.1662350));
  initVals.push_back(makePoint( {   2,   4,  1 },  0.0103089));
  initVals.push_back(makePoint( {   3,   4,  1 },  0.0191507));
  initVals.push_back(makePoint( { 125,   4,  1 },  0.0191507));
  initVals.push_back(makePoint( { 126,   4,  1 }, -0.0103089));
  initVals.push_back(makePoint( { 127,   4,  1 },  0.1662350));
  initVals.push_back(makePoint( {   1,   5,  1 },  0.0224927));
  initVals.push_back(makePoint( {   2,   5,  1 },  0.0102342));
  initVals.push_back(makePoint( { 126,   5,  1 },  0.0102342));
  initVals.push_back(makePoint( { 127,   5,  1 }, -0.0224927));
  initVals.push_back(makePoint( {   1, 123,  1 }, -0.0224927));
  initVals.push_back(makePoint( {   2, 123,  1 },  0.0102342));
  initVals.push_back(makePoint( { 126, 123,  1 },  0.0102342));
  initVals.push_back(makePoint( { 127, 123,  1 },  0.0224927));
  initVals.push_back(makePoint( {   1, 124,  1 },  0.1662350));
  initVals.push_back(makePoint( {   2, 124,  1 }, -0.0103089));
  initVals.push_back(makePoint( {   3, 124,  1 },  0.0191507));
  initVals.push_back(makePoint( { 125, 124,  1 },  0.0191507));
  initVals.push_back(makePoint( { 126, 124,  1 },  0.0103089));
  initVals.push_back(makePoint( { 127, 124,  1 },  0.1662350));
  initVals.push_back(makePoint( {   0, 125,  1 }, -0.0127882));
  initVals.push_back(makePoint( {   1, 125,  1 },  0.0340823));
  initVals.push_back(makePoint( {   2, 125,  1 }, -0.0145262));
  initVals.push_back(makePoint( {   3, 125,  1 },  0.1606770));
  initVals.push_back(makePoint( {   4, 125,  1 },  0.0191507));
  initVals.push_back(makePoint( { 124, 125,  1 },  0.0191507));
  initVals.push_back(makePoint( { 125, 125,  1 }, -0.1606770));
  initVals.push_back(makePoint( { 126, 125,  1 }, -0.0145262));
  initVals.push_back(makePoint( { 127, 125,  1 }, -0.0340823));
  initVals.push_back(makePoint( {   2, 126,  1 },  0.0142726));
  initVals.push_back(makePoint( {   3, 126,  1 }, -0.0145262));
  initVals.push_back(makePoint( {   4, 126,  1 }, -0.0103089));
  initVals.push_back(makePoint( {   5, 126,  1 },  0.0102342));
  initVals.push_back(makePoint( { 123, 126,  1 },  0.0102342));
  initVals.push_back(makePoint( { 124, 126,  1 },  0.0103089));
  initVals.push_back(makePoint( { 125, 126,  1 }, -0.0145262));
  initVals.push_back(makePoint( { 126, 126,  1 }, -0.0142726));
  initVals.push_back(makePoint( {   3, 127,  1 },  0.0340823));
  initVals.push_back(makePoint( {   4, 127,  1 },  0.1662350));
  initVals.push_back(makePoint( {   5, 127,  1 }, -0.0224927));
  initVals.push_back(makePoint( { 123, 127,  1 },  0.0224927));
  initVals.push_back(makePoint( { 124, 127,  1 },  0.1662350));
  initVals.push_back(makePoint( { 125, 127,  1 }, -0.0340823));
  initVals.push_back(makePoint( {   0,   0,  2 },  0.1558300));
  initVals.push_back(makePoint( {   2,   0,  2 },  0.1327290));
  initVals.push_back(makePoint( { 126,   0,  2 },  0.1327290));
  initVals.push_back(makePoint( {   1,   1,  2 }, -0.0358995));
  initVals.push_back(makePoint( {   2,   1,  2 },  0.1359910));
  initVals.push_back(makePoint( {   3,   1,  2 },  0.0529617));
  initVals.push_back(makePoint( { 125,   1,  2 },  0.0529617));
  initVals.push_back(makePoint( { 126,   1,  2 }, -0.1359910));
  initVals.push_back(makePoint( { 127,   1,  2 }, -0.0358995));
  initVals.push_back(makePoint( {   0,   2,  2 },  0.1327290));
  initVals.push_back(makePoint( {   1,   2,  2 },  0.1359910));
  initVals.push_back(makePoint( {   2,   2,  2 }, -0.0743355));
  initVals.push_back(makePoint( {   3,   2,  2 }, -0.0204280));
  initVals.push_back(makePoint( { 125,   2,  2 },  0.0204280));
  initVals.push_back(makePoint( { 126,   2,  2 }, -0.0743355));
  initVals.push_back(makePoint( { 127,   2,  2 }, -0.1359910));
  initVals.push_back(makePoint( {   1,   3,  2 },  0.0529617));
  initVals.push_back(makePoint( {   2,   3,  2 }, -0.0204280));
  initVals.push_back(makePoint( {   5,   3,  2 },  0.0120951));
  initVals.push_back(makePoint( { 123,   3,  2 },  0.0120951));
  initVals.push_back(makePoint( { 126,   3,  2 },  0.0204280));
  initVals.push_back(makePoint( { 127,   3,  2 },  0.0529617));
  initVals.push_back(makePoint( {   3,   5,  2 },  0.0120951));
  initVals.push_back(makePoint( { 125,   5,  2 },  0.0120951));
  initVals.push_back(makePoint( {   3, 123,  2 },  0.0120951));
  initVals.push_back(makePoint( { 125, 123,  2 },  0.0120951));
  initVals.push_back(makePoint( {   1, 125,  2 },  0.0529617));
  initVals.push_back(makePoint( {   2, 125,  2 },  0.0204280));
  initVals.push_back(makePoint( {   5, 125,  2 },  0.0120951));
  initVals.push_back(makePoint( { 123, 125,  2 },  0.0120951));
  initVals.push_back(makePoint( { 126, 125,  2 }, -0.0204280));
  initVals.push_back(makePoint( { 127, 125,  2 },  0.0529617));
  initVals.push_back(makePoint( {   0, 126,  2 },  0.1327290));
  initVals.push_back(makePoint( {   1, 126,  2 }, -0.1359910));
  initVals.push_back(makePoint( {   2, 126,  2 }, -0.0743355));
  initVals.push_back(makePoint( {   3, 126,  2 },  0.0204280));
  initVals.push_back(makePoint( { 125, 126,  2 }, -0.0204280));
  initVals.push_back(makePoint( { 126, 126,  2 }, -0.0743355));
  initVals.push_back(makePoint( { 127, 126,  2 },  0.1359910));
  initVals.push_back(makePoint( {   1, 127,  2 }, -0.0358995));
  initVals.push_back(makePoint( {   2, 127,  2 }, -0.1359910));
  initVals.push_back(makePoint( {   3, 127,  2 },  0.0529617));
  initVals.push_back(makePoint( { 125, 127,  2 },  0.0529617));
  initVals.push_back(makePoint( { 126, 127,  2 },  0.1359910));
  initVals.push_back(makePoint( { 127, 127,  2 }, -0.0358995));
  initVals.push_back(makePoint( {   0,   0, 62 },  0.1558300));
  initVals.push_back(makePoint( {   2,   0, 62 },  0.1327290));
  initVals.push_back(makePoint( { 126,   0, 62 },  0.1327290));
  initVals.push_back(makePoint( {   1,   1, 62 }, -0.0358995));
  initVals.push_back(makePoint( {   2,   1, 62 },  0.1359910));
  initVals.push_back(makePoint( {   3,   1, 62 },  0.0529617));
  initVals.push_back(makePoint( { 125,   1, 62 },  0.0529617));
  initVals.push_back(makePoint( { 126,   1, 62 }, -0.1359910));
  initVals.push_back(makePoint( { 127,   1, 62 }, -0.0358995));
  initVals.push_back(makePoint( {   0,   2, 62 },  0.1327290));
  initVals.push_back(makePoint( {   1,   2, 62 },  0.1359910));
  initVals.push_back(makePoint( {   2,   2, 62 }, -0.0743355));
  initVals.push_back(makePoint( {   3,   2, 62 }, -0.0204280));
  initVals.push_back(makePoint( { 125,   2, 62 },  0.0204280));
  initVals.push_back(makePoint( { 126,   2, 62 }, -0.0743355));
  initVals.push_back(makePoint( { 127,   2, 62 }, -0.1359910));
  initVals.push_back(makePoint( {   1,   3, 62 },  0.0529617));
  initVals.push_back(makePoint( {   2,   3, 62 }, -0.0204280));
  initVals.push_back(makePoint( {   5,   3, 62 },  0.0120951));
  initVals.push_back(makePoint( { 123,   3, 62 },  0.0120951));
  initVals.push_back(makePoint( { 126,   3, 62 },  0.0204280));
  initVals.push_back(makePoint( { 127,   3, 62 },  0.0529617));
  initVals.push_back(makePoint( {   3,   5, 62 },  0.0120951));
  initVals.push_back(makePoint( { 125,   5, 62 },  0.0120951));
  initVals.push_back(makePoint( {   3, 123, 62 },  0.0120951));
  initVals.push_back(makePoint( { 125, 123, 62 },  0.0120951));
  initVals.push_back(makePoint( {   1, 125, 62 },  0.0529617));
  initVals.push_back(makePoint( {   2, 125, 62 },  0.0204280));
  initVals.push_back(makePoint( {   5, 125, 62 },  0.0120951));
  initVals.push_back(makePoint( { 123, 125, 62 },  0.0120951));
  initVals.push_back(makePoint( { 126, 125, 62 }, -0.0204280));
  initVals.push_back(makePoint( { 127, 125, 62 },  0.0529617));
  initVals.push_back(makePoint( {   0, 126, 62 },  0.1327290));
  initVals.push_back(makePoint( {   1, 126, 62 }, -0.1359910));
  initVals.push_back(makePoint( {   2, 126, 62 }, -0.0743355));
  initVals.push_back(makePoint( {   3, 126, 62 },  0.0204280));
  initVals.push_back(makePoint( { 125, 126, 62 }, -0.0204280));
  initVals.push_back(makePoint( { 126, 126, 62 }, -0.0743355));
  initVals.push_back(makePoint( { 127, 126, 62 },  0.1359910));
  initVals.push_back(makePoint( {   1, 127, 62 }, -0.0358995));
  initVals.push_back(makePoint( {   2, 127, 62 }, -0.1359910));
  initVals.push_back(makePoint( {   3, 127, 62 },  0.0529617));
  initVals.push_back(makePoint( { 125, 127, 62 },  0.0529617));
  initVals.push_back(makePoint( { 126, 127, 62 },  0.1359910));
  initVals.push_back(makePoint( { 127, 127, 62 }, -0.0358995));
  initVals.push_back(makePoint( {   3,   0, 63 }, -0.0127882));
  initVals.push_back(makePoint( { 125,   0, 63 }, -0.0127882));
  initVals.push_back(makePoint( {   3,   1, 63 }, -0.0340823));
  initVals.push_back(makePoint( {   4,   1, 63 },  0.1662350));
  initVals.push_back(makePoint( {   5,   1, 63 },  0.0224927));
  initVals.push_back(makePoint( { 123,   1, 63 }, -0.0224927));
  initVals.push_back(makePoint( { 124,   1, 63 },  0.1662350));
  initVals.push_back(makePoint( { 125,   1, 63 },  0.0340823));
  initVals.push_back(makePoint( {   2,   2, 63 }, -0.0142726));
  initVals.push_back(makePoint( {   3,   2, 63 }, -0.0145262));
  initVals.push_back(makePoint( {   4,   2, 63 },  0.0103089));
  initVals.push_back(makePoint( {   5,   2, 63 },  0.0102342));
  initVals.push_back(makePoint( { 123,   2, 63 },  0.0102342));
  initVals.push_back(makePoint( { 124,   2, 63 }, -0.0103089));
  initVals.push_back(makePoint( { 125,   2, 63 }, -0.0145262));
  initVals.push_back(makePoint( { 126,   2, 63 },  0.0142726));
  initVals.push_back(makePoint( {   0,   3, 63 }, -0.0127882));
  initVals.push_back(makePoint( {   1,   3, 63 }, -0.0340823));
  initVals.push_back(makePoint( {   2,   3, 63 }, -0.0145262));
  initVals.push_back(makePoint( {   3,   3, 63 }, -0.1606770));
  initVals.push_back(makePoint( {   4,   3, 63 },  0.0191507));
  initVals.push_back(makePoint( { 124,   3, 63 },  0.0191507));
  initVals.push_back(makePoint( { 125,   3, 63 },  0.1606770));
  initVals.push_back(makePoint( { 126,   3, 63 }, -0.0145262));
  initVals.push_back(makePoint( { 127,   3, 63 },  0.0340823));
  initVals.push_back(makePoint( {   1,   4, 63 },  0.1662350));
  initVals.push_back(makePoint( {   2,   4, 63 },  0.0103089));
  initVals.push_back(makePoint( {   3,   4, 63 },  0.0191507));
  initVals.push_back(makePoint( { 125,   4, 63 },  0.0191507));
  initVals.push_back(makePoint( { 126,   4, 63 }, -0.0103089));
  initVals.push_back(makePoint( { 127,   4, 63 },  0.1662350));
  initVals.push_back(makePoint( {   1,   5, 63 },  0.0224927));
  initVals.push_back(makePoint( {   2,   5, 63 },  0.0102342));
  initVals.push_back(makePoint( { 126,   5, 63 },  0.0102342));
  initVals.push_back(makePoint( { 127,   5, 63 }, -0.0224927));
  initVals.push_back(makePoint( {   1, 123, 63 }, -0.0224927));
  initVals.push_back(makePoint( {   2, 123, 63 },  0.0102342));
  initVals.push_back(makePoint( { 126, 123, 63 },  0.0102342));
  initVals.push_back(makePoint( { 127, 123, 63 },  0.0224927));
  initVals.push_back(makePoint( {   1, 124, 63 },  0.1662350));
  initVals.push_back(makePoint( {   2, 124, 63 }, -0.0103089));
  initVals.push_back(makePoint( {   3, 124, 63 },  0.0191507));
  initVals.push_back(makePoint( { 125, 124, 63 },  0.0191507));
  initVals.push_back(makePoint( { 126, 124, 63 },  0.0103089));
  initVals.push_back(makePoint( { 127, 124, 63 },  0.1662350));
  initVals.push_back(makePoint( {   0, 125, 63 }, -0.0127882));
  initVals.push_back(makePoint( {   1, 125, 63 },  0.0340823));
  initVals.push_back(makePoint( {   2, 125, 63 }, -0.0145262));
  initVals.push_back(makePoint( {   3, 125, 63 },  0.1606770));
  initVals.push_back(makePoint( {   4, 125, 63 },  0.0191507));
  initVals.push_back(makePoint( { 124, 125, 63 },  0.0191507));
  initVals.push_back(makePoint( { 125, 125, 63 }, -0.1606770));
  initVals.push_back(makePoint( { 126, 125, 63 }, -0.0145262));
  initVals.push_back(makePoint( { 127, 125, 63 }, -0.0340823));
  initVals.push_back(makePoint( {   2, 126, 63 },  0.0142726));
  initVals.push_back(makePoint( {   3, 126, 63 }, -0.0145262));
  initVals.push_back(makePoint( {   4, 126, 63 }, -0.0103089));
  initVals.push_back(makePoint( {   5, 126, 63 },  0.0102342));
  initVals.push_back(makePoint( { 123, 126, 63 },  0.0102342));
  initVals.push_back(makePoint( { 124, 126, 63 },  0.0103089));
  initVals.push_back(makePoint( { 125, 126, 63 }, -0.0145262));
  initVals.push_back(makePoint( { 126, 126, 63 }, -0.0142726));
  initVals.push_back(makePoint( {   3, 127, 63 },  0.0340823));
  initVals.push_back(makePoint( {   4, 127, 63 },  0.1662350));
  initVals.push_back(makePoint( {   5, 127, 63 }, -0.0224927));
  initVals.push_back(makePoint( { 123, 127, 63 },  0.0224927));
  initVals.push_back(makePoint( { 124, 127, 63 },  0.1662350));
  initVals.push_back(makePoint( { 125, 127, 63 }, -0.0340823));

  // populate non-zero elements of array
  for (std::vector<point>::iterator it = initVals.begin(); it != initVals.end(); it++) {
    std::vector<int> k = std::get<0>(*it);	  
    double	     u = std::get<1>(*it);
                                                                                         
    int kz = k[2];
    int ky = k[1];
    int kx = k[0];
                                                                                         
    int index = kx + (gridSizes[2] * ky) + (gridSizes[2] * gridSizes[1] * kz);
                                                                                         
    data[index][0] = u * m_amplitude;
  }
} // end populateDataArray method
*/
