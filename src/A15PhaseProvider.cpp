#include "A15PhaseProvider.h"

#define _USE_MATH_DEFINES
#include <math.h>

FieldProvider A15PhaseProvider::generateInitialCondition(int gridSize) {

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

// method to reset field to initial A15 condition
void A15PhaseProvider::resetCondition(FieldProvider &field)
{
  // verify phaseID
  const int phaseID = field.getPhaseID();
  if (phaseID != m_phaseID)
    throw std::runtime_error("Could not reset FCC phase - incorrect phase ID");

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
void A15PhaseProvider::populateDataArray(double* dxVec, fftw_complex* realData, int numFieldElements, int* gridSizes)
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
  double coords[21][3] = {
    {0, 0, 0},
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1}, 
    {1, 1, 0},
    {1, 0, 1},
    {0, 1, 1},
    {1, 1, 1},
    {0.5, 0.5, 0.5},
    {0.75, 1, 0.5},
    {0.25, 1, 0.5},
    {0.5, 0.75, 0},
    {0.5, 0.25, 0},
    {1, 0.5, 0.25},
    {0, 0.5, 0.25},
    {1, 0.5, 0.75},
    {0, 0.5, 0.75},
    {0.5, 0.75, 1},
    {0.5, 0.25, 1},
    {0.75, 0, 0.5},
    {0.25, 0, 0.5}
  };

  double r0sqrd = m_period / 4.0;

  double sum = 0.0;

  for (int k = 0; k < Nz; k++) 
  for (int j = 0; j < Ny; j++) 
  for (int i = 0; i < Nx; i++) {

    int index = i + Nx * (j + Ny * k);

    double rx = i * dx;
    double ry = j * dx;
    double rz = k * dx;

    // loop over list of peaks:
    for (int m = 0; m < 21; m++) {
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

/*
void A15PhaseProvider::populateDataArray(double* dqVec, fftw_complex* cplxData, int numFieldElements, int* gridSizes)
{
  // reset grid spacing:
  double dq = 2 * M_PI / m_period;
  dqVec[0] = dq;
  dqVec[1] = dq;
  dqVec[2] = dq;

  // set all values in cplxData to zero
  for (int index = 0; index < numFieldElements; index++) {
    cplxData[index][0] = 0.0;
    cplxData[index][1] = 0.0;
  }

  // Fourier peaks describing A15 phase
  std::vector<point> initVals;
  initVals.push_back( makePoint( {  0,  0,  0 },    m_avDensity));          
  initVals.push_back( makePoint( { -2, -1,  0 },    0.0338277052131978)); //
  initVals.push_back( makePoint( {  2,  1,  0 },    0.0338277052131978)); //
  initVals.push_back( makePoint( {  2, -1,  0 },    0.0338277052131977)); //
  initVals.push_back( makePoint( { -2,  1,  0 },    0.0338277052131977)); //
  initVals.push_back( makePoint( { -1,  2,  0 },   -0.0338277052131977)); //
  initVals.push_back( makePoint( {  1, -2,  0 },   -0.0338277052131977)); //
  initVals.push_back( makePoint( {  1,  2,  0 },   -0.0338277052131977)); //
  initVals.push_back( makePoint( { -1, -2,  0 },   -0.0338277052131977)); //
  initVals.push_back( makePoint( {  0, -2,  1 },    0.0338228603092057)); //
  initVals.push_back( makePoint( { -2,  0,  1 },   -0.0338228603092056)); //
  initVals.push_back( makePoint( {  0,  2,  1 },    0.0338228603092055)); //
  initVals.push_back( makePoint( {  2,  0,  1 },   -0.0338228603092055)); //
  initVals.push_back( makePoint( { -1,  0,  2 },    0.0338139694707322)); //
  initVals.push_back( makePoint( {  0, -1,  2 },   -0.0338139694707320)); //
  initVals.push_back( makePoint( {  0,  1,  2 },   -0.0338139694707319)); //
  initVals.push_back( makePoint( {  1,  0,  2 },    0.0338139694707319)); //
  initVals.push_back( makePoint( {  0,  0,  2 },    0.0331853928814813)); //
  initVals.push_back( makePoint( {  0,  2,  0 },    0.0331847214985317)); //
  initVals.push_back( makePoint( {  2,  0,  0 },    0.0331847214985317)); //
  initVals.push_back( makePoint( {  0, -2,  0 },    0.0331847214985317)); //
  initVals.push_back( makePoint( { -2,  0,  0 },    0.0331847214985317)); //
  initVals.push_back( makePoint( {  1, -1,  2 },    0.0297069055330425)); //
  initVals.push_back( makePoint( { -1, -1,  2 },    0.0297069055330425)); //
  initVals.push_back( makePoint( {  1,  1,  2 },    0.0297069055330424)); //
  initVals.push_back( makePoint( { -1,  1,  2 },    0.0297069055330424)); //
  initVals.push_back( makePoint( { -1,  2,  1 },    0.0297065456817734)); //
  initVals.push_back( makePoint( { -2, -1,  1 },    0.0297065456817733)); //
  initVals.push_back( makePoint( { -1, -2,  1 },    0.0297065456817733)); //
  initVals.push_back( makePoint( {  2, -1,  1 },    0.0297065456817733)); //
  initVals.push_back( makePoint( { -2,  1,  1 },    0.0297065456817732)); //
  initVals.push_back( makePoint( {  2,  1,  1 },    0.0297065456817732)); //
  initVals.push_back( makePoint( {  1,  2,  1 },    0.0297065456817732)); //
  initVals.push_back( makePoint( {  1, -2,  1 },    0.0297065456817731)); //

  // load non-zero peaks into cplx data array
  for (std::vector<point>::iterator it = initVals.begin(); it != initVals.end(); it++) {
    std::vector<int> k = std::get<0>(*it);                                              
    double	     u = std::get<1>(*it);
                                                                                         
    int kz = k[2] < 0 ? k[2] + gridSizes[0] : k[2];
    int ky = k[1] < 0 ? k[1] + gridSizes[1] : k[1];
    int kx = k[0] < 0 ? k[0] + gridSizes[2] : k[0];
                                                                                         
    int index = kx + (gridSizes[2] * ky) + (gridSizes[2] * gridSizes[1] * kz);
                                                                                         
    cplxData[index][0] = u * m_amplitude;
  }
}
*/
