#include "FieldProvider.h"

#include <iomanip>
#include <stdexcept>
#include <stdlib.h>
#include <string>

// constructor to initialize fieldProvider from provided data
FieldProvider::FieldProvider(
  fftw_complex* data,
  const int dimension,
  int* gridSizes,
  double* dr,
  const bool real,
  int phaseID) 
{
  // dimension of field
  m_dimension = dimension;

  // grid sizes
  m_gridSizes = gridSizes;

  if (real) {
    m_dx = dr;
    m_dq = (double *) malloc(m_dimension * sizeof(double));
    updateDq();
  } else {
    m_dq = dr;
    m_dx = (double *) malloc(m_dimension * sizeof(double));
    updateDx();
  }

  // compute number of field elements:
  m_numFieldElements = 1;
  for (int d = 0; d < m_dimension; d++)
    m_numFieldElements *= gridSizes[d];

  // allocate memory for real and complex fields:
  m_realData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_numFieldElements);
  m_cplxData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_numFieldElements);
  m_cplxTemp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_numFieldElements);
  m_realTemp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_numFieldElements);

  // initialize plans
  // C2R plans use "BACKWARD" flag, R2C uses "FORWARD"
  m_fieldPlanR2C = fftw_plan_dft(m_dimension, m_gridSizes, m_realTemp, m_cplxData, FFTW_FORWARD,  FFTW_PATIENT);
  m_fieldPlanC2R = fftw_plan_dft(m_dimension, m_gridSizes, m_cplxTemp, m_realData, FFTW_BACKWARD, FFTW_PATIENT);

  // load data provided and transform:
  if (real) {
    for (int index = 0; index < m_numFieldElements; index++) {
      m_realData[index][0] = data[index][0];
      m_realData[index][1] = data[index][1];
    }

    transformR2C();
  } else {
    for (int index = 0; index < m_numFieldElements; index++) {
      m_cplxData[index][0] = data[index][0];
      m_cplxData[index][1] = data[index][1];
    }

    // perform C2R transform
    transformC2R();
  }

  // set phase ID
  m_phaseID = phaseID;
}

// constructor to initialize fieldProvider from file
FieldProvider::FieldProvider(std::string fileName) {

  std::ifstream inFile;

  // open file and verify
  inFile.open(fileName);
  if (!inFile) {
    std::string exceptionMessage = "unable to open file " + fileName;
    throw std::invalid_argument(exceptionMessage);
  }

  // read dimension of data (first line of file)
  inFile >> m_dimension;

  // read grid sizes - printed on next d lines
  m_gridSizes = (int*) malloc(sizeof(int) * m_dimension);
  m_numFieldElements = 1;
  for (int d = 0; d < m_dimension; d++) {
    inFile >> m_gridSizes[d];
    m_numFieldElements *= m_gridSizes[d];
  }

  // read grid cell spacing - next d lines
  double* dr = (double*) malloc(sizeof(double) * m_dimension);
  for (int d = 0; d < m_dimension; d++)
    inFile >> dr[d];

  // read phase ID
  inFile >> m_phaseID;

  // real or complex data ? - use to set grid spacing
  m_dx = (double*) malloc(sizeof(double) * m_dimension);
  m_dq = (double*) malloc(sizeof(double) * m_dimension);

  bool real;
  inFile >> real;
  if (real) setDx(dr);
  else setDq(dr);

  // allocate memory for real/cplx data and for temp data
  m_realData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_numFieldElements);
  m_cplxData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_numFieldElements);
  m_realTemp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_numFieldElements);
  m_cplxTemp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_numFieldElements);

  // initialize fftw plans
  m_fieldPlanR2C = fftw_plan_dft(m_dimension, m_gridSizes, m_realTemp, m_cplxData, FFTW_FORWARD,  FFTW_PATIENT);
  m_fieldPlanC2R = fftw_plan_dft(m_dimension, m_gridSizes, m_cplxTemp, m_realData, FFTW_BACKWARD, FFTW_PATIENT);

  // load data from file:
  double y;
  int index = 0;
  int pairIndex = 0;
  while (inFile >> y) {
    // make sure we don't have too many data points:
    if (index == m_numFieldElements) {
      std::string exceptionMessage = "too many data points in initializer file " + fileName;
      throw std::length_error(exceptionMessage);
    }

    if (pairIndex == 0) {
      if (real) m_realData[index][0] = y;
      else m_cplxData[index][0] = y;
      pairIndex++;
    } else {
      if (real) m_realData[index++][1] = y;
      else m_cplxData[index++][1] = y;
      pairIndex = 0;
    }

  }
  // make sure we have enough data points
  if (index != m_numFieldElements) {
    std::string exceptionMessage = "not enough data points in initializer file " + fileName;
    throw std::runtime_error(exceptionMessage);
  }

  // execute transform
  if (real) transformR2C();
  else transformC2R();

} // end of constructor

/*
 * =====================================
 *      methods to write to file:
 * =====================================
 */
void FieldProvider::saveAsInitializer(std::string fileName, bool real) {
  std::ofstream stream(fileName);
  stream << std::scientific << std::setprecision(10);

  // print dimension:
  stream << m_dimension << std::endl;

  // print grid sizes
  for (int d = 0; d < m_dimension; d++)
    stream << m_gridSizes[d] << std::endl;

  // print grid cell spacings
  for (int d = 0; d < m_dimension; d++) {
    if (real)
      stream << m_dx[d] << std::endl;
    else
      stream << m_dq[d] << std::endl;
  }

  // print phaseID
  stream << m_phaseID << std::endl;

  // print bool indicating whether data is real or complex
  stream << real << std::endl;

  // print data
  for (int index = 0; index < m_numFieldElements; index++) {
    if (real)
      stream << m_realData[index][0] << " " << m_realData[index][1] << std::endl;
    else
      stream << m_cplxData[index][0] << " " << m_cplxData[index][1] << std::endl;
  }

  stream.close();
} // end 'saveAsInitializer' method

void FieldProvider::saveForPlotting(std::string fileName, bool real) {
  std::ofstream stream(fileName);
  stream << std::scientific << std::setprecision(10);

  // print data
  for (int index = 0; index < m_numFieldElements; index++) {
    int i = index % m_gridSizes[m_dimension - 1];
    double x;
    if (real)
      x = i * m_dx[m_dimension - 1];
    else
      x = i * m_dq[m_dimension - 1];
    stream << x << ", ";

    if (m_dimension > 1) {
      int j = ((index - i) / m_gridSizes[m_dimension - 1]) % m_gridSizes[m_dimension - 2];
      double y;
      if (real)
        y = j * m_dx[m_dimension - 2];
      else
        y = j * m_dq[m_dimension - 2];
      stream << y << ", ";

      if (m_dimension == 3) {
        int k = (index - i - j * m_gridSizes[m_dimension - 1]) / (m_gridSizes[m_dimension - 1] * m_gridSizes[m_dimension - 2]);
        double z;
        if (real)
          z = k * m_dx[m_dimension - 3];
        else
          z = k * m_dq[m_dimension - 3];
        stream << z << ", ";
      }
    }

    double u;
    if (real)
      u = m_realData[index][0];
    else
      u = cplxMagnitude(m_cplxData[index]);

    stream << u << std::endl;
  }

  stream.close();
} // end 'saveForPlotting' method

// compute laplacian
void FieldProvider::laplacian(double* laplacian)
{
  int d = m_dimension;

  // get grid sizes and spacings:
  int    Nx(0),    Ny(0),    Nz(0);
  double dx2(0.0), dy2(0.0), dz2(0.0);
  Nx  = m_gridSizes[d-1];
  dx2 = m_dq[d-1] * m_dq[d-1];
  if (d > 1) {
    Ny  = m_gridSizes[d-2];
    dy2 = m_dq[d-2] * m_dq[d-2];
    if (d > 2) {
      Nz  = m_gridSizes[d-3];
      dz2 = m_dq[d-3] * m_dq[d-3];
    }
  }

  //
  int iMod(0), jMod(0), kMod(0);
  for (int index = 0; index < m_numFieldElements; index++) {
    int i = index % Nx;
    iMod  = i < Nx / 2 ? i : i - Nx;

    if (d > 1) {
      int j = (index - i) / Nx % Ny;
      jMod  = j < Ny / 2 ? j : j - Ny;

      if (d > 2) {
        int k = ((index - i) / Nx - j) / Ny;
        kMod  = k < Nz / 2 ? k : k - Nz;
      }
    } // end set iMod ... kMod

    double q2 = (iMod * iMod) * dx2 + (jMod * jMod) * dy2 + (kMod * kMod) * dz2;
    laplacian[index] = -q2;
  } // end loop over index
}


