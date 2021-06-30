#ifndef _FIELD_PROVIDER_GUARD
#define _FIELD_PROVIDER_GUARD

#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>
#include <utility>

#include "fftw3.h"

// class represents the field stored on an d-dimensional rectangular grid
// field is a collection of (N1 * N2 * ... * Nd) complex numbers

class FieldProvider {

private:
    fftw_complex* m_realData; // pointer to array of complex data representing field in real space
    fftw_complex* m_cplxData; // pointer to array of complex data representing field in Fourier space
    fftw_complex* m_cplxTemp; // C2R transform overwrites cplx data - use this + memcpy so we don't lose cplx field
    fftw_complex* m_realTemp; // R2C transform apparently also requires this?

    fftw_plan m_fieldPlanC2R, m_fieldPlanR2C;   // plans that allow us to convert between real and complex fields

    int m_dimension;         // dimension of field - 1, 2, 3
    int *m_gridSizes;        // array of grid sizes for each dimension - length of vector is m_d
    int m_numFieldElements;  // total number of grid points (product of elements in m_N)
    // Note about grid sizes: they must be packed in this order (3D): (Nz, Ny, Nx) or (2D): (Ny, Nx)

    double *m_dx;   // array of real grid cell spacing - length is m_d
    double *m_dq;   // array of cplx grid cell spacing - length is m_d
    // Note about grid spacings: they must be packed in this order (3D): (dz, dy, dx) or (2D): (dy, dx)
    
    int m_phaseID;	// field knows which phase it's initialized to - used for reset

    double *m_dx0;    // array with initial grid spacing (for reset)
    double *m_data0;  // array with initial data (for reset)

public:
    /*
     * =====================================
     *              deleter
     * =====================================
     */

     ~FieldProvider() {
       free(m_gridSizes);
       
       free(m_dx);
       free(m_dq);

       fftw_free(m_realData);
       fftw_free(m_cplxData);
       fftw_free(m_cplxTemp);
       fftw_free(m_realTemp);
       
       if (m_fieldPlanC2R != NULL) fftw_destroy_plan(m_fieldPlanC2R);
       if (m_fieldPlanR2C != NULL) fftw_destroy_plan(m_fieldPlanR2C);

       free(m_dx0);
       free(m_data0);
     }

    /*
     * =====================================
     *             constructors
     * =====================================
     */

    // constructor to initialize FieldProvider directly
    FieldProvider(
      fftw_complex* data,
      const int dimension,
      int* gridSizes,
      double* dr,
      const bool real,
      int phaseId);

    // constructor to initialize FieldProvider from data stored in file
    FieldProvider(std::string fileName);

    /*
     * =====================================
     *      getter and setter functions
     * =====================================
     */

    // get/set pointer to data array:
    fftw_complex* getRealDataPointer() { return m_realData; }
    fftw_complex* getCplxDataPointer() { return m_cplxData; }
    fftw_complex* getRealTempPointer() { return m_realTemp; }
    fftw_complex* getCplxTempPointer() { return m_cplxTemp; }	
    void setData( fftw_complex* data, bool real )
    {
      if (real) {
	for (int i = 0; i < m_numFieldElements; i++) {
	  m_realData[i][0] = data[i][0];
	  m_realData[i][1] = data[i][1];
	}
        transformR2C();
      } else {
	for (int i = 0; i < m_numFieldElements; i++) {
	  m_cplxData[i][0] = data[i][0];
	  m_cplxData[i][1] = data[i][1];
	}
        transformC2R();
      }
    }

    // get number of dimensions, d
    int getDimension() { return m_dimension; }

    // get vector of length d containing grid sizes
    int* getGridSizes() {return m_gridSizes; }

    // get number of elements (equals length of data array)
    int getNumFieldElements() {return m_numFieldElements; }

    // get/set grid cell spacing (can change when we do box optimization)
    double* getDx() { return m_dx; }
    double* getDq() { return m_dq; }
    void setDx( double* dx ) {
      free(m_dx);
      m_dx = dx;
      updateDq();
    }
    void setDq( double* dq ) {
      free(m_dq);
      m_dq = dq;
      updateDx();
    }

    void updateDq() {
      for (int d = 0; d < m_dimension; d++) m_dq[d] = 2 * M_PI / (m_dx[d] * m_gridSizes[d]);
    }
                                                                                             
    void updateDx() {
      for (int d = 0; d < m_dimension; d++) m_dx[d] = 2 * M_PI / (m_dq[d] * m_gridSizes[d]);
    }

    // get phase ID
    int getPhaseID() { return m_phaseID; }

    /*
     * =====================================
     *      methods to write to file:
     * =====================================
     */
    void saveAsInitializer(std::string fileName, bool real);

    void saveForPlotting(std::string fileName, bool real);

    /*
     * =====================================
     *        transform methods:
     * =====================================
     */
    void transformR2C() {
      memcpy(m_realTemp, m_realData, sizeof(fftw_complex) * m_numFieldElements);
      fftw_execute(m_fieldPlanR2C); 
      for (int index = 0; index < m_numFieldElements; index++)
      {
        m_cplxData[index][0] /= m_numFieldElements;
	m_cplxData[index][1] /= m_numFieldElements;
      }
    }
    void transformC2R() {
      memcpy(m_cplxTemp, m_cplxData, sizeof(fftw_complex) * m_numFieldElements);
      fftw_execute(m_fieldPlanC2R);
      for (int i = 0; i < m_numFieldElements; i++)
        m_realData[i][1] = 0.0;
    }

    /*
     * =====================================
     *        compute laplacian
     * =====================================
     */
    void laplacian(double* laplacian);

    /*
     * =====================================
     *  reset field provider to initial data
     * =====================================
     */
    void reset();

}; // end FieldProvider class definition
#endif


