// Copyright (c) 2014-2015, CosmicPy Developers
// Licensed under CeCILL 2.1 - see LICENSE.rst
#ifndef BESSELWINDOW_H
#define BESSELWINDOW_H
#include <math.h>
#include <iostream>
#include <fstream>
#include <boost/python.hpp>
#include <boost/numpy.hpp>
#include <boost/math/tr1.hpp>
#include <stdexcept>
#include <map>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "NumPyArrayData.h"

#define ASSERT_THROW(a,msg) if (!(a)) throw std::runtime_error(msg);

namespace bp = boost::python;
namespace np = boost::numpy;

namespace bm = boost::math::tr1;

class BesselWindow
{
  
public:
  
    BesselWindow(int64_t Nmax, int64_t Lmax, double Kmax, const char *filename="qlnTable.dat");
    
    ~BesselWindow();
    
    /*! Computes all the values of the radial grid of order \a l .*/
    np::ndarray rgrid(int64_t l);
    
    /*! Computes all the values of the wavenumber grid of order \a l for a given \a rmax. */
    np::ndarray kgrid(int64_t l, double rmax);
    
    /*! Computes Bessel window using the k_ln sampling (based on \a rmax) for k1 and an optimal sampling for k2.
     * \a l specifies the order of the transform, \a s is the comoving distance computed from
     * observed redshifts for a given cosmology, \a phi is the window function*/
    bp::tuple tabulated(int64_t l, double rmax, double kmax, np::ndarray& s, np::ndarray& phi);
    
    /*! Computes Bessel window using a custom sampling \a k for k1 and an optimal sampling for k2.
     * \a l specifies the order of the transform, \a s is the comoving distance computed from
     * observed redshifts for a given cosmology, \a phi is the window function*/
    bp::tuple optimal(int64_t l, double rmax, np::ndarray& k, np::ndarray& s, np::ndarray& phi);
   
    np::ndarray custom(int64_t l, np::ndarray &k1, np::ndarray &k2, np::ndarray &s, np::ndarray &phi);
    
    double getPrecision(){
	return precision;
    }
    void setPrecision(double p){
	precision = p;
	clearBuffer();
    }
    
    int64_t getNpoints(){
	return npoints;
    }
    void setNpoints(int64_t n){
	npoints = n;
	clearBuffer();
    }

private:
    double precision;
    int64_t npoints;
  
    /*! Structure to store precomputed values to preform the fast discrete Bessel Transform */
    struct preCompArray {
        long Nt;
        long nk1;
        long nk2;
        double rmax;
        double kmax;
        long N1;
        double * factor;
        double * k2;
        double * k1;
        double * Matk1;
        double * Matk2;
	long * k_interp_sup;
	double * k_interp_ratio;
    };
    

    /*! Simple comparator to order the map of precomputed arrays */
    struct classcomp {
        bool operator() (const int64_t& lhs, const int64_t& rhs) const
        {
            return lhs<rhs;
        }
    };
    
    /*! Map to store precomputed arrays to preform the fast discrete Bessel Transform */
    std::map<int64_t, preCompArray*, classcomp> buffer;
    
    /*! Clears the buffer and releases the allocated memory */
    void clearBuffer();
    
    /*! Computes values to preform the fast discrete Bessel Transform */
    preCompArray *preComputeArray(int64_t l, double rmax, double kmax, NumPyArrayData<double> &phi);
    
    /*! Frees precomputed array*/
    void freePreCompArray(preCompArray *array) ;
    
    /*! Loads precomputed zeros of Bessel functions from file \a filename.
     Or if the file does not exist, computes the zeros and stores them in \a filename*/
    void loadBesselZeroes(const char* filename);
    
    /*! Computes the width of the window function for order \a l and at a given scale \a k*/
    double getWidth(int64_t l, double k, double *phi);
    
            
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    
     /*! Structure to store parameters of evalWindow */
    struct eval_params {
        int l;
        double k0;
        double *phi;
        double *matrix;
        double a;
	int N;
	double K;
	double *qln;
    };
    
    /*! Utility function used by getWidth */
    static double evalWindow(double k, void* params);
    
    double *qln;   /*!< Array to store the zeros of the Bessel functions. */

    int64_t N, L;
    double K;
      
};

#endif // BESSELWINDOW_H
