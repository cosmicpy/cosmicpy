// Copyright (c) 2014-2015, CosmicPy Developers
// Licensed under CeCILL 2.1 - see LICENSE.rst
#include "besselwindow.h"
#include <boost/math/special_functions/bessel.hpp>

BesselWindow::BesselWindow(int64_t Nmax, int64_t Lmax, double Kmax, const char* filename)
{
    N=Nmax;
    L=Lmax;
    K=Kmax;

    npoints = 100;
    precision = 1e-3;
    
    // Allocate the qln table and copy over the zeros
    qln = (double*) malloc((L+1)*N*sizeof(double));
    loadBesselZeroes(filename);
    
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
}

BesselWindow::~BesselWindow()
{
    free(qln);
    gsl_root_fsolver_free(s);
}

void BesselWindow::freePreCompArray(preCompArray *array) {
    free(array->Matk1);
    free(array->Matk2);
    free(array->factor);
    free(array->k2);
    free(array->k1);
    free(array->k_interp_sup);
    free(array->k_interp_ratio);
    free(array);
}

void BesselWindow::clearBuffer()
{
  
  for (std::map<int64_t, preCompArray*, classcomp>::iterator it=buffer.begin(); it!=buffer.end(); ++it){
      freePreCompArray(it->second);
  }
  buffer.clear();
}


np::ndarray BesselWindow::rgrid(int64_t l)
{
    // Test the value of the order l, must be below the limit Lmax
    ASSERT_THROW( (l < L) , "Requested Bessel order l above Lmax");

    // Create a new double array to store the radial grid
    np::ndarray grid = np::zeros(bp::make_tuple(N), np::dtype::get_builtin<double>());
    // Object to easily access the data in the array
    NumPyArrayData<double> grid_data(grid);

    for (int64_t n=0 ; n < N; ++n) {
        grid_data(n) = qln[N*l + n]/K;
    }

    return grid;
}

boost::numpy::ndarray BesselWindow::kgrid(int64_t l, double rmax)
{
    // Test the value of the order l, must be below the limit Lmax
    ASSERT_THROW( (l < L) , "Requested Bessel order l above Lmax");

    // Create a new double array to store the radial grid
    np::ndarray grid = np::zeros(bp::make_tuple(N), np::dtype::get_builtin<double>());
    // Object to easily access the data in the array
    NumPyArrayData<double> grid_data(grid);

    for (int64_t n=0 ; n < N; ++n) {
        grid_data(n) = qln[N*l + n]/rmax;
    }

    return grid;
}


boost::python::tuple BesselWindow::tabulated(int64_t l, double rmax, double kmax, np::ndarray& s, np::ndarray& phi)
{
  
    // Objects to easily access the data in the arrays
    NumPyArrayData<double> s_data(s);
    NumPyArrayData<double> phi_data(phi);

    // Normalisation factor for the discrete transform
    double normFactor= M_PI/pow(K,3);
    
    // Request a precomputed array
    preCompArray *arr = preComputeArray(l, rmax, kmax, phi_data);

    // Create two new double arrays to store the wave number grids
    np::ndarray k1 = np::zeros(bp::make_tuple(arr->nk1), np::dtype::get_builtin<double>());
    np::ndarray k2 = np::zeros(bp::make_tuple(arr->nk2), np::dtype::get_builtin<double>());

    // Create double array to store the output window function
    np::ndarray w = np::zeros(bp::make_tuple(arr->nk1, arr->nk2), np::dtype::get_builtin<double>());

    // Objects to easily access the data in the arrays
    NumPyArrayData<double> k1_data(k1);
    NumPyArrayData<double> k2_data(k2);
    NumPyArrayData<double> w_data(w);

    long nk2 = arr->nk2;
    long nk1 = arr->nk1;

    // First, interpolate the window function phi to match the k2 grid
    double * phi_interp = (double *) malloc(nk2*arr->Nt*sizeof(double));
    #pragma omp parallel for schedule(guided)
    for(long indk2=0; indk2 < nk2; indk2++) {
        for(long p=0; p < arr->Nt; p++) {
            phi_interp[indk2*arr->Nt + p] = phi_data(p, arr->k_interp_sup[indk2]-1) + arr->k_interp_ratio[indk2]*(phi_data(p, arr->k_interp_sup[indk2]) - phi_data(p, arr->k_interp_sup[indk2] -1 ));
        }
    }

    #pragma omp parallel for schedule(guided)
    for(long indk1=0; indk1 < nk1; indk1++) {
        long offset1 = arr->N1*indk1;
        k1_data(indk1) = arr->k1[indk1];

        // Recomputing Matk1 to take s into account
        for(int64_t p=0; p < arr->N1 ; ++p) {
            if(qln[N*l +p]/K >= arr->rmax)
                break;
            arr->Matk1[offset1 + p] = arr->factor[p]*bm::sph_bessel<double>(l, s_data(p) * arr->k1[indk1]);
        }

        for(long indk2=0; indk2 < nk2; indk2++) {
            long offset2 = arr->Nt*indk2;
            double temp = 0;
            for(int64_t p=0; p < arr->Nt ; ++p) {

                if(qln[N*l +p]/K >= arr->rmax)
                    break;

                temp += arr->Matk1[offset1 + p] * arr->Matk2[offset2 + p] * phi_interp[offset2 + p];
            }
            w_data(indk1, indk2) = temp*arr->k1[indk1];
        }
    }

    #pragma omp parallel for schedule(guided)
    for(long indk2=0; indk2 < nk2; indk2++) {
        k2_data(indk2) = arr->k2[indk2];
    }

    free(phi_interp);

    return bp::make_tuple(k1, k2, w);
}

boost::python::tuple BesselWindow::optimal(int64_t l, double rmax, boost::numpy::ndarray& k, boost::numpy::ndarray& s, boost::numpy::ndarray& phi)
{
    double normFactor = M_PI/pow(K,3);
  
    // Objects to easily access the data in the input arrays
    NumPyArrayData<double> k_data(k);
    NumPyArrayData<double> s_data(s);
    NumPyArrayData<double> phi_data(phi);
    
    long nk1 = k.shape(0);
    long nk2 = npoints;
    
    // Create two new double arrays to store the wave number grids
    np::ndarray k2 = np::zeros(bp::make_tuple(nk1, npoints), np::dtype::get_builtin<double>());

    // Create double array to store the output window function
    np::ndarray w = np::zeros(bp::make_tuple(nk1, npoints), np::dtype::get_builtin<double>());

    // Objects to easily access the data in the arrays
    NumPyArrayData<double> k2_data(k2);
    NumPyArrayData<double> w_data(w);
    
    double * Matk2 = (double *) malloc(npoints*N*sizeof(double));
        
    long nk2_glob = npoints/2;
    long nk2_loc = npoints - nk2_glob;

    double *Matk1 = (double *) malloc(nk1*N*sizeof(double));
    double *glob_Mat = (double *) malloc(nk2_glob*N*sizeof(double));
    double *loc_Mat = (double *) malloc(nk2_loc*N*sizeof(double));
    double *glob_k   = (double *) malloc(nk2_glob*sizeof(double));
    double *loc_k   = (double *) malloc(nk2_loc*sizeof(double));
    double *vec  = (double *) malloc(nk2_glob*sizeof(double));
    double *kappaln = (double *) malloc(N*sizeof(double));

    for (int64_t q=0; q < nk2_glob ; ++q) {
        glob_k[q] = k_data(nk1-1)/((double) (nk2_glob+1.0))*(q+1);
    }
    
    // Computing base transform matrix for each k1
    #pragma omp parallel for schedule(guided)
    for (int64_t p=0; p < N ; ++p) {
        double kpln = 1.0/bm::sph_bessel<double>(l+1, qln[N*l + p]);
        kappaln[p] = kpln;
        kpln *= kpln*normFactor;
        for (int64_t q=0; q < nk1 ; ++q) {
            Matk1[q*N + p] = kpln* bm::sph_bessel<double>(l, s_data(p)*k_data(q));
        }
    }

    // Tabulating global matrix
    long kinterp=0;
    for (int64_t q=0; q < nk2_glob ; ++q) {
	while( qln[N*l + kinterp]/rmax  < glob_k[q] && kinterp < N-1){
	  kinterp++;
	}
        double tempVal =0;
	#pragma omp parallel for schedule(guided), reduction(+:tempVal)
        for (int64_t p=0; p < N ; ++p) {
            glob_Mat[q*N + p] = bm::sph_bessel<double>(l, qln[N*l + p]*glob_k[q]/K);
            tempVal += phi_data(p,kinterp)* normFactor*pow(glob_Mat[q*N + p]* kappaln[p],2);
        }
        vec[q] = tempVal;
    }

    // Computing maximum k
    double km =0;
    double val=0;
    long jm = 0;
    for(int j=0; j <nk2_glob; j++) {
        if(glob_k[j]* vec[j] >= val) {
            km = glob_k[j];
	    jm = j;
            val = glob_k[j]* vec[j];
        }
    }
    
    double * phi_arr = (double *) malloc(sizeof(double)*N);
    for (long i =0; i < N; i++) {
        phi_arr[i] = phi_data(i, jm);
    }
    double width = 1.5*getWidth(l, km, phi_arr);
    free(phi_arr);

    for(int64_t indk1=0; indk1 < nk1; indk1++) {

        double k = k_data(indk1);
        // Computing local matrix
        // First, compute the range of k values to probe
        double kmin = std::max(k - width/2.0, 0.0);
        double kmax = std::min(k + width/2.0, K);

        // Test if we overlap with the previously computed values
        long ind=0;

        // Computing missing elements for the local matrix
        for(int64_t indk2=0; indk2 < nk2_loc; indk2++) {
            if(indk2 < nk2_loc/2) {
                loc_k[indk2] = kmin + 2.0*(k - kmin)/((double) (nk2_loc+1))*(indk2+1);
            } else if(indk2 > nk2_loc/2) {
                loc_k[indk2] = k + 2.0*(kmax - k)/((double) nk2_loc)*(indk2 - nk2_loc/2) ;
            } else {
                loc_k[indk2] = k;
            }

            #pragma omp parallel for
            for (int64_t p=0; p < N ; ++p) {
                loc_Mat[indk2*N + p] = bm::sph_bessel<double>(l, qln[N*l + p]*loc_k[indk2]/K);
            }
        }

        // Ok, we have everything we need
        long ind_loc =0;
        long ind_glob =0;
        for(int64_t indk2=0; indk2 < nk2; indk2++) {

            if(ind_glob >= nk2_glob || (ind_loc < nk2_loc && loc_k[ind_loc] <= glob_k[ind_glob])) {
                k2_data(indk1, indk2) = loc_k[ind_loc];

                for(int64_t p=0; p < N ; ++p) {
                    Matk2[N*indk2 + p] = loc_Mat[ind_loc*N + p];
                }

                ind_loc++;
            } else {
                k2_data(indk1, indk2) = glob_k[ind_glob];
                for(int64_t p=0; p < N ; ++p) {
                    Matk2[N*indk2 + p] = glob_Mat[ind_glob*N + p] ;
                }
                ind_glob++;
            }
        }
        //TODO: Simple nearest value done here, a proper interpolation would be better.
        kinterp=0;
        for(int64_t indk2=0; indk2 < nk2; indk2++) {
            long offsetk2 = indk2*N;
	    while( qln[N*l + kinterp]/rmax  < glob_k[indk2] && kinterp < N - 1){
	      kinterp++;
	    }
	    
            double temp =0;
            for(int64_t p=0; p < N ; ++p) {
                temp += Matk2[indk2*N + p] * phi_data(p, kinterp)* Matk1[indk1*N + p];
            }
            w_data(indk1, indk2) = k_data(indk1) * temp;
        }
    }

    free(Matk2);
    free(Matk1);
    free(glob_Mat);
    free(loc_Mat);
    free(glob_k);
    free(loc_k);
    free(kappaln);
    free(vec);
    
    return bp::make_tuple(k2, w);
    
}

boost::numpy::ndarray BesselWindow::custom(int64_t l, boost::numpy::ndarray& k1, boost::numpy::ndarray& k2, boost::numpy::ndarray& s, boost::numpy::ndarray& phi)
{
    // Objects to easily access the data in the arrays
    NumPyArrayData<double> s_data(s);
    NumPyArrayData<double> phi_data(phi);
    NumPyArrayData<double> k1_data(k1);
    NumPyArrayData<double> k2_data(k2);
    
    double normFactor = M_PI/pow(K,3);
    
    long nk1 = k1.shape(0);
    long nk2 = k2.shape(0);
    
    // Create double array to store the output window function
    np::ndarray w = np::zeros(bp::make_tuple(nk1, nk2), np::dtype::get_builtin<double>());

    // Objects to easily access the data in the arrays
    NumPyArrayData<double> w_data(w);

    double *matrix1 = (double *) malloc(N*nk1*sizeof(double));
    double *matrix2 = (double *) malloc(N*nk2*sizeof(double));
    
    #pragma omp parallel
    {
        #pragma omp for schedule(guided)
        for (int64_t p=0; p < N ; ++p) {
            double tempValue =  normFactor/pow(bm::sph_bessel<double>(l+1, qln[N*l + p]),2);

            for (int64_t q=0; q< nk1 ; ++q) {
                matrix1[q*N + p] = tempValue* bm::sph_bessel<double>(l, s_data(p)*k1_data(q));
            }

            for (int64_t q=0; q< nk2 ; ++q) {
                matrix2[q*N + p] = phi_data(p, q) * bm::sph_bessel<double>(l, qln[N*l + p]*k2_data(q)/K);
            }
        }

        #pragma omp for schedule(guided)
        for (int64_t indk1=0; indk1<nk1; indk1++) {
            for (int64_t indk2=0; indk2<nk2; indk2++) {
                double tempValue = 0;
                for (int64_t p=0; p<N; p++) {
                    tempValue += matrix2[indk2*N + p] * matrix1[indk1*N + p];
                }
                w_data(indk1,indk2) = k1_data(indk1)*tempValue;
            }
        }
    }
    free(matrix1);
    free(matrix2);
   
   return w;  
}


BesselWindow::preCompArray* BesselWindow::preComputeArray(int64_t l, double rmax, double kmax, NumPyArrayData< double >& phi)
{
    double normFactor = M_PI/pow(K,3);
    preCompArray *arr = NULL;
    bool recompute = false;

    // Attempts to retrieve a precomputed array for multipole l
    try {
        arr = buffer.at(l);
        // Check if it matches the input parameters, if not it is recomputed
        if(arr->rmax != rmax || arr->kmax != kmax) {
            freePreCompArray(arr);
            buffer.erase(l);
            recompute = true;
        }
    } catch (std::out_of_range& oor) {
        // The array does not exists, in which case it is computed
        recompute = true;
    }
    
    // If the array was successfully recovered, it is returned...
    if(!recompute && arr != NULL) {
        return arr;
    }

    // Otherwise, it is recomputed
    ///////////////////////////////////////////////////////////////
    arr = new preCompArray;
    arr->kmax = kmax;
    arr->rmax = rmax;
    arr->Nt = 0;
    arr->nk1 = 0;

    // Computing the number of required points nk1
    for(long p=2; p < N; p++) {
        if(  arr->nk1 == 0 && qln[N*l +p]/rmax >= kmax) {
            arr->nk1 = p;
        }

        if(arr->Nt == 0 && qln[N*l +p]/K >= rmax) {
            arr->Nt = p;
        }
    }
    if(arr->Nt == 0)
        arr->Nt = N;
    arr->N1 = std::max(arr->Nt,arr->nk1);

    // Allocating the first Transform array, Just the triangular matrix
    arr->k1 = (double *) malloc(sizeof(double)*arr->nk1);
    arr->Matk1 = (double *) malloc(sizeof(double)*arr->N1*arr->N1);
    arr->factor = (double *) malloc(sizeof(double)*arr->N1);

    // Now, compute the first matrix
    #pragma omp parallel for schedule(guided)
    for (int64_t p=0; p < arr->N1 ; ++p) {
        if(p < arr->nk1)
            arr->k1[p]= qln[N*l + p]/rmax;

        double r = qln[N*l + p]/K;
        double bess =0;
        for (int64_t q=0; q < p ; ++q) {
            bess = bm::sph_bessel<double>(l, r*qln[N*l + q]/rmax);
            arr->Matk1[q*arr->N1 + p] = bess;
            arr->Matk1[p*arr->N1 + q] = bess;
        }
        bess = bm::sph_bessel<double>(l, r*qln[N*l + p]/rmax);
        arr->Matk1[p*arr->N1 + p] = bess;
        arr->factor[p] = normFactor/pow(bm::sph_bessel<double>(l+1, qln[N*l + p]),2);
    }

    // Locating maximum of the bessel window to compute the width
    double km =0;
    double val=0;
    long jm =0;
    for(int j=0; j < arr->N1; j++) {
        double vec =0;
        for (int64_t q=0; q < arr->N1 ; ++q) {
            vec += phi(q, j)*arr->Matk1[j*arr->N1 + q]*arr->Matk1[j*arr->N1 + q]*arr->factor[q];
            arr->Matk1[j*arr->N1 + q] *= arr->factor[q];
        }
        if(j < arr->nk1) {
            vec *= arr->k1[j];
            if(vec >= val) {
                km = qln[N*l + j]/rmax;
                jm = j;
                val = vec;
            }
        }
    }

    // Now computing optimal width for this l
    double * phi_arr = (double *) malloc(sizeof(double)*N);
    for (long i =0; i < N; i++) {
        phi_arr[i] = phi(i, jm);
    }
    double width = getWidth(l, km, phi_arr);
    free(phi_arr);
    double k2max = std::min(K, kmax+width);
    double step = width/((double ) npoints);
    arr->nk2 = k2max/step;
        
    // Allocating second array
    arr->k2 = (double *) malloc(sizeof(double) *arr->nk2);
    arr->Matk2 = (double *) malloc(sizeof(double) *arr->nk2*arr->Nt);

    double k2min = std::max(step/2.0, arr->k1[0]-width/2);
    
    #pragma omp parallel for schedule(guided)
    for(int64_t p=0; p < arr->nk2; p++) {
        arr->k2[p] = step*p+k2min;
        for(int64_t q = 0; q < arr->Nt; q++) {
            arr->Matk2[p*arr->Nt + q] = bm::sph_bessel<double>(l, qln[N*l + q]/K * arr->k2[p]);
        }
    }
    
    // Computing linear interpolation arrays
    arr->k_interp_sup = (long *) malloc(arr->nk2*sizeof(long));
    arr->k_interp_ratio = (double *) malloc(arr->nk2*sizeof(double));
    
    long ksup = 1;
    for (long i = 0; i < arr->nk2 ; i++) {
        while(arr->k2[i] > qln[N*l + ksup]/rmax && ksup < N - 1 ){
	  ksup++;
	}
	arr->k_interp_sup[i] = ksup;

	arr->k_interp_ratio[i] = (arr->k2[i] - qln[N*l + ksup - 1]/rmax)/(qln[N*l + ksup]/rmax - qln[N*l + ksup - 1]/rmax);
    }
    
    // Inserting the precomputed arrays in the buffer
    buffer[l] = arr;

    return arr;
}

double BesselWindow::evalWindow(double k, void* params)
{
    eval_params *par = (eval_params* ) params;
    
    double tempValue = 0;
    
    # pragma omp parallel for schedule(guided) reduction(+:tempValue)
    for (int64_t p=0; p < par->N; ++p) {
        tempValue += par->matrix[p] * bm::sph_bessel<double>(par->l, par->qln[par->N*par->l + p]*k/par->K);
    }
    
    return par->k0*tempValue - par->a;    
}


double BesselWindow::getWidth(int64_t l, double k, double *phi) {
    gsl_function F;
    
    double normFactor = sqrt(M_PI/pow(K,3));
    
    // Computes the transformation matrix
    double * matrix = (double *) malloc(N*sizeof(double));
    double tempValue;
    
    #pragma omp parallel for private(tempValue), schedule(guided)
    for (int64_t p=0; p < N ; ++p) {
        tempValue  = normFactor/bm::sph_bessel<double>(l+1, qln[N*l + p]);
        matrix[p] = tempValue * tempValue * bm::sph_bessel<double>(l, qln[N*l + p]*k/K)*phi[p];
    }
    
    struct eval_params params = {l, k, phi, matrix, 0, N, K, qln};
    
    double value = evalWindow(k,&params);
    
    params.a = value * precision;

    F.function = &evalWindow;
    F.params = &params;

    gsl_root_fsolver_set(s, &F, k, K);

    int status;
    int iter =0;
    double x_lo;
    double x_hi;
    double x_mean;
    double r;
    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        x_lo = gsl_root_fsolver_x_lower(s);
        x_hi = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(x_lo,x_hi,0,0.01);
        if( status == GSL_SUCCESS) {
            x_mean = 0.5*(x_hi + k);
            
            if(evalWindow(x_mean,&params) < 0) {
                gsl_root_fsolver_set(s, &F, k, x_mean);
                status = GSL_CONTINUE;
            }
        }
    } while (status == GSL_CONTINUE);
    
    return 2.0*(x_hi - k);
    
    free(matrix);
}


void BesselWindow::loadBesselZeroes(const char* filename) {

    int Nmax=0, Lmax=0;
    double val=0;
    bool fileOk = false;
    // try to open the file
    std::ifstream ifile(filename, std::ifstream::binary);
    if(ifile.good()) {
        // Let's try to read the number of zeros
        fileOk = true;
        ifile.read((char *) &Nmax,sizeof(int));
        ifile.read((char *) &Lmax,sizeof(int));
        if(Lmax >= L && Nmax >= N) {
            // Let's read all we can from this file
            for(int l=0; l <= Lmax; l++) {
                for(int p=0; p <Nmax; p++) {
                    if(! ifile.good()) {
                        fileOk = false;
                        break;
                    }
                    ifile.read((char *) &val,sizeof(double));
                    if( l <= L && p < N) qln[N*l + p] = val;
                }
            }
        } else {
            fileOk = false;
        }
    }
    // close the file
    ifile.close();

    if(fileOk == false) {
        std::cout << "Computing zeros of Bessel functions..."<< std::endl;

        #pragma omp parallel for schedule(guided)
        for(int l=0; l <= L; l++) {
            std::vector<long double> roots;
            boost::math::cyl_bessel_j_zero((double) (l+0.5), 1, N, std::back_inserter(roots));
            for(int p=0; p <N; p++) {
                qln[N*l + p] = roots[p];
            }
        }
        std::cout << "Done !" << std::endl;

        // Now we need to save this table so that it can be reused later.
        std::ofstream ofile(filename,  std::ofstream::binary | std::ofstream::trunc);
        if(ofile.good()) {
            Nmax = N;
            Lmax = L;
            ofile.write((char*) &Nmax, sizeof(int));
            ofile.write((char*) &Lmax, sizeof(int));
            for(int l=0; l <= Lmax; l++) {
                for(int p=0; p <Nmax; p++) {
                    val =  qln[N*l + p];
                    ofile.write((char *) &val, sizeof(double));
                }
            }
            ofile.close();
        } else {
            std::cout<< "Error, could not save bessel zeros in file " << filename << std::endl;
        }

    }
}