/******************************************************************************
 * Eigen-solver running on a local computer
 *
 * modal-eigen <K-mat-file> <M-mat-file> <fmin> <fmax> <# eig> <outfile>
 * It reads the K matrix and M matrix in given files, solve the generalized 
 * eigen-value problem, Kv = lambda Mv, and store the eigenvector and eigenvalues
 * into the given file.
 *****************************************************************************/
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <boost/program_options.hpp>
#ifdef USE_MKL
#   include <mkl.h>
#else
#   include "feast.h"
#   include "feast_sparse.h"
#endif
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

#include "utils/math.hpp"
#include "utils/term_msg.h"
#include "utils/nano_timer.h"

using namespace std;

static string stiffMFile, massMFile, outFile;
static int    numEigv = 200;
static double density, freqLow, freqHigh;
static bool   verbose = false;

static vector<int> ptrrow[2];
static vector<int> idxcol[2];
static vector<double> data[2];

/*
 * The matrix file only stores the upper triangle part
 */
static uint8_t read_csc_dmatrix(const char* file, 
        std::vector<int>& ptrrow,
        std::vector<int>& idxcol,
        std::vector<double>& data,
        int& nrow, int& ncol)
{
    if ( verbose ) 
    {
        PRINT_MSG("Load CSC matrix file: %s\n", file);
    }

    ifstream fin(file, ios::binary);
    if ( fin.fail() )
    {
        PRINT_ERROR("read_csc_dmatrix:: Cannot open file [%s] to read\n", file);
        return 255;
    }

    uint8_t ret;
    fin.read((char *)&ret, sizeof(uint8_t));
    if ( ret != 1 )
    {
        PRINT_ERROR("read_csc_dmatrix:: The matrix data should be in double format\n");
        return 255;
    }

    int n;
    fin.read((char *)&ret, sizeof(uint8_t));
    fin.read((char *)&nrow, sizeof(int));
    fin.read((char *)&ncol, sizeof(int));
    fin.read((char *)&n,    sizeof(int));

    if ( (ret & 1) && (nrow != ncol) ) // symmetric
    {
        PRINT_ERROR("read_csc_dmatrix:: Symmetric matrix should be square\n");
        return 255;
    }

    ptrrow.resize(nrow+1);
    idxcol.resize(n);
    data.resize(n);
    fin.read((char *)(ptrrow.data()), sizeof(int)*(nrow+1));
    fin.read((char *)(idxcol.data()), sizeof(int)*n);
    fin.read((char *)(data.data()),   sizeof(double)*n);

    fin.close();
    return ret;
}

static void parse_cmd(int argc, char* argv[])
{
    namespace po = boost::program_options;
    po::options_description genericOpt("Generic options");
    genericOpt.add_options()
            ("help,h", "display help information");
    po::options_description configOpt("Configuration");
    configOpt.add_options()
            ("neig,n", po::value<int>(&numEigv)->default_value(200), 
                    "Maximum number of smallest eigenvalues to compute")
            ("stiff,s", po::value<string>(&stiffMFile)->default_value(""),
                    "Name of the stiffness matrix file")
            ("mass,m", po::value<string>(&massMFile)->default_value(""),
                    "Name of the mass matrix file")
            ("out,o", po::value<string>(&outFile)->default_value(""),
                    "Name of the output modes file")
            ("density,d", po::value<double>(&density)->default_value(1.),
                    "Reference density value")
            ("fmin", po::value<double>(&freqLow)->default_value(5.),
                    "Lowest frequency value based on the estimated density")
            ("fmax", po::value<double>(&freqHigh)->default_value(15000.),
                    "Highest frequency value based on the estimated density")
            ("verbose,v", "Display details");
    // use configure file to specify the option
    po::options_description cfileOpt("Configure file");
    cfileOpt.add_options()
            ("cfg-file", po::value<string>(), "configuration file");

    po::options_description cmdOpts;
    cmdOpts.add(genericOpt).add(configOpt).add(cfileOpt);

    po::variables_map vm;
    store(po::parse_command_line(argc, argv, cmdOpts), vm);
    if ( vm.count("cfg-file") )
    {
        ifstream ifs(vm["cfg-file"].as<string>().c_str());
        store(parse_config_file(ifs, configOpt), vm);
    }
    po::notify(vm);

    if ( vm.count("help") )
    {
        printf("Usage: %s [options] \n", argv[0]);
        printf("       This executable takes as input the stiffness and mass matrices\n");
        printf("       of a tet. mesh, and computes the eigenvectors and eigenvalues\n");
        printf("       using the eigensolvers provided in Intel MKL\n");
        cout << cmdOpts << endl;
        exit(0);
    }
    verbose = vm.count("verbose") > 0;

    if ( massMFile.empty() )
    {
        PRINT_ERROR("Specify mass matrix file\n");
        exit(1);
    }

    if ( stiffMFile.empty() ) 
    {
        PRINT_ERROR("Specify stiffness matrix file\n");
        exit(1);
    }

    if ( density <= 0. )
    {
        PRINT_ERROR("Density value must be positive [d=%g now]\n", density);
        exit(1);
    }
    
    if ( freqLow >= freqHigh )
    {
        PRINT_ERROR("Frequency lower bound must be smaller than the upper bound\n");
        exit(1);
    }

    if ( outFile.empty() ) 
    {
        PRINT_ERROR("Specify the output file\n");
        exit(1);
    }

    if ( verbose )
    {
        PRINT_MSG("=============== Problem Summary ===============\n");
        PRINT_MSG("Mass Matrix:                %s\n", massMFile.c_str());
        PRINT_MSG("Stiffness Matrix:           %s\n", stiffMFile.c_str());
        PRINT_MSG("Output file:                %s\n", outFile.c_str());
        PRINT_MSG("# of eigenvalues est.:      %d\n", numEigv);
        PRINT_MSG("Reference density:          %g\n", density);
        PRINT_MSG("Low frequency:              %g\n", freqLow);
        PRINT_MSG("High frequency:             %g\n", freqHigh);
        PRINT_MSG("===============================================\n");
    }
}

struct EIGV_CMP_
{
    const double* eigv; // eigenvalue
    EIGV_CMP_(const double* ev):eigv(ev) {}

    inline bool operator ()(int a, int b)
    {   return eigv[a] < eigv[b]; }
};

struct _Parallel_Diag_M
{
    const int     N;
    const double* data;     // M matrix 
    const int*    ptrow;
    const int*    idxcol;
    const double* ev;       // un-normalized eigenvector
    double*       ret;      // return value (M scalars)
    double*       aux;      // aux to store M*v vectors

    _Parallel_Diag_M(int N, const double* d, const int* p, const int* i, const double* e, double* r, double* a):
                N(N), data(d), ptrow(p), idxcol(i), ev(e), ret(r), aux(a) { }

    void operator() (const tbb::blocked_range<int>& r) const
    {
        for(int vi = r.begin();vi != r.end();++ vi)
            eigenvec_quadratic(vi);
    }

    // compute v'*M*v
    void eigenvec_quadratic(int vi) const
    {
        const double* evi = ev + (vi*N);
        double* ptr = aux + (vi*N);

        for(int i = 0;i < N;++ i) ptr[i] = 0.;
        for(int i = 0;i < N;++ i)
        {
            for(int j = ptrow[i]-1;j < ptrow[i+1]-1;++ j)
            {
                int col = idxcol[j] - 1;
                ptr[i] += evi[col]*data[j];
                if ( col != i )
                    ptr[col] += evi[i]*data[j];
            }
        }

        ret[vi] = 0;
        for(int i = 0;i < N;++ i) ret[vi] += ptr[i]*evi[i];
    }
};

/*
 * File Format:
 * <int>: size of the eigen problem
 * <int>: # of eigenvalues
 * <eigenvalues>
 * <eigenvec_1>
 * ...
 * <eigenvec_n>
 */
static void write_eigenvalues(int nev, int nsz, 
        const int* ids, const double* eval, const double* evec, 
        const char* file)
{
    ofstream fout(file, std::ios::binary);
    if ( !fout.good() ) 
    {
        cerr << "write_eigenvalues::Cannot open file " << file << " to write" << endl;
        return;
    }

    // size of the eigen-problem. Here square matrix is assumed.
    fout.write((char *)&nsz, sizeof(int));
    fout.write((char *)&nev, sizeof(int));

    // output eigenvalues
    for(int vid = 0;vid < nev;++ vid)
    {
        fout.write((const char*)&eval[ids[vid]], sizeof(double));
        printf("ev#%3d:  %lf %lfHz\n", vid, eval[ids[vid]], sqrt(eval[ids[vid]]/density)*0.5*M_1_PI);
    }
    // output eigenvectors
    for(int vid = 0;vid < nev;++ vid)
        fout.write((const char*)&evec[ids[vid]*nsz], sizeof(double)*nsz);

    fout.close();
}

int main(int argc, char* argv[])
{
    int n = tbb::task_scheduler_init::default_num_threads();
    PRINT_MSG("%d threads are enabled\n", n);
    tbb::task_scheduler_init tbbInit;

    parse_cmd(argc, argv);

    //// load the data
    int nrowK, ncolK;
    uint8_t c = read_csc_dmatrix(stiffMFile.c_str(), ptrrow[0],
                            idxcol[0], data[0], nrowK, ncolK);
    if ( c == 255 )
    {
        PRINT_ERROR("Fail to read the matrix file: %s\n", stiffMFile.c_str());
        exit(1);
    }
    if ( !(c & 1) )
    {
        PRINT_ERROR("Symmetric matrix is excepted from file [%s]\n", stiffMFile.c_str());
        exit(1);
    }

    int nrowM, ncolM;
    c = read_csc_dmatrix(massMFile.c_str(), ptrrow[1], 
                         idxcol[1], data[1], nrowM, ncolM);
    if ( c == 255 )
    {
        PRINT_ERROR("Fail to read the matrix file: %s\n", massMFile.c_str());
        exit(1);
    }
    if ( !(c & 1) || !(c & 2) )
    {
        PRINT_ERROR("Symmetric positive definite matrix is excepted from file [%s]\n", massMFile.c_str());
        exit(1);
    }
    if ( nrowK != nrowM )
    {
        PRINT_ERROR("Two matrices should have the same size.\n");
        exit(1);
    }
    if ( numEigv <= 0 || numEigv > nrowK-2 )
    {
        PRINT_ERROR("number of eigenvalues is out of range: maximum=%d\n", nrowK-2);
        exit(1);
    }

    if ( verbose ) PRINT_MSG("Sparse matrix size: %d\n", nrowK);

    // solve the generalized eigen problem using FEAST
    const double eMin = M_SQR(freqLow * 2. * M_PI) * density;
    const double eMax = M_SQR(freqHigh * 2. * M_PI) * density;
    if ( verbose ) PRINT_MSG("Search frequency interval: [%g, %g]\n", eMin, eMax);
    int M, info = 0, M0 = numEigv;

    vector<double> eval(M0), evec(M0*nrowK), res(M0);

    MKL_INT feastparam[128]; 
    double  epsout;
    int     loop;
    const char UPLO = 'U';

    info = 3;
    feastinit(feastparam);
    if ( verbose ) feastparam[0]=1;  /*change from default value */

    dfeast_scsrgv(&UPLO, &nrowK, data[0].data(), ptrrow[0].data(), idxcol[0].data(),
                  data[1].data(), ptrrow[1].data(), idxcol[1].data(),
                  feastparam, &epsout, &loop, &eMin, &eMax, &M0, 
                  eval.data(), evec.data(), &M, res.data(), &info);

    switch (info)
    {
        case 0:
            break;
        case 1:
            PRINT_ERROR("Failed to run FEAST. [code = %d]\n", info);
            PRINT_ERROR("no eigenvalue found in the search interval\n");
            return 1;
        case 3:
            PRINT_ERROR("Failed to run FEAST. [code = %d]\n", info);
            PRINT_ERROR("Size of M0 (=%d) is too small\n", numEigv);
            return 1;
        default:
            PRINT_ERROR("Failed to run FEAST. [code = %d]\n", info);
            return 1;
    }

    if ( verbose ) PRINT_MSG("%d eigenvalues are found in the interval\n", M);

    vector<double> dv(M);
    vector<double> aux(M*nrowK);
    //// check weather or not the eigen vectors are normalized
    _Parallel_Diag_M pp(nrowK, data[1].data(), ptrrow[1].data(), idxcol[1].data(),
                        evec.data(), dv.data(), aux.data());
    if ( M < 10 )
        pp(tbb::blocked_range<int>(0, M));
    else
        tbb::parallel_for(tbb::blocked_range<int>(0,M), pp);

    for(int i = 0;i < M;++ i)
    {
        if ( fabs(dv[i] - 1.0) > 1E-8 )
        {
            if ( verbose ) PRINT_WARNING("normalize eigenvector %d: nv = %g\n", i, dv[i]);
            double ss = 1. / sqrt(dv[i]);
            cblas_dscal(nrowK, ss, &(evec[i*nrowK]), 1);
        }
    }

    //// sort the eigen pair in increasing order
    vector<int> sortids(M);
    for(int i = 0;i < M;++ i) sortids[i] = i;
    sort(sortids.begin(), sortids.end(), EIGV_CMP_(eval.data()));

    //// save the data
    if ( verbose ) PRINT_MSG("Write eigenvalues & eigenvectors into file: %s\n", outFile.c_str());
    write_eigenvalues(M, nrowK, sortids.data(), eval.data(), evec.data(), outFile.c_str());
    return 0;
}
