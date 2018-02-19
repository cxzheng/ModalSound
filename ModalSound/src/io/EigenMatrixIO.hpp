#ifndef EIGEN_MATRIX_IO_INC
#   define EIGEN_MATRIX_IO_INC

#ifdef MKL_DOMAIN_BLAS
#define EIGEN_MKL_DOMAIN_BLAS MKL_DOMAIN_BLAS
#else
#define EIGEN_MKL_DOMAIN_BLAS MKL_BLAS
#endif
#define EIGEN_MKL_DOMAIN_BLAS MKL_BLAS

#include <stdio.h>
#include <stdint.h>
#include "io/IOEndian.h"
#include <Eigen/Dense>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/static_assert.hpp>
#include "utils/macros.h"

namespace internal 
{

template<typename _Scalar, typename _EigenMatrix>
struct EigenMatrixIOImpl
{
    BOOST_STATIC_ASSERT((boost::is_floating_point<_Scalar>::value));

    static int write_ma_file(const char* filename, const _EigenMatrix* mat)
    {
        fprintf(stderr, "Only floating point matrix is supported now\n");
        return 0;
    }

    static int write_ma_file(FILE* fd, const _EigenMatrix* mat)
    {
        fprintf(stderr, "Only floating point matrix is supported now\n");
        return 0;
    }

    static int load_ma_file(const char* filename, _EigenMatrix* mat)
    {
        fprintf(stderr, "Only floating point matrix is supported now\n");
        return 0;
    }

    static int load_ma_file(FILE* fd, _EigenMatrix* mat)
    {
        fprintf(stderr, "Only floating point matrix is supported now\n");
        return 0;
    }
};

static const char MA_MAGIC_STR[] = "MA";

/*
 * write matrix with Double precision
 */
template<typename _EigenMatrix>
struct EigenMatrixIOImpl<double, _EigenMatrix>
{
    static int write_ma_file(const char* filename, const _EigenMatrix* mat)
    {
#ifdef IO_LITTLE_ENDIAN
        char  mtype[] = "\x9\x12";
#else
        char  mtype[] = "\x9\x21";
#endif
        int   nds[] = { 2, (int)mat->rows(), (int)mat->cols() };
        int   nele = nds[1] * nds[2];
        FILE* fd;

        fd = fopen(filename, "wb");
        if ( !fd ) return ERROR_RETURN;
        if ( fwrite(MA_MAGIC_STR, sizeof(char), 2, fd) != 2 )   // magic str
            goto ERR_RET;
        if ( fwrite(mtype, sizeof(char), 2, fd) != 2 )          // data type + endianess
            goto ERR_RET;
        if ( fwrite(nds, sizeof(int), 3, fd) != 3 )             // 2D data + nrow + ncolumn
            goto ERR_RET;

        if ( fwrite((const void*)mat->data(), sizeof(double), nele, fd) != nele )
            goto ERR_RET;

        fclose(fd);
        return SUCC_RETURN;

ERR_RET:
        fclose(fd);
        return -2;
    }

    static int write_ma_file(FILE* fd, const _EigenMatrix* mat)
    {
        if ( !fd ) return ERROR_RETURN;

#ifdef IO_LITTLE_ENDIAN
        char  mtype[] = "\x9\x12";
#else
        char  mtype[] = "\x9\x21";
#endif
        int   nds[] = { 2, (int)mat->rows(), (int)mat->cols() };
        int   nele = nds[1] * nds[2];

        if ( fwrite(MA_MAGIC_STR, sizeof(char), 2, fd) != 2 )   // magic str
            goto ERR_RET;
        if ( fwrite(mtype, sizeof(char), 2, fd) != 2 )          // data type + endianess
            goto ERR_RET;
        if ( fwrite(nds, sizeof(int), 3, fd) != 3 )             // 2D data + nrow + ncolumn
            goto ERR_RET;

        if ( fwrite((const void*)mat->data(), sizeof(double), nele, fd) != nele )
            goto ERR_RET;

        return SUCC_RETURN;

ERR_RET:
        return -2;
    }

    static int load_ma_file(const char* filename, _EigenMatrix* mat)
    {
        char magic[4];
        int  vals[3];
        int  endian, nv;
        double* ptr;

        assert(mat);
        FILE* fd = fopen(filename, "rb");
        if ( !fd ) return ERROR_RETURN;
        fread((void *)magic, sizeof(char), 4, fd);
        if ( magic[0] != 'M' || magic[1] != 'A' )
        {
            fprintf(stderr, "Wrong file format; (.MA) is expected\n");
            goto ERR_RET;
        }
        if ( (int)magic[2] != 9 )   // double precision
        {
            fprintf(stderr, "Wrong matrix precision; double is expected\n");
            goto ERR_RET;
        }
        if ( (int)magic[3] != 18 && (int)magic[3] != 33 ) 
        {
            fprintf(stderr, "Wrong endianess value\n");
            goto ERR_RET;
        }
        endian = magic[3] == 18 ? 0 : 1;
        fread((void *)vals, sizeof(int), 3, fd);
        vals[0] = !endian ? (int)le32toh((uint32_t)vals[0]) :
                            (int)be32toh((uint32_t)vals[0]);
        if ( vals[0] != 2 ) 
        {
            fprintf(stderr, "2D matrix data is expected\n");
            goto ERR_RET;
        }
        vals[1] = !endian ? (int)le32toh((uint32_t)vals[1]) :    // # of rows
                            (int)be32toh((uint32_t)vals[1]);
        vals[2] = !endian ? (int)le32toh((uint32_t)vals[2]) :    // # of columns
                            (int)be32toh((uint32_t)vals[2]);
        mat->resize(vals[1], vals[2]);
        ptr = mat->data();
        nv = vals[1] * vals[2];
        if ( fread((void *)ptr, sizeof(double), nv, fd) != nv ) 
        {
            fprintf(stderr, "Cannot read enough data\n");
            goto ERR_RET;
        }

#ifdef  IO_BIG_ENDIAN
        if ( !endian )
        {
            for(int i = 0;i < nv;++ i) 
                ptr[i] = (double)le64toh((uint64_t)ptr[i]);
        }
#else   // little-endian
        if ( endian )
        {
            for(int i = 0;i < nv;++ i) 
                ptr[i] = (double)be64toh((uint64_t)ptr[i]);
        }
#endif
        fclose(fd);
        return SUCC_RETURN;

ERR_RET:
        fclose(fd);
        return -2;
    }

    static int load_ma_file(FILE* fd, _EigenMatrix* mat)
    {
        if ( !fd ) return ERROR_RETURN;

        char magic[4];
        int  vals[3];
        int  endian, nv;
        double* ptr;

        assert(mat);
        fread((void *)magic, sizeof(char), 4, fd);
        if ( magic[0] != 'M' || magic[1] != 'A' )
        {
            fprintf(stderr, "Wrong file format; (.MA) is expected\n");
            goto ERR_RET;
        }
        if ( (int)magic[2] != 9 )   // double precision
        {
            fprintf(stderr, "Wrong matrix precision; double is expected\n");
            goto ERR_RET;
        }
        if ( (int)magic[3] != 18 && (int)magic[3] != 33 ) 
        {
            fprintf(stderr, "Wrong endianess value\n");
            goto ERR_RET;
        }
        endian = magic[3] == 18 ? 0 : 1;
        fread((void *)vals, sizeof(int), 3, fd);
        vals[0] = !endian ? (int)le32toh((uint32_t)vals[0]) :
                            (int)be32toh((uint32_t)vals[0]);
        if ( vals[0] != 2 ) 
        {
            fprintf(stderr, "2D matrix data is expected\n");
            goto ERR_RET;
        }
        vals[1] = !endian ? (int)le32toh((uint32_t)vals[1]) :    // # of rows
                            (int)be32toh((uint32_t)vals[1]);
        vals[2] = !endian ? (int)le32toh((uint32_t)vals[2]) :    // # of columns
                            (int)be32toh((uint32_t)vals[2]);
        mat->resize(vals[1], vals[2]);
        ptr = mat->data();
        nv = vals[1] * vals[2];
        if ( fread((void *)ptr, sizeof(double), nv, fd) != nv ) 
        {
            fprintf(stderr, "Cannot read enough data\n");
            goto ERR_RET;
        }

#ifdef  IO_BIG_ENDIAN
        if ( !endian )
        {
            for(int i = 0;i < nv;++ i) 
                ptr[i] = (double)le64toh((uint64_t)ptr[i]);
        }
#else   // little-endian
        if ( endian )
        {
            for(int i = 0;i < nv;++ i) 
                ptr[i] = (double)be64toh((uint64_t)ptr[i]);
        }
#endif
        return SUCC_RETURN;

ERR_RET:
        return -2;
    }
};

/*
 * write matrix with float precision (single precision)
 */
template<typename _EigenMatrix>
struct EigenMatrixIOImpl<float, _EigenMatrix>
{
    static int write_ma_file(const char* filename, const _EigenMatrix* mat)
    {
#ifdef IO_LITTLE_ENDIAN
        char  mtype[] = "\x7\x12";                              // 0x7: float
#else
        char  mtype[] = "\x7\x21";                              // 0x7: float
#endif
        int   nds[] = { 2, (int)mat->rows(), (int)mat->cols() };
        int   nele = nds[1] * nds[2];
        FILE* fd;

        fd = fopen(filename, "wb");
        if ( !fd ) return ERROR_RETURN;
        if ( fwrite(MA_MAGIC_STR, sizeof(char), 2, fd) != 2 )   // magic str
            goto ERR_RET;
        if ( fwrite(mtype, sizeof(char), 2, fd) != 2 )          // data type + endianess
            goto ERR_RET;
        if ( fwrite(nds, sizeof(int), 3, fd) != 3 )             // 2D data + nrow + ncolumn
            goto ERR_RET;

        if ( fwrite((const void*)mat->data(), sizeof(float), nele, fd) != nele )
            goto ERR_RET;

        fclose(fd);
        return SUCC_RETURN;

ERR_RET:
        fclose(fd);
        return -2;
    }

    static int write_ma_file(FILE* fd, const _EigenMatrix* mat)
    {
        if ( !fd ) return ERROR_RETURN;
#ifdef IO_LITTLE_ENDIAN
        char  mtype[] = "\x7\x12";                              // 0x7: float
#else
        char  mtype[] = "\x7\x21";                              // 0x7: float
#endif
        int   nds[] = { 2, (int)mat->rows(), (int)mat->cols() };
        int   nele = nds[1] * nds[2];

        if ( fwrite(MA_MAGIC_STR, sizeof(char), 2, fd) != 2 )   // magic str
            goto ERR_RET;
        if ( fwrite(mtype, sizeof(char), 2, fd) != 2 )          // data type + endianess
            goto ERR_RET;
        if ( fwrite(nds, sizeof(int), 3, fd) != 3 )             // 2D data + nrow + ncolumn
            goto ERR_RET;

        if ( fwrite((const void*)mat->data(), sizeof(float), nele, fd) != nele )
            goto ERR_RET;

        return SUCC_RETURN;
ERR_RET:
        return -2;
    }

    static int load_ma_file(const char* filename, _EigenMatrix* mat)
    {
        char magic[4];
        int  vals[3];
        int  endian, nv;
        float* ptr;

        assert(mat);
        FILE* fd = fopen(filename, "rb");
        if ( !fd ) return ERROR_RETURN;
        fread((void *)magic, sizeof(char), 4, fd);
        if ( magic[0] != 'M' || magic[1] != 'A' )
        {
            fprintf(stderr, "Wrong file format; (.MA) is expected\n");
            goto ERR_RET;
        }
        if ( (int)magic[2] != 7 )   // double precision
        {
            fprintf(stderr, "Wrong matrix precision; float is expected\n");
            goto ERR_RET;
        }
        if ( (int)magic[3] != 18 && (int)magic[3] != 33 ) 
        {
            fprintf(stderr, "Wrong endianess value\n");
            goto ERR_RET;
        }
        endian = magic[3] == 18 ? 0 : 1;
        fread((void *)vals, sizeof(int), 3, fd);
        vals[0] = !endian ? (int)le32toh((uint32_t)vals[0]) :
                            (int)be32toh((uint32_t)vals[0]);
        if ( vals[0] != 2 ) 
        {
            fprintf(stderr, "2D matrix data is expected\n");
            goto ERR_RET;
        }
        vals[1] = !endian ? (int)le32toh((uint32_t)vals[1]) :    // # of rows
                            (int)be32toh((uint32_t)vals[1]);
        vals[2] = !endian ? (int)le32toh((uint32_t)vals[2]) :    // # of columns
                            (int)be32toh((uint32_t)vals[2]);
        mat->resize(vals[1], vals[2]);
        ptr = mat->data();
        nv = vals[1] * vals[2];
        if ( fread((void *)ptr, sizeof(float), nv, fd) != nv ) 
        {
            fprintf(stderr, "Cannot read enough data\n");
            goto ERR_RET;
        }

#ifdef  IO_BIG_ENDIAN
        if ( !endian ) 
        {
            for(int i = 0;i < nv;++ i) 
                ptr[i] = (float)le32toh((uint32_t)ptr[i]);
        }
#else
        if ( endian )
        {
            for(int i = 0;i < nv;++ i) 
                ptr[i] = (float)be32toh((uint32_t)ptr[i]);
        }
#endif
        fclose(fd);
        return SUCC_RETURN;

ERR_RET:
        fclose(fd);
        return -2;
    }
    
    static int load_ma_file(FILE* fd, _EigenMatrix* mat)
    {
        if ( !fd ) return ERROR_RETURN;

        char magic[4];
        int  vals[3];
        int  endian, nv;
        float* ptr;

        assert(mat);
        fread((void *)magic, sizeof(char), 4, fd);
        if ( magic[0] != 'M' || magic[1] != 'A' )
        {
            fprintf(stderr, "Wrong file format; (.MA) is expected\n");
            goto ERR_RET;
        }
        if ( (int)magic[2] != 7 )   // double precision
        {
            fprintf(stderr, "Wrong matrix precision; float is expected\n");
            goto ERR_RET;
        }
        if ( (int)magic[3] != 18 && (int)magic[3] != 33 ) 
        {
            fprintf(stderr, "Wrong endianess value\n");
            goto ERR_RET;
        }
        endian = magic[3] == 18 ? 0 : 1;
        fread((void *)vals, sizeof(int), 3, fd);
        vals[0] = !endian ? (int)le32toh((uint32_t)vals[0]) :
                            (int)be32toh((uint32_t)vals[0]);
        if ( vals[0] != 2 ) 
        {
            fprintf(stderr, "2D matrix data is expected\n");
            goto ERR_RET;
        }
        vals[1] = !endian ? (int)le32toh((uint32_t)vals[1]) :    // # of rows
                            (int)be32toh((uint32_t)vals[1]);
        vals[2] = !endian ? (int)le32toh((uint32_t)vals[2]) :    // # of columns
                            (int)be32toh((uint32_t)vals[2]);
        mat->resize(vals[1], vals[2]);
        ptr = mat->data();
        nv = vals[1] * vals[2];
        if ( fread((void *)ptr, sizeof(float), nv, fd) != nv ) 
        {
            fprintf(stderr, "Cannot read enough data\n");
            goto ERR_RET;
        }

#ifdef  IO_BIG_ENDIAN
        if ( !endian ) 
        {
            for(int i = 0;i < nv;++ i) 
                ptr[i] = (float)le32toh((uint32_t)ptr[i]);
        }
#else
        if ( endian )
        {
            for(int i = 0;i < nv;++ i) 
                ptr[i] = (float)be32toh((uint32_t)ptr[i]);
        }
#endif
        return SUCC_RETURN;
ERR_RET:
        return -2;
    }
};

}

// --------------------------------------------------------------------------------------

template<class EigenMatrix>
int write_ma_eigen_matrix(const char* filename, const EigenMatrix* mat)
{
    return internal::EigenMatrixIOImpl<
        typename EigenMatrix::Scalar,
        EigenMatrix>::write_ma_file(filename, mat);
}

template<class EigenMatrix>
int write_ma_eigen_matrix(FILE* fd, const EigenMatrix* mat)
{
    return internal::EigenMatrixIOImpl<
        typename EigenMatrix::Scalar,
        EigenMatrix>::write_ma_file(fd, mat);
}

template<class EigenMatrix>
int read_ma_eigen_matrix(const char* filename, EigenMatrix* mat)
{
    return internal::EigenMatrixIOImpl<
        typename EigenMatrix::Scalar,
        EigenMatrix>::load_ma_file(filename, mat);
}

template<class EigenMatrix>
int read_ma_eigen_matrix(FILE* fd, EigenMatrix* mat)
{
    return internal::EigenMatrixIOImpl<
        typename EigenMatrix::Scalar,
        EigenMatrix>::load_ma_file(fd, mat);
}
#endif
