#ifndef MATRIX_IO_HPP
#   define MATRIX_IO_HPP

#include <iostream>
#include <fstream>
#include <stdint.h>
#include "utils/macros.h"
#include "sc/Matrix.hpp"
#include "sc/PardisoMatrix.hpp"

struct PardisoMatrixIO
{
    /* write as matlab sparse matrix 
     * the format is like:
     * m n 0
     * r1 c1 v1
     * r2 c2 v2
     * ...
     * rn cn vn
     */
    template <typename T>
    static void write_matlab_sparse(const PardisoMatrix<T>& mat, std::ostream& out)
    {
        out << mat.num_rows() << ' ' << mat.num_cols() << " 0" << std::endl;
        for(int i = 0;i < (int)mat.rows_.size();++ i)
        {
            typename std::map<int, T*>::const_iterator end = mat.rows_[i].end();
            for(typename std::map<int, T*>::const_iterator it = mat.rows_[i].begin();
                    it != end;++ it)
            {
                out << i+1 << ' ' << it->first+1 << ' ' << *(it->second) << std::endl;
                if ( mat.is_symmetric() && i != it->first )
                    out << it->first+1 << ' ' << i+1 << ' ' << *(it->second) << std::endl;
            }
        }
    }

    /*
     * binary file format for the Matlab sparse matrix structure
     * <unsigned char> :: 0: float 1: double 2: complex float 3: complex double
     * <int> <int> <int>:: height x width ; # of nonzeros
     * <int> <int> <value>:: row index, col index, value
     */
    static int write_matlab_sparse_bin(const PardisoMatrix<double>& mat, const char* file)
    {
        std::ofstream fout(file, std::ios::binary);
        if ( !fout.good() )
        {
            std::cerr << "PardisoMatrixIO::Fail to open file: " << file << " to write" << std::endl;
            return ERROR_RETURN;
        }

        uint8_t label = 1;  // double
        fout.write((char *)&label, sizeof(uint8_t));

        // --------------------------------------------------------------

        int h = mat.num_rows(), w = mat.num_cols(); 
        int n = 0;
        for(int i = 0;i < (int)mat.rows_.size();++ i)
        {
            std::map<int, double*>::const_iterator end = mat.rows_[i].end();
            for(std::map<int, double*>::const_iterator it = mat.rows_[i].begin();
                    it != end;++ it)
            {
                ++ n;
                if ( mat.is_symmetric() && i != it->first ) ++ n;
            }
        }

        // --------------------------------------------------------------

        fout.write((char *)&h, sizeof(int));    // # of rows
        fout.write((char *)&w, sizeof(int));    // # of cols
        fout.write((char *)&n, sizeof(int));    // # of nonzeros

        for(int i = 0;i < (int)mat.rows_.size();++ i)
        {
            std::map<int, double*>::const_iterator end = mat.rows_[i].end();
            for(std::map<int, double*>::const_iterator it = mat.rows_[i].begin();
                    it != end;++ it)
            {
                h = i+1;            // row index
                w = it->first + 1;  // col index (1-based)
                fout.write((char *)&h, sizeof(int));    // row index
                fout.write((char *)&w, sizeof(int));    // col index
                fout.write((char *)(it->second), sizeof(double));
                if ( mat.is_symmetric() && i != it->first ) 
                {
                    fout.write((char *)&w, sizeof(int));    // row index
                    fout.write((char *)&h, sizeof(int));    // col index
                    fout.write((char *)(it->second), sizeof(double));
                }
            } // end for
        }

        fout.close();
        return SUCC_RETURN;
    }

    /*
     * My own CSC binary format for sparse matrix
     * <unsigned char> :: 0: float 1: double 2: complex float 3: complex double
     * <unsigned char> :: v & 0b00000001 != 0: symmetric ; v & 0b00000010 != 0: positive def. 
     * <int> <int> <int>:: height x width ; # of nonzeros
     * <int,int...> rowidx (1-based)
     * <int,int...> colidx (1-based)
     * <...> data
     */
    static int write_csc_format(const PardisoMatrix<double>& mat, const char* file)
    {
        std::ofstream fout(file, std::ios::binary);
        if ( !fout.good() )
        {
            std::cerr << "PardisoMatrixIO::Fail to open file: " << file << " to write" << std::endl;
            return ERROR_RETURN;
        }

        uint8_t label = 1;  // double
        fout.write((char *)&label, sizeof(uint8_t));

        label = mat.is_symmetric() ? (uint8_t)1 : (uint8_t)0;
        if ( mat.is_positive_definite() ) label |= 2;
        fout.write((char *)&label, sizeof(uint8_t));

        int h = mat.num_rows(), w = mat.num_cols(), n = mat.num_nonzeros();
        fout.write((char *)&h, sizeof(int));
        fout.write((char *)&w, sizeof(int));
        fout.write((char *)&n, sizeof(int));

        const int* rowidx = mat.row_indices();
        const int* colidx = mat.col_indices();
        const double* data = mat.data();
        fout.write((const char *)rowidx, sizeof(int)*(h+1));
        fout.write((const char *)colidx, sizeof(int)*n);
        fout.write((const char *)data, sizeof(double)*n);

        fout.close();
        return SUCC_RETURN;
    }
};

// ======================================================================================

static const char MA_MAGIC_STR[] = "MA";
/* 
 * Load 2D matrix from file in double
 *
 * 2D matrix file format:
 * first two bytes: 'MA'
 * 3rd byte:        [data precision type] 9: double
 * 4th byte:        [endianess]
 */
Matrix<double>* load_ma_matrixd(const char* file);
int             write_ma_matrixd(const char* filename, const Matrix<double>* mat);
int             write_ma_matrixd(const char* filename, int m, int n, const double* ele);

// ======================================================================================

#endif
