/*
 * =====================================================================================
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * -------------------------------------------------------------------------------------
 *
 *       Filename:  MatrixIO.cpp
 *
 *        Version:  1.0
 *        Created:  03/04/11 14:23:34
 *       Revision:  none
 *       Compiler:  gcc/intel compiler
 *
 *         Author:  Changxi Zheng (cz), cxzheng@cs.cornell.edu
 *                  Cornell University
 *
 * =====================================================================================
 */
#include "MatrixIO.hpp"
#include "utils/macros.h"

#ifndef __APPLE__

#include <endian.h>
#if !defined(le32toh) || !defined(be32toh) || !defined(le64toh) || !defined(be64toh)
# include <byteswap.h>
#endif 

#ifndef le32toh
#   if __BYTE_ORDER == __LITTLE_ENDIAN
#       define le32toh(x) (x)
#   else
#       define le32toh(x) __bswap_32 (x)
#   endif
#endif

#ifndef be32toh
#   if __BYTE_ORDER == __LITTLE_ENDIAN 
#       define be32toh(x) __bswap_32 (x)
#   else
#       define be32toh(x) (x)
#   endif       
#endif

#ifndef le64toh
#   if __BYTE_ORDER == __LITTLE_ENDIAN 
#       define le64toh(x) (x)
#   else
#       define le64toh(x) __bswap_64 (x)
#   endif       
#endif

#ifndef be64toh
#   if __BYTE_ORDER == __LITTLE_ENDIAN 
#       define be64toh(x) __bswap_64 (x)
#   else
#       define be64toh(x) (x)
#   endif       
#endif

#else   // fix Missing byteswap.h and endian.h on Mac OS X

#include <machine/endian.h>
#include <libkern/OSByteOrder.h>

#define htobe16(x) OSSwapHostToBigInt16(x)
#define htole16(x) OSSwapHostToLittleInt16(x)
#define be16toh(x) OSSwapBigToHostInt16(x)
#define le16toh(x) OSSwapLittleToHostInt16(x)

#define htobe32(x) OSSwapHostToBigInt32(x)
#define htole32(x) OSSwapHostToLittleInt32(x)
#define be32toh(x) OSSwapBigToHostInt32(x)
#define le32toh(x) OSSwapLittleToHostInt32(x)

#define htobe64(x) OSSwapHostToBigInt64(x)
#define htole64(x) OSSwapHostToLittleInt64(x)
#define be64toh(x) OSSwapBigToHostInt64(x)
#define le64toh(x) OSSwapLittleToHostInt64(x)

#define __BIG_ENDIAN    BIG_ENDIAN
#define __LITTLE_ENDIAN LITTLE_ENDIAN
#define __BYTE_ORDER    BYTE_ORDER

#endif


Matrix<double>* load_ma_matrixd(const char* file)
{
    char magic[4];
    int  vals[3];
    int  endian, nv;
    Matrix<double>* mret = NULL;
    double* ptr;

    FILE* fd = fopen(file, "rb");
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
    mret = new Matrix<double>(vals[1], vals[2]);
    nv = vals[1] * vals[2];
    ptr = mret->data();
    if ( fread((void *)ptr, sizeof(double), nv, fd) != nv ) 
    {
        fprintf(stderr, "Cannot read enough data\n");
        delete mret;
        goto ERR_RET;
    }

    if ( !endian && __BYTE_ORDER != __LITTLE_ENDIAN ) 
    {
        for(int i = 0;i < nv;++ i) 
                ptr[i] = (double)le64toh((uint64_t)ptr[i]);
    }
    else if ( endian && __BYTE_ORDER == __LITTLE_ENDIAN) 
    {
        for(int i = 0;i < nv;++ i) 
            ptr[i] = (double)be64toh((uint64_t)ptr[i]);
    }

    fclose(fd);
    return mret;

ERR_RET:
    fclose(fd);
    return NULL;
}

int write_ma_matrixd(const char* filename, const Matrix<double>* mat)
{
    char  mtype[] = "\x9\x12";
    int   nds[] ={ 2, (int)mat->shape()[0], (int)mat->shape()[1] };
    int   nele = nds[1] * nds[2];
    FILE* fd;

    if ( __BYTE_ORDER == __BIG_ENDIAN ) mtype[1] = '\x21';

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

int write_ma_matrixd(const char* filename, int m, int n, const double* ele)
{
    char  mtype[] = "\x9\x12";
    int   nds[] ={ 2, m, n };
    int   nele = nds[1] * nds[2];
    FILE* fd;

    if ( __BYTE_ORDER == __BIG_ENDIAN ) mtype[1] = '\x21';

    fd = fopen(filename, "wb");
    if ( !fd ) return ERROR_RETURN;
    if ( fwrite(MA_MAGIC_STR, sizeof(char), 2, fd) != 2 )   // magic str
        goto ERR_RET;
    if ( fwrite(mtype, sizeof(char), 2, fd) != 2 )          // data type + endianess
        goto ERR_RET;
    if ( fwrite(nds, sizeof(int), 3, fd) != 3 )             // 2D data + nrow + ncolumn
        goto ERR_RET;

    if ( fwrite((const void*)ele, sizeof(double), nele, fd) != nele )
        goto ERR_RET;
    
    fclose(fd);
    return SUCC_RETURN;

ERR_RET:
    fclose(fd);
    return -2;
}

