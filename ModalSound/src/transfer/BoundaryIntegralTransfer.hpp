#ifndef BOUNDARY_INTEGRAL_TRANSFER_HPP
#   define BOUNDARY_INTEGRAL_TRANSFER_HPP

#include <complex>
#include <vector>
#include <sstream>
#include <fstream>
#include <cctype>
#include <iomanip>

#include "utils/term_msg.h"
#include "utils/strings.hpp"
#include "geometry/Triangle.hpp"
#include "HelmholtzSolutions.hpp"

/*
 * Evaluate the transfer function value (P) using the boundary element 
 * integration
 * 
 * The input files are the input and output files for FastBEM solver. Then 
 * it evaluates the transfer value (P) at given location by computing the boundary 
 * integral.
 */
template <typename T>
class BoundaryIntegralTransfer
{
    public:
        typedef std::complex<T>     TComplex;

        BoundaryIntegralTransfer(const char* fileinput, const char* fileoutput);

        /* 
         * Kirchoff integral
         */
        virtual TComplex eval(const Point3<T>& pt) const
        {
            const T ONE_THIRD = 1. / 3.;
            TComplex ret = 0;
        
            for(int i = 0;i < tgl_.size();++ i)
            {
                // center point of a triangle
                Point3<T> c = (vtx_[tgl_[i][0]] + vtx_[tgl_[i][1]] + vtx_[tgl_[i][2]]) * ONE_THIRD;
                TComplex val, deri;
                // Green's function value and its derivative
                helmholtz_green_func_with_dir_deri(waveNum_, pt, c, tglNml_[i], val, deri);
                // Integral
                ret += iOmegaRho_ * v_[i] * val * tglArea_[i];  // i*omega*rho*vel * G(x,y) * s_i
                ret -= p_[i] * deri * tglArea_[i];
            }
            return ret;
        }

        inline T wave_number() const
        {   return waveNum_;    }

        inline size_t num_elements() const
        {   return tgl_.size(); }

    protected:
        std::vector< TComplex >     v_;             // normal velocity as boundary condition
        std::vector< TComplex >     p_;             // pressure value on boundary
                                                    // load from FastBEM output file

        //// rigid object mesh geometry
        std::vector< Point3<T> >    vtx_;           // nodal positions
        std::vector< Tuple3i >      tgl_;           // surface triangles
        std::vector< T >            tglArea_;
        std::vector< Vector3<T> >   tglNml_;
        std::vector< Point3d >      fieldPts_;      // output field points

        TComplex                    iOmegaRho_;
        T                           waveNum_;       // wave number
        T                           speed_;
        T                           density_;
        T                           freq_;
        T                           omega_;
};

///////////////////////////////////////////////////////////////////////////////

template <typename T>
BoundaryIntegralTransfer<T>::BoundaryIntegralTransfer(
        const char* fileinput, const char* fileoutput)  //.
{
    char text[1024];
    int ntgl, nvtx, nfp;

    // ==================================================================
    //// NOTE: assume no empty line
    std::ifstream fin(fileinput);
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot load file: %s\n", fileinput);
        exit(1);
    }

    fin.getline(text, 1024);            // first line
    fin.getline(text, 1024);            // 2nd line "Complete 1/3"
    fin.getline(text, 1024);            // 3rd line "Full 0 0.d0"
    fin.getline(text, 1024);            // 4th line "# B.E. # Nodes # F.P. # F.C."
    std::istringstream(text) >> ntgl >> nvtx >> nfp;

    fin.getline(text, 1024);            // 5th line
    int ia, ib;
    std::istringstream(text) >> ia >> ib;
    if ( ia || ib )
    {
        PRINT_ERROR("Unsupported problem setting at L. 5\n");
        exit(1);
    }
    fin.getline(text, 1024);            // 6th line
    std::istringstream(text) >> ia;
    if ( ia )
    {
        PRINT_ERROR("Unsupported problem setting at L. 6\n");
        exit(1);
    }

    fin.getline(text, 1024);            // 7th line: speed and density
    std::istringstream(text) >> speed_ >> density_;

    fin.getline(text, 1024);                // L.8
    fin.getline(text, 1024);                // L.9
    fin.getline(text, 1024);                // L.10
    if ( text[0] != '$' || text[2] != 'N' )  // Start to read number of nodes
    {
        PRINT_ERROR("incorrect file format at L. 7\n");
        exit(1);
    }

    //// load vertices
    int id;
    T xx, yy, zz;
    vtx_.resize(nvtx);
    for(int i = 0;i < nvtx;++ i)
    {
        fin.getline(text, 1024);
        std::istringstream(text) >> id >> xx >> yy >> zz;
        vtx_[id-1].set(xx, yy, zz);
    }

    //// load boundary conditions
    fin.getline(text, 1024);
    if ( text[0] != '$' || text[2] != 'E' )
    {
        PRINT_ERROR("incorrect file format for boundary condition\n");
        exit(1);
    }
    v_.resize(ntgl);
    tgl_.resize(ntgl);
    tglArea_.resize(ntgl);
    tglNml_.resize(ntgl);
    int idx, idy, idz;
    char c1, c2;
    T v1, v2;
    for(int i = 0;i < ntgl;++ i)
    {
        fin.getline(text, 1024);
        std::istringstream(text) >> id >> idx >> idy >> idz >> ia >> c1 >> v1 >> c2 >> v2;
        if (ia != 2)
        {
            PRINT_ERROR("Only Neumann boundary condition is supported now\n");
            exit(1);
        }

        v_[id-1] = TComplex(v1, v2);
        tgl_[id-1].set(idx-1, idy-1, idz-1);    // 1-based --> 0-based
        tglNml_[id-1]  = Triangle<T>::normal(vtx_[idx-1], vtx_[idy-1], vtx_[idz-1]);
        tglArea_[id-1] = tglNml_[id-1].normalize2() * 0.5;
    } // end for

    //// Load the field points (evaluation locations) if any
    fin.getline(text, 1024);
    if ( text[0] != '$' || text[2] != 'F' )
    {
        PRINT_ERROR("incorrect file format for boundary condition\n");
        exit(1);
    }
    fieldPts_.resize(nfp);
    for(int i = 0;i < nfp;++ i)
    {
        fin.getline(text, 1024);
        std::istringstream(text) >> id >> xx >> yy >> zz;
        fieldPts_[id-1].set(xx, yy, zz);
    }

    if ( fin.fail() ) 
    {
        PRINT_ERROR("unexpected ending of file\n");
        exit(1);
    }
    fin.close();

    // ==================================================================
    /////////////////// load the boundary pressure values
    fin.open(fileoutput);
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot load file: %s\n", fileoutput);
        exit(1);
    }

    p_.resize(ntgl);
    int ptr;

    do {
        fin.getline(text, 1024);
        if ( fin.fail() ) 
        {
            PRINT_ERROR("Cannot load boundary pressure values\n");
            exit(1);
        }
        // trim the leading space
        ptr = 0;
        while ( text[ptr] != '\n' && text[ptr] != '\0' && isspace(text[ptr]) ) ++ ptr;
    } while (text[ptr] != 'F' || text[ptr+1] != 'r' || text[ptr+2] != 'e' || text[ptr+3] != 'q'); 

    PRINT_MSG("read solution: %s\n", text);
    std::vector<std::string> tks = sploosh::tokenize(std::string(text));
    for(size_t i = 0;i < tks[6].length();++ i) if ( tks[6][i] == 'D' ) tks[6][i] = 'E';
    std::istringstream(tks[6]) >> freq_;

    //// go find the beginning of boundary element values
    do {
        fin.getline(text, 1024);
        if ( fin.fail() ) 
        {
            PRINT_ERROR("Cannot load boundary pressure values\n");
            exit(1);
        }

        // trim the leading space
        ptr = 0;
        while ( text[ptr] != '\n' && text[ptr] != '\0' && isspace(text[ptr]) ) ++ ptr;
    } while (text[ptr] != 'E' || text[ptr+1] != 'l' || text[ptr+2] != 'e' );

    //// load the boundary pressure values
    for(int i = 0;i < ntgl;++ i)
    {
        fin.getline(text, 1024);
        std::istringstream(text) >> id >> c1 >> v1 >> c2 >> v2;
        p_[id-1] = TComplex(v1, v2);    // load boundary Pressure values
    }

    if ( fin.fail() )
    {
        PRINT_ERROR("not enough data to read\n");
        exit(1);
    }
    fin.close();

    //////// initialize
    omega_   = freq_ * 2 * M_PI;
    waveNum_ = omega_ / speed_;
    iOmegaRho_ = TComplex(0, omega_*density_);

    printf("============ Equivalent Source Simulation ============\n");
    printf("Number of elements      = %d\n", (int)tgl_.size());
    printf("Number of nodes         = %d\n", (int)vtx_.size());
    printf("Number of field points  = %d\n", (int)fieldPts_.size());
    printf("Frequency               = %lf\n", freq_);
    printf("Wave number (K)         = %lf\n", waveNum_);
    printf("Speed of sound (c)      = %lf\n", speed_);
    printf("Mass density (rho)      = %lf\n", density_);
}

#endif

/*
template <typename T> 
BoundaryIntegralTransfer<T>::TComplex 
BoundaryIntegralTransfer<T>::eval(const Point3<T>& pt) const
{
    const T ONE_THIRD = 1. / 3.;
    TComplex ret = 0;

    for(int i = 0;i < tgl_.size();++ i)
    {
        Point3<T> c = (vtx_[tgl_[i][0]] + vtx_[tgl_[i][1]] + vtx_[tgl_[i][2]]) * ONE_THIRD;
        TComplex val, deri;
        helmholtz_green_func_with_dir_deri(waveNum_, pt, c, tglNml_[i], val, deri);
        ret += iOmegaRho_ * v_[i] * val * tglArea_[i];
        ret -= p_[i] * deri * tglArea_[i];
    }
    return ret;
}
*/

