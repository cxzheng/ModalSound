#ifndef FMM_BOUNDARY_INTEGRAL_TRANSFER_HPP
#   define FMM_BOUNDARY_INTEGRAL_TRANSFER_HPP

#include <string.h>
#include <fstream>
#include <limits>

#include "utils/term_msg.h"
#include "HelmholtzBasis.hpp"
#include "BoundaryIntegralTransfer.hpp"
#include "io/EigenMatrixIO.hpp"

/*
 * Evaluate the transfer function value (P) using the boundary element 
 * integration with Fast Multipole Method
 * 
 * The input files are the input and output files for FastBEM solver. Then 
 * it evaluates the transfer value (P) at given location by computing the boundary 
 * integral.
 */
template <typename T>
class FMMBoundIntTransfer : public BoundaryIntegralTransfer<T>
{
    public:
        typedef BoundaryIntegralTransfer<T>     TParent;
        typedef typename TParent::TComplex      TComplex;

        FMMBoundIntTransfer(const char* fileinput, const char* fileoutput);
        /*
         * Initialize the integrator with given mass center
         */
        FMMBoundIntTransfer(const char* fileinput, const char* fileoutput, 
                            const Point3<T>& mc);

        /*
         * precompute moments
         */
        void compute_moments();
        /*
         * write the precomputed moments into file
         * make sure moments are computed before calling this method
         */
        void store_moments(const char* file);

        /*
         * Evaluate the sound pressure at a given position
         */
        TComplex eval(const Point3<T>& pt) const
        {
            TComplex ret = 0;
        
            Point3<T> sc;
            cartesian_to_spherical(objCenter_, pt, sc);
        
            for(int n = 0, cptr = 0;n < nexpan_;++ n)
            for(int m = -n;m <= n;++ m)
                ret += preMoments_[cptr ++] * HelmholtzBasis::singular_basis(m, n, waveNum_, sc);
        
            return TComplex(-waveNum_*ret.imag(), waveNum_*ret.real()); // ret * ik
        }

        inline const Point3<T>& obj_center() const
        {   return objCenter_;  }
        inline int expansion_num() const
        {   return nexpan_;     }
        inline const std::vector<TComplex>& moments() const
        {   return preMoments_; }
        //void set_expansion_num(int nexpan) 
        //{   nexpan_ = nexpan; }

    private:
        void init();

    protected:
        using   TParent::v_;
        using   TParent::p_;
        using   TParent::vtx_;
        using   TParent::tgl_;
        using   TParent::tglArea_;
        using   TParent::tglNml_;
        using   TParent::iOmegaRho_;
        using   TParent::waveNum_;
        using   TParent::speed_;
        using   TParent::density_;
        using   TParent::omega_;

        int                         nexpan_;         // expansion level
        Point3<T>                   objCenter_;
        std::vector< TComplex >     preMoments_;     // precomputed multipole moment 
};

///////////////////////////////////////////////////////////////////////////////

template <typename T>
FMMBoundIntTransfer<T>::FMMBoundIntTransfer(
        const char* fileinput, 
        const char* fileoutput):
        BoundaryIntegralTransfer<T>(fileinput, fileoutput)
{
    //// find the center of object, it should not be exactly the same with any triangle center
    Point3<T> minPt(std::numeric_limits<T>::infinity(), 
                    std::numeric_limits<T>::infinity(),
                    std::numeric_limits<T>::infinity()); 
    Point3<T> maxPt(-std::numeric_limits<T>::infinity(), 
                    -std::numeric_limits<T>::infinity(),
                    -std::numeric_limits<T>::infinity()); 
    // finding a bounding box
    for(size_t i = 0;i < vtx_.size();++ i)
    {
        minPt.x = fmin(vtx_[i].x, minPt.x);
        minPt.y = fmin(vtx_[i].y, minPt.y);
        minPt.z = fmin(vtx_[i].z, minPt.z);

        maxPt.x = fmax(vtx_[i].x, maxPt.x);
        maxPt.y = fmax(vtx_[i].y, maxPt.y);
        maxPt.z = fmax(vtx_[i].z, maxPt.z);
    }

    objCenter_ = (maxPt + minPt) * 0.5;
    
    init();
}

template <typename T>
FMMBoundIntTransfer<T>::FMMBoundIntTransfer(const char* fileinput, 
        const char* fileoutput, const Point3<T>& mc):
        objCenter_(mc),
        BoundaryIntegralTransfer<T>(fileinput, fileoutput)
{
    init();
}

template <typename T>
void FMMBoundIntTransfer<T>::init()
{
    //// NOTE: To make it more stable, make sure this objCenter_ is not overlap 
    //   with any triangle center
    for(size_t i = 0;i < vtx_.size();++ i)
        if ( vtx_[i].distance_sqr(objCenter_) < PrecisionType<T>::EPS )
            PRINT_WARNING("FMMBoundIntTransfer >> The given object center is too close to a vertex\n");
    
    //// estimate number of expansion
    if ( omega_ < 125.6637 || omega_ > 125663.7 )
    {
        PRINT_WARNING("The given frequency(%lf) is out of range. Ignore it\n", omega_);
        nexpan_ = 0;
    }
    else
    {
        const int MAX_NEXPAN = 40;
        const int MIN_NEXPAN = 16;
        nexpan_ = (int)((omega_ - 125.6637) * (MAX_NEXPAN - MIN_NEXPAN + 1) / 
                        (125663.7 - 125.6637) + MIN_NEXPAN);
        nexpan_ = std::min(MAX_NEXPAN, std::max(nexpan_, MIN_NEXPAN));
    }

    printf("# of FMM expansions     = %d\n", nexpan_);
    printf("======================================================\n");
} // end init()

template <typename T>
void FMMBoundIntTransfer<T>::compute_moments()
{
    const T ONE_THIRD = 1. / 3.;

    preMoments_.resize(nexpan_*nexpan_);
    memset(&preMoments_[0], 0, sizeof(TComplex)*preMoments_.size());

    for (int i = 0;i < tgl_.size();++ i)
    {
        Point3<T> sc;
        Point3<T> tc = (vtx_[tgl_[i][0]] + vtx_[tgl_[i][1]] + vtx_[tgl_[i][2]]) * ONE_THIRD;
        cartesian_to_spherical(objCenter_, tc, sc);

        for(int n = 0, cptr = 0;n < nexpan_;cptr += ((++ n)<<1))
        {
            //// compute the R^0_n
            TComplex A = HelmholtzBasis::regular_basis(0, n, waveNum_, sc); 
            TComplex B = HelmholtzBasis::regular_basis_dir_deriv(0, n, waveNum_, sc, tglNml_[i]);
            TComplex t = iOmegaRho_ * v_[i];
            preMoments_[cptr] += tglArea_[i] * (t*A - p_[i]*B);

            //// compute the R^m_n utilizing R^{-m}_n = conjugate(R^m_n)
            for(int m = 1;m <= n;++ m)
            {
                A = HelmholtzBasis::regular_basis(m, n, waveNum_, sc); 
                B = HelmholtzBasis::regular_basis_dir_deriv(m, n, waveNum_, sc, tglNml_[i]);
                preMoments_[cptr - m] += tglArea_[i] * (t*A - p_[i]*B);
                preMoments_[cptr + m] += tglArea_[i] * (t*conj(A) - p_[i]*conj(B));
            }
        }
    }
}

template <typename T>
void FMMBoundIntTransfer<T>::store_moments(const char* file)
{
    Eigen::Matrix<T, Eigen::Dynamic, 2> mat(nexpan_*nexpan_, 2);

    for(int n = 0, cptr = 0;n < nexpan_;++ n)
    for(int m = -n;m <= n;++ m)
    {
        mat(cptr, 0) = preMoments_[cptr].real();
        mat(cptr, 1) = preMoments_[cptr].imag();
        ++ cptr;
    }

    if ( _FAILED(write_ma_eigen_matrix(file, &mat)) )
    {
        PRINT_ERROR("Cannot write to file: %s\n", file);
    }

#if 0
    std::ofstream fout(file, std::ios::binary);
    if ( fout.fail() )
    {
        PRINT_ERROR("Cannot open file: %s\n", file);
        return;
    }

    fout.write((const char*)&nexpan_, sizeof(int));
    for(int n = 0, cptr = 0;n < nexpan_;++ n)
    for(int m = -n;m <= n;++ m)
        fout.write((const char*)&preMoments_[cptr ++], sizeof(TComplex));

    fout.flush();
    fout.close();
#endif
}

#endif

