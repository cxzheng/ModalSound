#ifndef FMM_TRANSFER_EVAL_INC
#   define FMM_TRANSFER_EVAL_INC

#include <complex>
#include <vector>
#include "sc/Matrix.hpp"
#include "sc/SphericalFunc.hpp"
#include "HelmholtzBasis.hpp"
#include "sploosh.pb.h"
#include "config.h"

class FMMTransferEval
{
    public:
        typedef std::complex<REAL>  TComplex;

        FMMTransferEval(const char* dir);
        FMMTransferEval(const sploosh::FMMoments& ms, const std::string& dir);

        ~FMMTransferEval()
        {   delete moments_; }

        /*
         * Evaluate the sound pressure at a given position
         */
        TComplex eval(const Point3<REAL>& pt) const
        {
            TComplex ret = 0;
        
            Point3<REAL> sc;
            cartesian_to_spherical(center_, pt, sc);
        
            for(int n = 0, cptr = 0;n < nexpan_;++ n)
            for(int m = -n;m <= n;++ m)
            {
                TComplex c((*moments_)[cptr][0], (*moments_)[cptr][1]);
                ++ cptr;
                ret += c * HelmholtzBasis::singular_basis(m, n, waveNum_, sc);
            }
        
            return TComplex(-waveNum_*ret.imag(), waveNum_*ret.real()); // ret * ik
        }

        inline int expansion_num() const
        {   return nexpan_; }

        inline REAL wave_number() const
        {   return waveNum_; }

    private:
        int                         nexpan_;        // expansion level
        REAL                        waveNum_;
        Point3<REAL>                center_;
        Matrix<REAL>*               moments_;       // precomputed multipole moment 
                                                    // nexp^2 x 2 matrix
};

#endif
