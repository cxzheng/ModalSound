#ifndef HELMHOLTZ_SOLUTION_HPP
#   define HELMHOLTZ_SOLUTION_HPP

#include "sc/SpecialFunc.hpp"
#include "sc/SphericalFunc.hpp"

/*
 * Evaluate the Helmholtz Green's function
 *
 * \param k  wave number
 * \param r  distance to the center
 */
template <typename T>
static inline std::complex<T> helmholtz_green_func(T k, T r)
{
    const T GS = 0.25 * M_1_PI / r;  // 1/(4*pi*r)
    return std::complex<T>(cos(k*r)*GS, sin(k*r)*GS);
}

/*
 * Evaluate the directional derivative of Helmholtz Green's function
 */
template <typename T>
static inline std::complex<T> helmholtz_green_func_dir_deri(T wn, 
        const Point3<T>& src, const Point3<T>& pos, const Vector3<T>& dir)
{
    const Vector3<T> rv = pos - src;
    const T r = rv.length();
    const T invr = 1. / r;
    const T kr = wn * r;

    const T prpx = rv.x * invr, prpy = rv.y * invr, prpz = rv.z * invr;
    const T sfr  = prpx*dir.x + prpy*dir.y + prpz*dir.z;

    return sfr * 0.25 * M_1_PI * std::complex<T>(-invr*invr, wn*invr) * 
            std::complex<T>(cos(kr), sin(kr));
}

/*
 * Evaluate both the Helmholtz Green's function value and its directional
 * derivative
 */
template <typename T>
static inline void helmholtz_green_func_with_dir_deri(T wn, 
        const Point3<T>& src, const Point3<T>& pos, const Vector3<T>& dir,
        std::complex<T>& val, std::complex<T>& deri)
{
    const Vector3<T> rv = pos - src;
    const T r = rv.length();
    const T invr = 1. / r;
    const T kr = wn * r;
    const T GS = 0.25 * M_1_PI / r;  // 1/(4*pi*r)

    val = std::complex<T>(cos(kr)*GS, sin(kr)*GS);

    const T prpx = rv.x * invr, prpy = rv.y * invr, prpz = rv.z * invr;
    const T sfr  = prpx*dir.x + prpy*dir.y + prpz*dir.z;

    deri = sfr * 0.25 * M_1_PI * std::complex<T>(-invr*invr, wn*invr) * 
            std::complex<T>(cos(kr), sin(kr));
}

/*
 * compute j_n * sph(m, n), the basic solution of Helmholtz equation
 */
template <typename T, int Order>
class SphericalBesselHarmonics
{
    public:
        /*
         * \param wn    the wave number
         * \param src   the source position
         * \param pos   the receiver position
         * \param output store the output data
         * \param incl  increment of the output pointer
         */
        void eval(T wn, const Point3<T>& src, const Point3<T>& pos, 
                  std::complex<T>* output, int incl)
        {
            Point3<T> scoord;

            cartesian_to_spherical(src, pos, scoord);
            const T kr = wn * scoord.r;

            // compute j0(kr) * sph(0, 0)
            output[0] = spb_.j0(kr) * sph_.Y0();
            int idx = incl;

            if ( Order > 0 )
            {
                const T j1 = spb_.j1(kr);
                for(int i = -1;i <= 1;++ i, idx += incl)
                    output[idx] = j1 * sph_.Y1(i, scoord.theta, scoord.phi);
            }

            if ( Order > 1 )
            {
                const T j2 = spb_.j2(kr);
                for(int i = -2;i <= 2;++ i, idx += incl)
                    output[idx] = j2 * sph_.Y2(i, scoord.theta, scoord.phi);
            }
        }

        /* 
         * compute directional derivative 
         */
        void dir_deriv(T wn, const Point3<T>& src, const Point3<T>& pos,
                const Vector3<T>& dir, std::complex<T>* output, int incl)
        {
            Point3<T> scoord;
            const Vector3<T> rv = pos - src;

            const T xysqr = cartesian_to_spherical(rv, scoord);
            const T xyd   = sqrt(xysqr);
            const T invxysqr = 1. / xysqr;
            const T invxyd   = 1. / xyd;
            const T kr = wn * scoord.r;

            const T invr = 1. / scoord.r;
            const T invr2 = M_SQR(invr);

            const T prpx = rv.x * invr, prpy = rv.y * invr, prpz = rv.z * invr;
            const T pthetapx = rv.x*rv.z*invr2*invxyd;
            const T pthetapy = rv.y*rv.z*invr2*invxyd;
            const T pthetapz = -xyd*invr2;
            const T pphipx   = -rv.y * invxysqr;
            const T pphipy   =  rv.x * invxysqr;
            // const T pphipz   = 0;

            const T sfr     = prpx*dir.x + prpy*dir.y + prpz*dir.z;
            const T sftheta = pthetapx*dir.x + pthetapy*dir.y + pthetapz*dir.z; 
            const T sfphi   = pphipx*dir.x + pphipy*dir.y; // + pphipz*dir.z;

            // order 0:
            output[0] = sfr*wn*spb_.j0_deriv(kr)*sph_.Y0();
            int idx = incl;

            if ( Order > 0 )
            {
                const T j1      = spb_.j1(kr);
                const T j1deriv = wn * spb_.j1_deriv(kr);

                for(int i = -1;i <= 1;++ i, idx += incl)
                    output[idx] =  sfr * j1deriv * sph_.Y1(i, scoord.theta, scoord.phi) +
                                   sftheta * j1 * sph_.Y1_deriv_theta(i, scoord.theta, scoord.phi) +
                                   sfphi * j1 * sph_.Y1_deriv_phi(i, scoord.theta, scoord.phi);
            }

            if ( Order > 1 )
            {
                const T j2      = spb_.j2(kr);
                const T j2deriv = wn * spb_.j2_deriv(kr);
                for(int i = -2;i <= 2;++ i, idx += incl)
                    output[idx] = sfr * j2deriv * sph_.Y2(i, scoord.theta, scoord.phi) + 
                                  sftheta * j2 * sph_.Y2_deriv_theta(i, scoord.theta, scoord.phi) + 
                                  sfphi * j2 * sph_.Y2_deriv_phi(i, scoord.theta, scoord.phi);
            }
        }

    private:
        SphericalBessel<T>     spb_;
        SphericalHarmonics<T>  sph_;
};

template <typename T, int Order>
class SphericalHankelHarmonics
{
    public:
        /*
         * \param wn    the wave number
         * \param src   the source position
         * \param pos   the receiver position
         * \param output store the output data
         * \param incl  increment of the output pointer
         */
        void eval_2nd(T wn, const Point3<T>& src, const Point3<T>& pos, 
                     std::complex<T>* output, int incl) const
        {
            Point3<T> scoord;
            cartesian_to_spherical(src, pos, scoord);
            const T kr = wn * scoord.r;

            // compute h0(kr)*sph(0,0)
            output[0] = spk_.h2nd0(kr) * sph_.Y0();
            int idx = incl;

            if ( Order > 0 )
            {
                const std::complex<T> h1 = spk_.h2nd1(kr);
                for(int i = -1;i <=1;++ i, idx += incl)
                    output[idx] = h1 * sph_.Y1(i, scoord.theta, scoord.phi);
            }

            if ( Order > 1 )
            {
                const std::complex<T> h2 = spk_.h2nd2(kr);
                for(int i = -2;i <= 2;++ i, idx += incl)
                    output[idx] = h2 * sph_.Y2(i, scoord.theta, scoord.phi);
            }
        }

        /*
         * m = 0 -->    l = 0, m = 0
         * m =[1,3] --> l = 1, m = m-2
         * m =[4,8] --> l = 2, m = m-6
         */
        std::complex<T> eval_2nd(T wn, const Point3<T>& sc, int m) const
        {
            const T kr = wn * sc.r;

            // compute h0(kr)*sph(0,0)
            if ( m == 0 ) return spk_.h2nd0(kr) * sph_.Y0();
            if ( m < 4 ) 
            {
                assert(Order > 0);
                const std::complex<T> h1 = spk_.h2nd1(kr);
                return h1 * sph_.Y1(m-2, sc.theta, sc.phi);
            }
            assert(Order > 1 && m < 9);
            const std::complex<T> h2 = spk_.h2nd2(kr);
            return h2 * sph_.Y2(m-6, sc.theta, sc.phi);
        }

        /*
         * compute the normal derivative of h0*sph
         * See the gradient of spherical coordinate at
         * http://mathworld.wolfram.com/SphericalCoordinates.html
         *
         * \param wn  wave number K
         * \param src source position
         * \param pos surface sampling point position
         * \param dir normal direction
         */
        void dir_deriv_2nd(T wn, const Point3<T>& src, const Point3<T>& pos,
                const Vector3<T>& dir, std::complex<T>* output, int incl) const
        {
            Point3<T> scoord;
            const Vector3<T> rv = pos - src;

            const T xysqr = cartesian_to_spherical(rv, scoord);
            const T xyd   = sqrt(xysqr);
            const T invxysqr = 1. / xysqr;
            const T invxyd   = 1. / xyd;
            const T kr = wn * scoord.r;

            const T invr = 1. / scoord.r;
            const T invr2 = M_SQR(invr);

            const T prpx = rv.x * invr, prpy = rv.y * invr, prpz = rv.z * invr;
            const T pthetapx = rv.x*rv.z*invr2*invxyd;
            const T pthetapy = rv.y*rv.z*invr2*invxyd;
            const T pthetapz = -xyd*invr2;
            const T pphipx   = -rv.y * invxysqr;
            const T pphipy   =  rv.x * invxysqr;
            // const T pphipz   = 0;

            const T sfr     = prpx*dir.x + prpy*dir.y + prpz*dir.z;
            const T sftheta = pthetapx*dir.x + pthetapy*dir.y + pthetapz*dir.z; 
            const T sfphi   = pphipx*dir.x + pphipy*dir.y; // + pphipz*dir.z;

            // order 0:
            output[0] = sfr*wn*spk_.h2nd0_deriv(kr) * sph_.Y0();
            int idx = incl;

            if ( Order > 0 )
            {
                const std::complex<T> h1      = spk_.h2nd1(kr);
                const std::complex<T> h1deriv = wn * spk_.h2nd1_deriv(kr);

                for(int i = -1;i <= 1;++ i, idx += incl)
                    output[idx] = sfr     * h1deriv * sph_.Y1(i, scoord.theta, scoord.phi) +
                                  sftheta * h1      * sph_.Y1_deriv_theta(i, scoord.theta, scoord.phi) +
                                  sfphi   * h1      * sph_.Y1_deriv_phi(i, scoord.theta, scoord.phi);
            }

            if ( Order > 1 )
            {
                const std::complex<T> h2      = spk_.h2nd2(kr);
                const std::complex<T> h2deriv = wn * spk_.h2nd2_deriv(kr);

                for(int i = -2;i <= 2;++ i, idx += incl)
                    output[idx] = sfr     * h2deriv * sph_.Y2(i, scoord.theta, scoord.phi) + 
                                  sftheta * h2      * sph_.Y2_deriv_theta(i, scoord.theta, scoord.phi) +
                                  sfphi   * h2      * sph_.Y2_deriv_phi(i, scoord.theta, scoord.phi);
            }
        }

        std::complex<T> dir_deriv_2nd(T wn, const Point3<T>& src, 
                const Point3<T>& pos, const Vector3<T>& dir, int m) const
        {
            Point3<T> scoord;
            const Vector3<T> rv = pos - src;

            const T xysqr = cartesian_to_spherical(rv, scoord);
            const T xyd   = sqrt(xysqr);
            const T invxysqr = 1. / xysqr;
            const T invxyd   = 1. / xyd;
            const T kr = wn * scoord.r;

            const T invr = 1. / scoord.r;
            const T invr2 = M_SQR(invr);

            const T prpx = rv.x * invr, prpy = rv.y * invr, prpz = rv.z * invr;
            const T sfr  = prpx*dir.x + prpy*dir.y + prpz*dir.z;

            if ( m == 0 ) return sfr * wn * spk_.h2nd0_deriv(kr) * sph_.Y0();

            const T pthetapx = rv.x*rv.z*invr2*invxyd;
            const T pthetapy = rv.y*rv.z*invr2*invxyd;
            const T pthetapz = -xyd*invr2;
            const T pphipx   = -rv.y * invxysqr;
            const T pphipy   =  rv.x * invxysqr;
            // const T pphipz   = 0;
            
            const T sftheta = pthetapx*dir.x + pthetapy*dir.y + pthetapz*dir.z; 
            const T sfphi   = pphipx*dir.x + pphipy*dir.y; // + pphipz*dir.z;

            if ( m < 4 )
            {
                assert(Order > 0);
                const std::complex<T> h1      = spk_.h2nd1(kr);
                const std::complex<T> h1deriv = wn * spk_.h2nd1_deriv(kr);
                return sfr     * h1deriv * sph_.Y1(m-2, scoord.theta, scoord.phi) +
                       sftheta * h1      * sph_.Y1_deriv_theta(m-2, scoord.theta, scoord.phi) +
                       sfphi   * h1      * sph_.Y1_deriv_phi(m-2, scoord.theta, scoord.phi);
            }

            assert(Order > 1 && m < 9);
            const std::complex<T> h2      = spk_.h2nd2(kr);
            const std::complex<T> h2deriv = wn * spk_.h2nd2_deriv(kr);

            return sfr     * h2deriv * sph_.Y2(m-6, scoord.theta, scoord.phi) + 
                   sftheta * h2      * sph_.Y2_deriv_theta(m-6, scoord.theta, scoord.phi) +
                   sfphi   * h2      * sph_.Y2_deriv_phi(m-6, scoord.theta, scoord.phi);
        }


    private:
        SphericalHankel<T>      spk_;
        SphericalHarmonics<T>   sph_;
};
#endif
