#ifndef DEFORMABLE_LINEAR_H
#   define DEFORMABLE_LINEAR_H

#include "config.h"
#include "geometry/Tet.hpp"
#include "linearalgebra/Matrix.hpp"

class LinearMaterial
{
    public:
        LinearMaterial() { }
        LinearMaterial(REAL lam, REAL mu):lambda_(lam), mu_(mu)
        { }

        void set_parameters(REAL l, REAL m)
        {
            lambda_ = l;
            mu_ = m;
        }

        /*! 
         * evaluate the stiffness matrix for a single tetrahedron 
         * AT REST POSITION
         * Use the formula given in paper Brien et al. 2002 (Synthesizing Sounds from Rigid-Body Simulations)
         */
        void rest_stiffness_matrix(const Tet<REAL>& tet,
                Matrix<REAL>& out) const;
    private:
        REAL    lambda_;
        REAL    mu_;
};

#endif
