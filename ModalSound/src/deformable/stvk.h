/******************************************************************************
 *  File: stvk.hpp
 *  Implement the St. Venant-Kirchhoff model
 *  Copyright (c) 2007 by Changxi Zheng
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
 ******************************************************************************/
#ifndef DEFORMABLE_STVK_H
#   define DEFORMABLE_STVK_H

#include "config.h"
#include "geometry/Tet.hpp"
#include "sc/Matrix3.hpp"
#include "sc/Matrix.hpp"

class StVKMaterial
{
    public:
        StVKMaterial() { }
        StVKMaterial(REAL lam, REAL mu):lambda_(lam), mu_(mu)
        { }
        StVKMaterial(const StVKMaterial& mat):
                lambda_(mat.lambda_), mu_(mat.mu_)
        { }
        StVKMaterial(const StVKMaterial* mat):
                lambda_(mat->lambda_), mu_(mat->mu_)
        { }

        void set_parameters(REAL l, REAL m)
        {
            lambda_ = l;
            mu_ = m;
        }

        /*! first Piola Kirchhoff stress tensor */
        void first_piola_kirchhoff(const Matrix3<REAL>& F, 
                Matrix3<REAL>& out) const;
        /*! second Piola Kirchhoff stress tensor */
        void second_piola_kirchhoff(const Matrix3<REAL>& F, 
                Matrix3<REAL>& out) const;      // $$TESTED
        void cauchy_stress_tensor(const Matrix3<REAL>& F,
                Matrix3<REAL>& out) const;
        /*! evaluate the stiffness matrix for a single tetrahedron */
        void stiffness_matrix(const Tet<REAL>& tet,
                Matrix<REAL>& out) const;       // $$TESTED
        /*! strain energy density in undeformed volumn */
        REAL strain_energy_density(const Matrix3<REAL>& F) const;

    protected:
        REAL    lambda_;
        REAL    mu_;
};

#endif
