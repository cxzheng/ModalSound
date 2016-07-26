/******************************************************************************
 *  File: SimpleDistFunc.h
 *
 *  This file is part of isostuffer
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
#ifndef SIMPLE_DIST_FUNC
#   define SIMPLE_DIST_FUNC

#include "linearalgebra/Vector3.hpp"
#include "geometry/Point3.hpp"
#include "utils/math.hpp"

class CylinderDistFunc
{
    public:
        CylinderDistFunc(double rad, double h):bc_(0.,0.,0.), 
                height_(h), rad_(rad), orient_(0, 1, 0)
        {}
        CylinderDistFunc(const Point3d& c, double rad, double h):
                bc_(c), height_(h), rad_(rad), orient_(0, 1, 0)
        {}
        CylinderDistFunc(const Point3d& c, const Vector3d& ori, double rad, double h):
                bc_(c), height_(h), rad_(rad), orient_(ori)
        {}

        inline double operator() (const Point3d& pt) const
        {
            Vector3d dvec = pt - bc_;
            double dh = dvec.dotProduct(orient_);
            double dr = sqrt(dvec.lengthSqr() - dh*dh);
            if ( dh < 0 || dh > height_ )
            {
                if ( dr < rad_ ) 
                    return dh < 0 ? -dh : dh - height_;
                else
                    return dh < 0 ? sqrt(dh*dh + M_SQR(dr-rad_)) :
                                    sqrt(M_SQR(dh-height_) + M_SQR(dr-rad_));
            }
            else
                return dr - rad_;
        }

    private:
        Point3d     bc_;            // bottom center
        double      height_;
        double      rad_;
        Vector3d    orient_;
};

#endif
