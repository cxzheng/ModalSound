/******************************************************************************
 *  File: RigidBody.hpp
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
#ifndef RIGID_BODY_HPP
#   define RIGID_BODY_HPP

#include <stdlib.h>
#include "geometry/Point3.hpp"
#include "linearalgebra/Vector3.hpp"
#include "linearalgebra/Quaternion.hpp"

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

/*!
 * The basic rigid body
 */
template <typename T, class TMesh>
class RigidBody
{
    public:
        const TMesh* mesh() const
        {  return mp_mesh; }

        RigidBody(TMesh* mesh, T density):m_fixed(false), 
                m_density(density), mp_mesh(mesh)
        {
            init();
        }

        /*
         * given the initial position of a point on this rigid body,
         * return the current position of that point due to the movement 
         * of this rigid body
         */
        inline Point3<T> current_position(const Point3<T>& pt) const
        {   return m_x + m_q.rotate(pt - m_x0);  }
        
        void advance_velocity(T dt)
        {
            m_v.scaleAdd(dt, m_acc);
            m_omega.scaleAdd(dt, m_angAcc);     // advance angular velocity
        }

        void advance_state(T dt)
        {
            m_x.scaleAdd(dt, m_v);

            // the rotation in dt is [cos(|omega|*dt*0.5),
            // sin(|omega|*dt*0.5)*omega/|omega|)]
            T omeganorm = m_omega.length();
            if ( omeganorm > EPS )
            {
                const T dang = omeganorm * dt * 0.5;
                Quaternion<T> rot(cos(dang), m_omega * (sin(dang) / omeganorm));
                m_q = rot * m_q;
                m_q.normalize();
            }
            m_invq = m_q.conjugate();

            // m_invI = R * m_invI0 * R^T
            const Matrix3<T> R = m_q.toMatrix3();
            m_invI = R * m_invI0 * R.transpose();
        }

        /*
         * Return the velocity of a given point on this rigid body
        inline Vector3<T> current_velocity(const Point3<T>& pt) const
        {
            return m_v + m_omega.crossProduct(pt - m_x);
        }
         */

        inline T mass_inverse() const
        {   return m_invMass; }

        const Point3<T>& initial_mass_center() const
        {  return m_x0; }
        const Point3<T>& mass_center() const
        {  return m_x; }
        const Quaternion<T>& rotation() const
        {  return m_q; }
        T density() const
        {  return m_density; }
        const Vector3<T>& angular_velocity() const
        {  return m_omega; }
        bool is_fixed() const
        {  return m_fixed; }

        void set_fixed(bool f)
        {  m_fixed = f; }

    protected:
        void init();

    protected:
        bool            m_fixed;
        T               m_density;

        /*======= constant quantities =======*/
        T               m_mass;
        T               m_invMass;
        Point3<T>       m_x0;       // center of mass
        Matrix3<T>      m_invI0;    // inverse of inertia at initial configuration

        /*======= State Variables ======*/
        Point3<T>       m_x;        // current position
        Quaternion<T>   m_q;        // rotation
        Quaternion<T>   m_invq;     // inverse rotation: rotate from current state to initial state
        Matrix3<T>      m_invI;     // inverse of current inertia

        /*======= Derived quantities =======*/
        Vector3<T>      m_v;        // current velocity
        Vector3<T>      m_omega;    // current angular velocity

        /*======= Computed quantities =======*/
        Vector3<T>      m_acc;      // acceleration
        Vector3<T>      m_angAcc;   // angular acceleration

        TMesh*          mp_mesh;
};

//////////////////////////////////////////////////////////////////////////////////////

template <typename T, class TMesh>
void RigidBody<T, TMesh>::init()
{
    const std::vector<T>&           ms  = mp_mesh->masses();
    const std::vector< Point3<T> >& vtx = mp_mesh->vertices();

    Matrix3<T> I0;  // inertia
    m_mass = 0;
    m_x0.zero();
    for(size_t i = 0;i < ms.size();++ i)
    {
        m_mass += ms[i];
        m_x0 += ms[i] * vtx[i];
    }

    if ( m_mass < EPS )
    {
        fprintf(stderr, "ERROR: object has zero mass\n");
        exit(1);
    }

    m_x0 /= m_mass;
    m_mass *= m_density;
    m_invMass = 1. / m_mass;

    m_x = m_x0;

    // update inertia
    for(size_t i = 0;i < ms.size();++ i)
    {
        Vector3<T> ri = vtx[i] - m_x0;
        I0 += Matrix3<T>(
                ms[i]*(M_SQR(ri.y)+M_SQR(ri.z)), -ms[i]*ri.x*ri.y, -ms[i]*ri.x*ri.z,
                -ms[i]*ri.y*ri.x, ms[i]*(M_SQR(ri.x)+M_SQR(ri.z)), -ms[i]*ri.y*ri.z,
                -ms[i]*ri.z*ri.x, -ms[i]*ri.z*ri.y, ms[i]*(M_SQR(ri.x)+M_SQR(ri.y)));
    }
    I0 *= m_density;
    m_invI0 = I0.inverse();
    m_invI  = m_invI0;
}

#ifdef USE_NAMESPACE
}
#endif
#endif
