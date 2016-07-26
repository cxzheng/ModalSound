/******************************************************************************
 *  File: LSCollisionRigidBody.hpp
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
#ifndef LSCOLLISION_RIGID_BODY_HPP
#   define LSCOLLISION_RIGID_BODY_HPP

#include <vector>
#include <map>

#include "RigidBody.hpp"
#include "LSCollisionDetect.hpp"
#include "linearalgebra/Matrix3.hpp"
#include "linearalgebra/Vector3.hpp"

#ifdef USE_RECORDER
#   include "geometry/KDTree.hpp"
#   include "geometry/Point3KDStruct.hpp"
#endif

#ifdef USE_NAMESPACE
namespace carbine
{
#endif

class RigidBodySimulator;
class ContactGraph;

/*!
 * rigid body with Level-set based collision detection
 */
template <typename T, class _TMesh>
class LSCollisionRigidBody : public RigidBody<T, _TMesh>
{
    friend class RigidBodySimulator;
    friend class ContactGraph;

    public:
        typedef LSCollisionRigidBody<T, _TMesh>         TSelf;
        typedef CollisionProcessor<T, TSelf>            TCollProc;
        typedef _TMesh                                  TMesh;

    public:
        LSCollisionRigidBody(int id, TMesh* mesh, T density, T restCoeff, T mu):
                RigidBody<T, TMesh>(mesh, density), m_id(id), m_predTs(0),
                m_restCoeff(restCoeff), m_mu(mu), m_cproc(this, mesh), 
                m_pinned(false)
        { 
#ifdef USE_RECORDER
            init_kdtree();
#endif
        }

        /* 
         * apply external forces to the rigid body 
         */
        void apply_forces()
        {
            if ( m_fixed )
            {
                m_acc.zero();
                m_angAcc.zero();
                return;
            }

            m_acc.set((T)0, -(T)GRAVITY, (T)0);
            m_acc += m_extF * m_invMass;

            //m_acc.set((T)GRAVITY*sin(M_DEG2RAD(20.)), -(T)GRAVITY*cos(M_DEG2RAD(20.)), 0.);
            m_angAcc.zero();
            m_angAcc += m_invI*m_extT;

            clean_external_force();
        }

        void clean_external_force()
        {
            m_extF.zero();
            m_extT.zero();
        }

        /* input the force from external */
        void apply_external_force(const Vector3<T>& f, const Point3<T>& p)
        {
            m_extF += f;
            m_extT += (p - m_x).crossProduct(f);
        }

        ///* the vertex position in predicted configuration */
        //inline Point3<T> predicted_position(int vtxid) const
        //{  return m_predx + m_predq.rotate(mp_mesh->vertex(vtxid) - m_x0); }

        inline Point3<T> predicted_position(const Point3<T>& pt) const
        {  return m_predx + m_predq.rotate(pt - m_x0); }

        /*
         * given the predicted position of a point on this rigid body,
         * return the initial position of that point
         */
        Point3<T> initial_predicted_position(const Point3<T>& pt) const
        {  return m_x0 + m_predinvq.rotate(pt - m_predx); }

        /*
         * given a point in predicted configuration, return the position
         * of that point in current configuration
        inline Point3<T> pred2cur_position(const Point3<T>& pt) const
        {
            Vector3<T> r = pt - m_predx;    // r in pred config
            m_predinvq.rotate_vector(r);    // r in initial config
            m_q.rotate_vector(r);           // r in current config
            return m_x + r;
        }
         */

        inline Vector3<T> predicted_velocity(const Point3<T>& pt) const
        { return m_predv + m_predw.crossProduct(pt - m_predx); }

        /*
         * NOTE: make sure nml has been normalized
         */
        inline Vector3<T> predicted_normal(const Vector3<T>& nml) const
        {  return m_predq.rotate(nml); }

        /*
         * apply impulse on this object, changing the current velocity
         */
        void apply_impulse_to_prediction(
                const Vector3<T>& impl, const Vector3<T>& r);

        /*
         * update predicted states by using current force
         * predv = v + dt*a
         * predx = x + dt(v + dt*a)
         */
        void update_force_predicted_position(REAL dt);

        /*
         * update predicted states by using current velocity
         * predx = x + dt*v
         */
        void update_velocity_predicted_position(REAL dt);

        const Quaternion<T>& predicted_inverse_rotation() const
        {  return m_predinvq; }
        TCollProc* collision_processor() 
        {  return &m_cproc; }
        const TCollProc* collision_processor() const
        {  return &m_cproc; }
        inline T friction_coeff() const
        {  return m_mu; }
        inline T rest_coeff() const
        {  return m_restCoeff; } 
        int id() const
        {  return m_id; }

    protected:
        using RigidBody<T, TMesh>::m_fixed;
        using RigidBody<T, TMesh>::m_x0;
        using RigidBody<T, TMesh>::m_v;
        using RigidBody<T, TMesh>::m_acc;
        using RigidBody<T, TMesh>::m_omega;
        using RigidBody<T, TMesh>::m_x;
        using RigidBody<T, TMesh>::m_q;
        using RigidBody<T, TMesh>::m_angAcc;
        using RigidBody<T, TMesh>::m_invI0;
        using RigidBody<T, TMesh>::m_invI;
        using RigidBody<T, TMesh>::mp_mesh;
        using RigidBody<T, TMesh>::m_invMass;

        int                 m_id;           // obj ID
        /*======= Predicted State Variables ======*/
        int                 m_predTs;
        Point3<T>           m_predx;        // predicted mass center
        Quaternion<T>       m_predq;        // rotation
        Quaternion<T>       m_predinvq;     // predicted inverse of quaternion
        Matrix3<T>          m_predinvI;     // inverse of predicted inertia
        Vector3<T>          m_predv;        // predicted velocity
        Vector3<T>          m_predw;        // predicted angular velocity

        REAL                m_restCoeff;    // restitution coefficient
        REAL                m_mu;           // friction coefficient
        TCollProc           m_cproc;        // collision processor

        Vector3<T>          m_extF;         // external force
        Vector3<T>          m_extT;         // external torque
        /*
         * Indicate if this object is pinned temporarily. This value is used 
         * in shock propagation.
         */
        bool                m_pinned;

#ifdef USE_RECORDER
        /* 
         * KD tree is used to find the closest vertex from the collision point.
         *
         * It is useful when recording the collision impulse at each timestep
         */
        KDTree< 3, Point3<T>, Point3DistSqr<T>, Point3Acc<T> >  m_kdtree;

    public:
        const KDTree< 3, Point3<T>, Point3DistSqr<T>, Point3Acc<T> >& kdtree() const
        {  return m_kdtree; }
        /*
         * initialize KD-tree. This method should be called explicitly before using 
         * KD-tree. 
         */
        void init_kdtree();
#endif
};

///////////////////////////////////////////////////////////////////////////////

/*
 * apply impulse on this object, changing the current velocity
 */
template <typename T, class _TMesh>
void LSCollisionRigidBody<T, _TMesh>::apply_impulse_to_prediction(
        const Vector3<T>& impl, const Vector3<T>& r)
{
    if ( m_fixed ) return;

    // change the current unpredicted velocity and angular velocity
    m_v += impl * m_invMass;    

    // when computing the angular velocity changes, use the predicted
    // position and inertia
    // m_omega += m_predinvI * (pt - m_predx).crossProduct(impl);
    m_omega += m_predinvI * r.crossProduct(impl);
}

/*
 * update predicted states by using current force
 * Once m_v/m_omega is changed, make sure pred_timestamp get increased,
 * such that m_predx, m_predinvq can get updated
 *
 * predv = v + dt*a
 * predx = x + dt(v + dt*a)
 */
template <typename T, class _TMesh>
void LSCollisionRigidBody<T, _TMesh>::update_force_predicted_position(REAL dt)
{
    if ( m_predTs >= m_cproc.pred_timestamp() ) return;
    m_predTs = m_cproc.pred_timestamp();

    // advance velocity
    m_predv.scaleAdd(dt, m_acc, m_v);
    m_predw.scaleAdd(dt, m_angAcc, m_omega);

    // advance state using the predicted velocity
    m_predx.scaleAdd(dt, m_predv, m_x);
    T omeganorm = m_predw.length();
    if ( omeganorm > EPS )
    {
        const T dang = omeganorm * dt * 0.5;
        Quaternion<T> rot(cos(dang), m_predw * (sin(dang) / omeganorm));
        m_predq = rot * m_q;
        m_predq.normalize();
    }
    else
        m_predq = m_q;

    m_predinvq = m_predq.conjugate();

    const Matrix3<T> R = m_predq.toMatrix3();
    m_predinvI = R * m_invI0 * R.transpose();
}

/*
 * update predicted states by using current velocity
 */
template <typename T, class _TMesh>
void LSCollisionRigidBody<T, _TMesh>::update_velocity_predicted_position(REAL dt)
{
    if ( m_predTs >= m_cproc.pred_timestamp() ) return;
    m_predTs = m_cproc.pred_timestamp();

    m_predv = m_v;
    m_predw = m_omega;

    // advance state using current velocity
    m_predx.scaleAdd(dt, m_predv, m_x);     // position
    T omeganorm = m_predw.length();
    if ( omeganorm > EPS )
    {   // has some rotation velocity
        const T dang = omeganorm * dt * 0.5;
        Quaternion<T> rot(cos(dang), m_predw * (sin(dang) / omeganorm));
        m_predq = rot * m_q;                // update rotation
        m_predq.normalize();
    }
    else
        m_predq = m_q;

    m_predinvq = m_predq.conjugate();

    const Matrix3<T> R = m_predq.toMatrix3();
    m_predinvI = R * m_invI0 * R.transpose();
}

#ifdef USE_RECORDER
/*
 * NOTE that the KD-tree is created based on the rest position of the
 * rigid body.
 */
template <typename T, class _TMesh>
void LSCollisionRigidBody<T, _TMesh>::init_kdtree()
{
    const std::vector< Point3<T> >& restPos = mp_mesh->rest_positions();
    m_kdtree.reinitialize(&restPos[0], restPos.size());
}
#endif

#ifdef USE_NAMESPACE
}
#endif
#endif
