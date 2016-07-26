/*
 * =====================================================================================
 *
 *       Filename:  ModalModel.h
 *
 *    Description:  Load & maintain the basic modal model
 *
 *        Version:  1.0
 *        Created:  07/10/16 17:15:47
 *       Revision:  none
 *       Compiler:  icpc
 *
 *         Author:  Changxi Zheng, cxz@cs.columbia.edu
 *                  Columbia University
 *
 * =====================================================================================
 */
#ifndef  MODAL_MODEL_INC
#   define  MODAL_MODEL_INC

#include <vector>
#include "sc/Vector3.hpp"

/*
 * Load, maintain and manipulate the linear modal model
 */
class ModalModel
{
    public:
        ModalModel(const std::string& filename, double density, double alpha, double beta);

        void accum_modal_impulse(int vid, const Vector3d* imp, double* out) const;

        /* ------- getter methods ------- */
        int num_modes() const
        {   return numModes_; }
        int len_eigvec() const      // length of eigen vector
        {   return n3_; }
        double inv_density() const
        {   return invDensity_; }
        const std::vector<double>& omega() const
        {   return omega_; }
        const std::vector<double>& damped_omega() const
        {   return omegaD_; }
        const std::vector<double>& damping_vector() const
        {   return c_; }
        const std::vector<double>& eigenvec() const
        {   return eigenvec_; }    // eigenvector in column order
        const double* eigenvec_ptr() const
        {   return &eigenvec_[0]; }

        // Return pointer to the shape vector of a mode at given vertex
        const double* shape_vec_ptr(int vid, int modeId) const
        {   return &eigenvec_[modeId*n3_ + vid*3]; }

        double omega(int mid) const
        {   
            assert(mid >= 0 && mid < numModes_);
            return omega_[mid];
        }

        double rayleigh_damping_coeff(int mid) const
        {   return c_[mid]; } 
        
    private:
        void load_eigenmodes(const char* file);

    protected:
        double      density_;
        double      invDensity_;    // 1./density
        double      alpha_;         // Rayleigh damping parameters
        double      beta_;

        int         n3_;            // length of each eigen-vector
        int         numModes_;      // number of modes

        /* v1_1 v1_2 ... v1_n v2_1 v2_2 ... */
        std::vector<double>     eigenvec_;      // eigenvector in column order
        std::vector<double>     eigenmodes_;

        std::vector<double>     omega_;         // natural frequency
        std::vector<double>     omegaD_;        // damped natural frequency
        std::vector<double>     freqs_;         // frequency values
        std::vector<double>     c_;
};
#endif   /* ----- #ifndef MODAL_MODEL_INC  ----- */

