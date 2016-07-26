#ifndef SC_DIAGONAL_MATRIX_HPP
#   define SC_DIAGONAL_MATRIX_HPP

#include <vector>

/*!
 * Diagonal matrix. i.e. The mass matrix in unreduced model
 */
template <typename T>
class DiagonalMatrix
{
    public:
        int size() const { return data_.size(); }

        DiagonalMatrix<T>& operator = (const std::vector<T>& vec)
        {
            data_ = vec;
            return *this;
        }

        void resize(int n) 
        {
            data_.resize(n); 
        }

        T operator [] (int n) const
        {
            assert(n >= 0 && n < (int)data_.size());
            return data_[n];
        }

        T& operator [] (int n) 
        {
            assert(n >= 0 && n < (int)data_.size());
            return data_[n];
        }

        void multiply(const std::vector<T>& in, std::vector<T>& out) const
        {
            assert(in.size() == data_.size());
            if ( out.size() != in.size() ) out.resize(in.size());

#ifdef USE_OPENMP
            #pragma omp parallel for default(none) schedule(dynamic, 5000) shared(out, in)
#endif
            for(int i = 0;i < (int)data_.size();++ i)
                out[i] = in[i] * data_[i];
        }

        void multiply(const T* in, T* out) const
        {
#ifdef USE_OPENMP
            #pragma omp parallel for default(none) schedule(dynamic, 5000) shared(out, in)
#endif
            for(int i = 0;i < (int)data_.size();++ i)
                out[i] = in[i] * data_[i];
        }

        const T* data() const
        { return &data_[0]; }

    private:
        std::vector<T>  data_;
};

#endif
