#ifndef MODAL_SHAPE_INC
#   define MODAL_SHAPE_INC

#include <vector>

class ModalShape
{

    public:
        ModalShape(const char* file);

        int num_modes() const
        {   return numModes_; }

        double eigen_mode(int i) const
        {   return eigenmodes_[i]; }

    private:
        int     n3_;
        int     numModes_;

        std::vector<double>     eigenvec_;
        std::vector<double>     eigenmodes_;

};

#endif
