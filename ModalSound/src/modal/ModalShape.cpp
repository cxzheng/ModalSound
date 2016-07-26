#include "ModalShape.h"
#include <fstream>
#include "utils/term_msg.h"
#include "utils/macros.h"

using namespace std;

/*
 * Load Eigen modes
 */
ModalShape::ModalShape(const char* file)
{
    ifstream fin(file, ios::binary);

    fin.read((char *)&n3_, sizeof(int));       // size of eigen problem
    fin.read((char *)&numModes_, sizeof(int)); // number of modes
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot read file: %s\n", file);
        SHOULD_NEVER_HAPPEN(2);
    }

    eigenmodes_.resize(numModes_);
    PRINT_MSG("Load %d eigenmodes\n", numModes_);
    fin.read((char *)&eigenmodes_[0], sizeof(double)*numModes_);

    //// all the eigen vectors are stored in a n3 x nModes matrix
    //// it is stored as v1 v2 ...
    eigenvec_.resize(n3_*numModes_);
    fin.read((char *)&eigenvec_[0], sizeof(double)*eigenvec_.size());
    if ( fin.fail() )
    {
        PRINT_ERROR("Cannot read file successfully: %s\n", file);
        SHOULD_NEVER_HAPPEN(2);
    }
}

