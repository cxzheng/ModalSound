#include "FMMTransferEval.h"
#include <fstream>
#include <boost/filesystem.hpp>
#include "utils/term_msg.h"
#include "utils/macros.h"
#include "sploosh.pb.h"
#include "io/MatrixIO.hpp"

using namespace std;

FMMTransferEval::FMMTransferEval(const char* d) : moments_(NULL)
{
    namespace fs = boost::filesystem;

    fs::path dir(d);
    fs::path protofile = dir / "moments.pbuf";
    fs::path datafile  = dir / "moments.dat";

    if ( !fs::exists(protofile) || !fs::exists(datafile) )
        SHOULD_NEVER_HAPPEN(-1);

    sploosh::FMMoments moments;
    ifstream fin(protofile.string().c_str(), ios::binary);
    if ( !moments.ParseFromIstream(&fin) )
        SHOULD_NEVER_HAPPEN(-1);

    nexpan_  = (int)moments.numexp();
    waveNum_ = moments.wavenum();
    center_.set( moments.center().x(), moments.center().y(), moments.center().z() );
    moments_ = load_ma_matrixd( datafile.string().c_str() );
}

/*
 * Load the moment data, assuming the filename is given in the FMMoments instance,
 * and the file is stored in the given directory.
 */
FMMTransferEval::FMMTransferEval(const sploosh::FMMoments& ms, const string& d) : moments_(NULL)
{
    namespace fs = boost::filesystem;
    fs::path dir(d);

    nexpan_  = (int)ms.numexp();
    waveNum_ = ms.wavenum();
    center_.set( ms.center().x(), ms.center().y(), ms.center().z() );

    fs::path datafile  = dir / ms.mfile();
    moments_ = load_ma_matrixd( datafile.string().c_str() );
    if ( !moments_ ) 
    {
        PRINT_ERROR("Failed to load the moment file: %s\n", datafile.string().c_str() );
        SHOULD_NEVER_HAPPEN(-4);
    }
    PRINT_MSG("load moment file: %s\n", datafile.string().c_str());
}
