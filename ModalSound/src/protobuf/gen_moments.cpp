#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "transfer/FMMBoundIntTransfer.hpp"
#include "protobuf/sploosh.pb.h"
#include "utils/term_msg.h"

using namespace std;

namespace fs = boost::filesystem;

static void generate_moments(const char* fileIn, const char* fileOut)
{
    char filename[256];
    FMMBoundIntTransfer<double> tran(fileIn, fileOut);
    tran.compute_moments();

    sprintf(filename, "moments.dat");
    tran.store_moments(filename);

    sploosh::FMMoments moments;
    sploosh::Vector3d  center;
    PRINT_MSG("OBJ Center: %.lf %.lf %.lf\n", tran.obj_center().x, tran.obj_center().y, tran.obj_center().z);

    moments.mutable_center()->set_x(tran.obj_center().x);
    moments.mutable_center()->set_y(tran.obj_center().y);
    moments.mutable_center()->set_z(tran.obj_center().z);
    moments.set_numexp( (uint32_t)tran.expansion_num() );
    moments.set_wavenum( tran.wave_number() );
    moments.set_mfile(filename);

    ofstream output("moments.pbuf", ios::trunc|ios::binary);
    if ( !moments.SerializeToOstream(&output) )
    {
        PRINT_ERROR("Failed to write to file moments.pbuf\n");
    }
}

static void generate_moments(const char* fileIn, const char* fileOut, 
        int id, sploosh::FMMoments* moment)
{
    char filein[256];
    char fileout[256];
    sprintf(filein, fileIn, id);
    sprintf(fileout, fileOut, id);

    FMMBoundIntTransfer<double> tran(filein, fileout);

    sprintf(fileout, "moments-%d.dat", id);
    tran.compute_moments();
    tran.store_moments(fileout);

    moment->mutable_center()->set_x(tran.obj_center().x);
    moment->mutable_center()->set_y(tran.obj_center().y);
    moment->mutable_center()->set_z(tran.obj_center().z);
    moment->set_numexp( (uint32_t)tran.expansion_num() );
    moment->set_wavenum( tran.wave_number() );
    moment->set_mfile( fileout );
}

int main(int argc, char* argv[])
{
    // Verify that the version of the library that we linked against is
    // compatible with the version of the headers we compiled against.
    GOOGLE_PROTOBUF_VERIFY_VERSION;

    if ( argc == 3 )
    {
        generate_moments(argv[1], argv[2]);
        // Optional:  Delete all global objects allocated by libprotobuf.
        google::protobuf::ShutdownProtobufLibrary();
    }
    else if ( argc == 5 )
    {
        int stId = atoi(argv[3]);
        int edId = atoi(argv[4]);

        sploosh::ModalMoments mms;
        for(int id = stId;id <= edId;++ id)
        {
            generate_moments(argv[1], argv[2], id, mms.add_moment());
        }

        ofstream output("moments.pbuf", ios::trunc|ios::binary);
        if ( !mms.SerializeToOstream(&output) )
        {
            PRINT_ERROR("Failed to write to file moments.pbuf\n");
        }
        // Optional:  Delete all global objects allocated by libprotobuf.
        google::protobuf::ShutdownProtobufLibrary();
    }
    else
    {
        cout << "Usage: " << fs::path(argv[0]).filename().string() << " [fbem input file] [fbem output file]" << endl;
        cout << "       --------- or --------- " << endl;
        cout << "       " << fs::path(argv[0]).filename().string() << " [fbem input file] [fbem output file] [start Id] [end Id]" << endl;
        cout << endl;
        cout << "This program generates output files in the current directory" << endl;
    }
    return 0;
}
