#ifndef __TEST_COMMON__
#define __TEST_COMMON__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "SIM_COMMON.h"

using namespace PhysBAM;

template<class TV>
struct TEST_COMMON
{
    typedef typename TV::SCALAR T;
    std::string output_directory;
    PARSE_ARGS parse_args;
    SIM_COMMON<TV> sim;

    TEST_COMMON();
    ~TEST_COMMON();
    void Init_1();
    void Init_2(int argc,char** argv);
    void Init_3();
};

#endif
