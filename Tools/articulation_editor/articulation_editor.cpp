#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include "ARTICULATION_VISUALIZATION.h"
using namespace PhysBAM;
// ==========================================================================
int main(int argc, char *argv[])
{
    bool type_double = false;   // float by default
    if (PARSE_ARGS::Find_And_Remove("-float", argc, argv))
        type_double = false;
    if (PARSE_ARGS::Find_And_Remove("-double", argc, argv))
        type_double = true;

    ANIMATED_VISUALIZATION *visualization = 0;
    if(!type_double) visualization=new ARTICULATION_VISUALIZATION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    //else visualization=new ARTICULATION_VISUALIZATION<double>;
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
    visualization->Initialize_And_Run(argc, argv);

    return 0;
}
