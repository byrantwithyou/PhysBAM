#ifndef __OBJECTS_COMMON__
#define __OBJECTS_COMMON__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include "HEADER.h"

using namespace PhysBAM;

template<class TV> class BOUNDARY_CONDITIONS;

template<class TV>
struct OBJECTS_COMMON
{
    OBJECTS_COMMON();
    ~OBJECTS_COMMON();

    GRID<TV> grid;
    BOUNDARY_CONDITIONS<TV>* bc;
    ACCURACY_INFO<TV::m> ai;
};

#endif
