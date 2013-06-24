#ifndef __PARAMETERS_COMMON__
#define __PARAMETERS_COMMON__

#include <Tools/Parsing/PARSE_ARGS.h>
#include "HEADER.h"

using namespace PhysBAM;

template<class T>
struct PARAMETERS_COMMON
{
    PARAMETERS_COMMON();

    T time,dt,rho,mu;
    T theta_threshold,cg_tolerance;
    bool print_matrix;

    void Init_1(PARSE_ARGS& parse_args);
    void Init_2(PARSE_ARGS& parse_args);
};

#endif
