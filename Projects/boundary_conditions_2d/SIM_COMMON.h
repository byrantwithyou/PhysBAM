#ifndef __SIM_COMMON__
#define __SIM_COMMON__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include "OBJECTS_COMMON.h"
#include "PARAMETERS_COMMON.h"

using namespace PhysBAM;

template<class TV>
struct SIM_COMMON
{
    typedef typename TV::SCALAR T;
    SIM_COMMON();

    PARAMETERS_COMMON<T> param;
    OBJECTS_COMMON<TV> obj;
    int resolution,steps;
    int base_resolution;
    proj_type proj_algo;
    bool check_leaks;
    bool use_accuracy_samples;
    bool use_extrapolation;
    bool use_viscosity,use_projection,use_advection;
    bool use_proj_slip,no_gibou;

    void Init_1(PARSE_ARGS& parse_args);
    void Init_2(PARSE_ARGS& parse_args);
    void Init_3();
};

#endif
