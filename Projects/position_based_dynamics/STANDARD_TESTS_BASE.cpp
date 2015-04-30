//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Parsing/PARSE_ARGS.h>
#include "STANDARD_TESTS_BASE.h"
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> STANDARD_TESTS_BASE<TV>::
STANDARD_TESTS_BASE(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
    :PBD_EXAMPLE<TV>(stream_type),test_number(0),resolution(32),user_resolution(false),stored_last_frame(0),
    user_last_frame(false),seed(1234),scale_mass(1),scale_stiffness(1),scale_speed(1)
{
    T framerate=24;
    parse_args.Extra(&test_number,"example number","example number to run");
    parse_args.Add("-restart",&restart,"frame","restart frame");
    parse_args.Add("-resolution",&resolution,&user_resolution,"resolution","grid resolution");
    parse_args.Add("-substeps",&write_substeps_level,"level","output-substep level");
    parse_args.Add("-substeps_delay",&substeps_delay_frame,"frame","delay substeps until after this frame");
    parse_args.Add("-last_frame",&last_frame,&user_last_frame,"frame","number of frames to simulate");
    parse_args.Add("-test_diff",&test_diff,"test analytic derivatives");
    parse_args.Add("-threads",&threads,"threads","Number of threads");
    parse_args.Add("-o",&output_directory,"dir","Output directory");
    parse_args.Add("-framerate",&framerate,"rate","Number of frames per second");
    parse_args.Add("-min_dt",&min_dt,"dt","Minimum time step size");
    parse_args.Add("-max_dt",&max_dt,"dt","Maximum time step size");
    parse_args.Add("-print_stats",&print_stats,"Print momentum/energy stats");
    parse_args.Add("-solver_iterations",&solver_iterations,"iter","Solver iterations");
    parse_args.Add("-seed",&seed,"seed","Random number seed");
    parse_args.Add("-scale_mass",&scale_mass,"scale","Scale mass of particles");
    parse_args.Add("-scale_stiffness",&scale_stiffness,"scale","Scale stiffness");
    parse_args.Add("-scale_speed",&scale_speed,"scale","Scale initial speed of simulated object");
    parse_args.Parse(true);

    frame_dt=1/framerate;

#ifdef USE_OPENMP
    omp_set_num_threads(threads);
#pragma omp parallel
#pragma omp single
    {
        if(omp_get_num_threads()!=threads) PHYSBAM_FATAL_ERROR();
        LOG::cout<<"Running on "<<threads<<" threads"<<std::endl;
    }
#endif

    stored_last_frame=last_frame;
    random.Set_Seed(seed);

}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STANDARD_TESTS_BASE<TV>::
~STANDARD_TESTS_BASE()
{}

template class STANDARD_TESTS_BASE<VECTOR<float,2> >;
template class STANDARD_TESTS_BASE<VECTOR<float,3> >;
template class STANDARD_TESTS_BASE<VECTOR<double,2> >;
template class STANDARD_TESTS_BASE<VECTOR<double,3> >;
}
