//#####################################################################
// Copyright 2017, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//##################################################################### 
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_STATE.h>
#include "COMMON_DATA.h"

namespace PhysBAM
{
template<class TV> void
Setup_Common(SHALLOW_WATER_STATE<TV>& st,PARSE_ARGS& parse_args,COMMON_DATA<TV>& cd)
{
    typedef typename TV::SCALAR T;
    T framerate=24;
    parse_args.Extra(&cd.test_number,"example number","example number to run");
    parse_args.Add("-restart",&st.restart,"frame","restart frame");
    parse_args.Add("-resolution",&cd.resolution,&cd.user_resolution,"resolution","grid resolution");
    parse_args.Add("-substeps",&st.write_substeps_level,"level","output-substep level");
    parse_args.Add("-substeps_delay",&st.substeps_delay_frame,"frame","delay substeps until after this frame");
    parse_args.Add("-last_frame",&st.last_frame,&cd.user_last_frame,"frame","number of frames to simulate");
    parse_args.Add("-threads",&cd.threads,"threads","Number of threads");
    parse_args.Add("-o",&st.output_directory,"dir","Output directory");
    parse_args.Add("-framerate",&framerate,"rate","Number of frames per second");
    parse_args.Add("-min_dt",&st.min_dt,"dt","Minimum time step size");
    parse_args.Add("-max_dt",&st.max_dt,"dt","Maximum time step size");
    parse_args.Add("-cfl",&st.cfl,"cfl","CFL number");
    parse_args.Add("-seed",&cd.seed,"seed","Random number seed");
    parse_args.Add("-m",&cd.m,"scale","meter scale");
    parse_args.Add("-s",&cd.s,"scale","second scale");
    parse_args.Add("-kg",&cd.kg,"scale","kilogram scale");
    parse_args.Add("-d",&st.data_directory,"dir","data directory");
    parse_args.Add("-test_output_prefix",&st.test_output_prefix,&st.use_test_output,"","prefix to use for test output");
    parse_args.Parse(true);
    st.frame_dt=cd.s/framerate;
}
template void Setup_Common<VECTOR<double,1> >(
    SHALLOW_WATER_STATE<VECTOR<double,1> >&,PARSE_ARGS&,
    COMMON_DATA<VECTOR<double,1> >&);
template void Setup_Common<VECTOR<double,2> >(
    SHALLOW_WATER_STATE<VECTOR<double,2> >&,PARSE_ARGS&,
    COMMON_DATA<VECTOR<double,2> >&);
template void Setup_Common<VECTOR<float,1> >(
    SHALLOW_WATER_STATE<VECTOR<float,1> >&,PARSE_ARGS&,
    COMMON_DATA<VECTOR<float,1> >&);
template void Setup_Common<VECTOR<float,2> >(
    SHALLOW_WATER_STATE<VECTOR<float,2> >&,PARSE_ARGS&,
    COMMON_DATA<VECTOR<float,2> >&);
}
