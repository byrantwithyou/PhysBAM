//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STRAWMAN_DRIVER
//#####################################################################
#include "STRAWMAN_DRIVER.h"
#include "STRAWMAN_EXAMPLE.h"

namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> STRAWMAN_DRIVER<TV>::
STRAWMAN_DRIVER(STRAWMAN_EXAMPLE<TV>& example)
: BASE(example),example(example)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STRAWMAN_DRIVER<TV>::
~STRAWMAN_DRIVER()
{}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void STRAWMAN_DRIVER<TV>::
Initialize()
{
    BASE::Initialize();
    example.After_Initialization();
    if(!example.restart) Write_Output_Files(example.first_frame);
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void STRAWMAN_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    if((current_frame-example.first_frame)%20==0){
        example.pls_evolution->Reseed_Particles(time);
        example.pls_evolution->Delete_Particles_Outside_Grid();}

    bool done=false;
    for(int substeps=1;!done;++substeps){
        T dt=example.Compute_Dt(time,target_time,done);
        example.Euler_Step(dt,time);
        LOG::cout<<"Substepping from "<<time<<" with dt = "<<dt<<std::endl;
        time+=dt;}
}
//#####################################################################
template class STRAWMAN_DRIVER<VECTOR<float,2> >;
template class STRAWMAN_DRIVER<VECTOR<double,2> >;
}
