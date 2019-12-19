//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_DRIVER.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_STATE.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SHALLOW_WATER_DRIVER<TV>::
SHALLOW_WATER_DRIVER(SHALLOW_WATER_STATE<TV>& state)
    :state(state)
{
    DEBUG_SUBSTEPS::writer=[=](const std::string& title){Write_Substep(title);};
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SHALLOW_WATER_DRIVER<TV>::
~SHALLOW_WATER_DRIVER()
{
    DEBUG_SUBSTEPS::writer=0;
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void SHALLOW_WATER_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(state.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void SHALLOW_WATER_DRIVER<TV>::
Initialize()
{
    LOG::cout<<std::setprecision(16)<<std::endl;
    DEBUG_SUBSTEPS::write_substeps_level=state.substeps_delay_frame<0?state.write_substeps_level:-1;

    // setup time
    current_frame=state.restart;
    state.viewer_dir.Set(state.restart);
    state.U.Resize(state.grid.Node_Indices(state.ghost));
    PHYSBAM_ASSERT(state.initialize);
    state.initialize();

    if(state.restart) state.Read_Output_Files();

    if(!state.restart) Write_Output_Files();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",1);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void SHALLOW_WATER_DRIVER<TV>::
Advance_One_Time_Step()
{
    if(state.begin_time_step) state.begin_time_step(state.time);

    if(state.end_time_step) state.end_time_step(state.time);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void SHALLOW_WATER_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    for(;current_frame<frame;current_frame++){
        LOG::SCOPE scope("FRAME","frame %d",current_frame+1);
        if(state.begin_frame) state.begin_frame(current_frame);
        if(state.substeps_delay_frame==current_frame)
            DEBUG_SUBSTEPS::write_substeps_level=state.write_substeps_level;
        T time_at_frame=state.time+state.frame_dt;
        bool done=false;
        for(int substep=0;!done;substep++){
            LOG::SCOPE scope("SUBSTEP","substep %d",substep+1);
            state.dt=Compute_Dt();
            state.dt=clamp(state.dt,state.min_dt,state.max_dt);
            T next_time=state.time+state.dt;
            if(next_time>time_at_frame){
                next_time=time_at_frame;
                done=true;}
            else if(next_time+state.dt>time_at_frame) next_time=(state.time+time_at_frame)/2;
            state.dt=next_time-state.time;
            LOG::cout<<"substep dt: "<<state.dt<<std::endl;

            Advance_One_Time_Step();
            state.time=next_time;}
        if(state.end_frame) state.end_frame(current_frame);
        Write_Output_Files();}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void SHALLOW_WATER_DRIVER<TV>::
Write_Substep(const std::string& title)
{
    state.frame_title=title;
    LOG::printf("Writing substep [%s]: output_number=%P, time=%g, frame=%i\n",
        state.frame_title,state.viewer_dir.frame_stack,state.time,current_frame);
    Write_Output_Files();
    state.frame_title="";
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void SHALLOW_WATER_DRIVER<TV>::
Write_Output_Files()
{
    LOG::SCOPE scope("Write_Output_Files");
    state.viewer_dir.Start_Directory(0,state.frame_title);
    state.frame_title="";
    if(state.viewer_dir.First_Frame())
        Write_To_File(state.stream_type,state.viewer_dir.output_directory+"/common/grid",state.grid);
    state.Write_Output_Files();

    ARRAY<T,TV_INT> h(state.U.domain);
    ARRAY<TV,TV_INT> v(state.U.domain);
    for(RANGE_ITERATOR<TV::m> it(state.U.domain);it.Valid();it.Next()){
        VECTOR<T,TV::m+1> UU=state.U(it.index);
        T a=UU(0);
        h(it.index)=a;
        v(it.index)=UU.template Slice<1,TV::m>()/a;}
    Write_To_File(state.stream_type,state.viewer_dir.current_directory+"/centered_velocities",v);
    Write_To_File(state.stream_type,state.viewer_dir.current_directory+"/heightfield",h);
    state.debug_particles.Write_Debug_Particles(state.stream_type,state.viewer_dir);
    state.viewer_dir.Finish_Directory();
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR SHALLOW_WATER_DRIVER<TV>::
Compute_Dt() const
{
    return state.cfl*state.shallow_water.CFL();
}
//#####################################################################
template class SHALLOW_WATER_DRIVER<VECTOR<float,1> >;
template class SHALLOW_WATER_DRIVER<VECTOR<double,1> >;
template class SHALLOW_WATER_DRIVER<VECTOR<float,2> >;
template class SHALLOW_WATER_DRIVER<VECTOR<double,2> >;
}
