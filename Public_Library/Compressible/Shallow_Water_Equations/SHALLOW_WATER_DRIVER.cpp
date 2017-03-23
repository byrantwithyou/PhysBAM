//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_DRIVER.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_STATE.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace PhysBAM{
namespace{
    template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
    {
        ((SHALLOW_WATER_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
    }
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SHALLOW_WATER_DRIVER<TV>::
SHALLOW_WATER_DRIVER(SHALLOW_WATER_STATE<TV>& state)
    :state(state)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SHALLOW_WATER_DRIVER<TV>::
~SHALLOW_WATER_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
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
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(state.substeps_delay_frame<0?state.write_substeps_level:-1);

    // setup time
    output_number=current_frame=state.restart;

    PHYSBAM_ASSERT(state.initialize);
    state.initialize();

    if(state.restart) state.Read_Output_Files(state.restart);

    if(!state.restart) Write_Output_Files(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",0,1);
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
            DEBUG_SUBSTEPS::Set_Write_Substeps_Level(state.write_substeps_level);
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
        Write_Output_Files(++output_number);}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void SHALLOW_WATER_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=state.write_substeps_level){
        state.frame_title=title;
        LOG::printf("Writing substep [%s]: output_number=%i, time=%g, frame=%i, substep=%i\n",
            state.frame_title,output_number+1,state.time,current_frame,substep);
        Write_Output_Files(++output_number);
        state.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void SHALLOW_WATER_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    LOG::SCOPE scope("Write_Output_Files");
    FILE_UTILITIES::Create_Directory(state.output_directory);
    FILE_UTILITIES::Create_Directory(state.output_directory+LOG::sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(state.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(state.output_directory+LOG::sprintf("/%d/frame_title",frame),state.frame_title);
    if(frame==0)
        FILE_UTILITIES::Write_To_Text_File(state.output_directory+"/common/first_frame",frame,"\n");
    state.Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_Text_File(state.output_directory+"/common/last_frame",frame,"\n");
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
