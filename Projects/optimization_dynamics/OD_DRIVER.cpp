#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Log/SCOPE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <iomanip>
#include "OD_DRIVER.h"
using namespace PhysBAM;
namespace{
    template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
    {
        ((OD_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
    }
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> OD_DRIVER<TV>::
OD_DRIVER(OD_EXAMPLE<TV>& example,OD_SOLVER<TV>& solver)
    :example(example),solver(solver)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> OD_DRIVER<TV>::
~OD_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void OD_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void OD_DRIVER<TV>::
Initialize()
{
    LOG::cout<<std::setprecision(16)<<std::endl;
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.substeps_delay_frame<0?example.write_substeps_level:-1);

    // setup time
    output_number=current_frame=example.restart;

    example.Initialize();
    if(example.restart)
        example.Read_Output_Files(example.restart);


    if(!example.restart) Write_Output_Files(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",0,1);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void OD_DRIVER<TV>::
Advance_One_Time_Step()
{
    auto& X=example.particles.X;
    auto& V=example.particles.V;
    auto& mass=example.particles.mass;
    auto P(X);
    auto dt=example.dt;
    example.Begin_Time_Step(example.time);

    Apply_External_Forces();

    PHYSBAM_DEBUG_WRITE_SUBSTEP("after external forces",0,1);
    solver.Solve(dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after projecting constraints",0,1);

    for(int i=0;i<X.m;i++)
        V(i)=(X(i)-P(i))/dt;

    TV p,c;
    auto L=p.Cross(c);
    T ke=0;
    for(int i=0;i<X.m;i++){
        TV pi=mass(i)*V(i);
        p+=pi;
        L+=pi.Cross(X(i));
        c+=mass(i)*X(i);
        ke+=mass(i)*V(i).Magnitude_Squared();}
    LOG::cout<<"linear momentum: "<<p<<std::endl;
    LOG::cout<<"angular momentum: "<<L<<std::endl;
    LOG::cout<<"center of mass: "<<c<<std::endl;
    LOG::cout<<"kinetic energy: "<<ke<<std::endl;

    example.End_Time_Step(example.time);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void OD_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    for(;current_frame<frame;current_frame++){
        LOG::SCOPE scope("FRAME","frame %d",current_frame+1);
        example.Begin_Frame(current_frame);
        if(example.substeps_delay_frame==current_frame)
            DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);
        T time_at_frame=example.time+example.frame_dt;
        bool done=false;
        for(int substep=0;!done;substep++){
            LOG::SCOPE scope("SUBSTEP","substep %d",substep+1);
            example.dt=clamp(example.dt,example.min_dt,example.max_dt);
            T next_time=example.time+example.dt;
            if(next_time>time_at_frame){
                next_time=time_at_frame;
                done=true;}
            else if(next_time+example.dt>time_at_frame) next_time=(example.time+time_at_frame)/2;
            example.dt=next_time-example.time;
            LOG::cout<<"substep dt: "<<example.dt<<std::endl;

            Advance_One_Time_Step();
            example.time=next_time;}
        example.End_Frame(current_frame);
        Write_Output_Files(++output_number);}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void OD_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        LOG::printf("Writing substep [%s]: output_number=%i, time=%g, frame=%i, substep=%i\n",
            example.frame_title,output_number+1,time,current_frame,substep);
        Write_Output_Files(++output_number);
        example.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void OD_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+LOG::sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+LOG::sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame==0)
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
// Function Apply_External_Forces
//#####################################################################
template<class TV> void OD_DRIVER<TV>::
Apply_External_Forces()
{
}

namespace PhysBAM{
template class OD_DRIVER<VECTOR<float,2> >;
template class OD_DRIVER<VECTOR<float,3> >;
template class OD_DRIVER<VECTOR<double,2> >;
template class OD_DRIVER<VECTOR<double,3> >;
}
