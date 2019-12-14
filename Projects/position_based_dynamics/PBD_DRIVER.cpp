#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include "PBD_CONSTRAINTS.h"
#include "PBD_DRIVER.h"
#include <iomanip>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PBD_DRIVER<TV>::
PBD_DRIVER(PBD_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::writer=[=](const std::string& title){Write_Substep(title);};
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PBD_DRIVER<TV>::
~PBD_DRIVER()
{
    DEBUG_SUBSTEPS::writer=0;
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void PBD_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void PBD_DRIVER<TV>::
Initialize()
{
    LOG::cout<<std::setprecision(16)<<std::endl;
    DEBUG_SUBSTEPS::write_substeps_level=example.substeps_delay_frame<0?example.write_substeps_level:-1;

    // setup time
    output_number=current_frame=example.restart;

    example.Initialize();
    if(example.restart)
        example.Read_Output_Files();

    example.P.Resize(example.X.m);

    if(!example.restart) Write_Output_Files();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",1);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void PBD_DRIVER<TV>::
Advance_One_Time_Step()
{
    auto& X=example.X;
    auto& P=example.P;
    auto& V=example.V;
    auto& w=example.w;
    auto dt=example.dt;
    example.Begin_Time_Step(example.time);

    Apply_External_Forces();

    for(int i=0;i<X.m;i++)
        P(i)=X(i)+dt*V(i);

    PHYSBAM_DEBUG_WRITE_SUBSTEP("after external forces",1);
    for(int i=0;i<example.solver_iterations;i++)
        Project_Constraints();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after projecting constraints",1);

    for(int i=0;i<X.m;i++){
        V(i)=(P(i)-X(i))/dt;
        X(i)=P(i);}

    TV p,c;
    auto L=p.Cross(c);
    T ke=0;
    for(int i=0;i<X.m;i++){
        TV pi=V(i)/w(i);
        p+=pi;
        L+=pi.Cross(X(i));
        c+=X(i)/w(i);
        ke+=V(i).Magnitude_Squared()/w(i);}
    LOG::cout<<"linear momentum: "<<p<<std::endl;
    LOG::cout<<"angular momentum: "<<L<<std::endl;
    LOG::cout<<"center of mass: "<<c<<std::endl;
    LOG::cout<<"kinetic energy: "<<ke<<std::endl;

    example.End_Time_Step(example.time);
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void PBD_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    for(;current_frame<frame;current_frame++){
        LOG::SCOPE scope("FRAME","frame %d",current_frame+1);
        example.Begin_Frame(current_frame);
        if(example.substeps_delay_frame==current_frame)
            DEBUG_SUBSTEPS::write_substeps_level=example.write_substeps_level;
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
        Write_Output_Files();}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void PBD_DRIVER<TV>::
Write_Substep(const std::string& title)
{
    example.frame_title=title;
    LOG::printf("Writing substep [%s]: output_number=%i, time=%g, frame=%i\n",
        example.frame_title,output_number+1,example.time,current_frame);
    Write_Output_Files();
    example.frame_title="";
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void PBD_DRIVER<TV>::
Write_Output_Files()
{
    example.viewer_dir.Start_Directory(0,example.frame_title);
    example.frame_title="";
    example.Write_Output_Files();
    example.viewer_dir.Finish_Directory();
}
//#####################################################################
// Function Apply_External_Forces
//#####################################################################
template<class TV> void PBD_DRIVER<TV>::
Apply_External_Forces()
{
}
//#####################################################################
// Function Project_Contraints
//#####################################################################
template<class TV> void PBD_DRIVER<TV>::
Project_Constraints()
{
    for(int i=0;i<example.constraints.m;i++){
        if(example.test_diff) example.constraints(i)->Test_Diff(example.P);
        example.constraints(i)->Project(example.P,example.w,example.solver_iterations);}
}

namespace PhysBAM{
template class PBD_DRIVER<VECTOR<float,2> >;
template class PBD_DRIVER<VECTOR<float,3> >;
template class PBD_DRIVER<VECTOR<double,2> >;
template class PBD_DRIVER<VECTOR<double,3> >;
}
