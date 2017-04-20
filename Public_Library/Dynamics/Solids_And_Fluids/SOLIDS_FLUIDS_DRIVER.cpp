//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_DRIVER
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_FLUIDS_DRIVER<TV>::
SOLIDS_FLUIDS_DRIVER(SOLIDS_FLUIDS_EXAMPLE<TV>& example_input)
    :BASE(example_input),example(example_input),project_at_frame_boundaries(true)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_FLUIDS_DRIVER<TV>::
~SOLIDS_FLUIDS_DRIVER()
{}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER<TV>::
Execute_Main_Program()
{
    {LOG::SCOPE scope("INITIALIZING","Initializing");
    Initialize();
    example.Post_Initialization();
    example.Log_Parameters();
    if(!example.restart) Write_Output_Files(example.first_frame);}
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER<TV>::
Initialize()
{
    if(example.auto_restart){
        std::string last_frame_file=example.output_directory+"/common/last_frame";
        int last_frame;Read_From_Text_File(last_frame_file,last_frame);
        example.restart=true;example.restart_frame=last_frame;
        LOG::cout<<"Auto Restart from frame "<<last_frame<<" (from file "<<last_frame_file<<")"<<std::endl;}
    if(example.restart){current_frame=example.restart_frame;Read_Time(current_frame);}else current_frame=example.first_frame;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);
    example.solid_body_collection.deformable_body_collection.mpi_solids=example.solid_body_collection.deformable_body_collection.mpi_solids;
}
//#####################################################################
// Function Simulate_To_Frame
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_DRIVER<TV>::
Simulate_To_Frame(const int frame_input)
{
    while(current_frame<frame_input){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        Preprocess_Frame(current_frame+1);
        example.solids_evolution->kinematic_evolution.Get_Current_Kinematic_Keyframes(example.Time_At_Frame(current_frame+1)-time,time);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        Postprocess_Frame(++current_frame);
        if(example.write_output_files && example.write_substeps_level==-1) Write_Output_Files(current_frame);
        else if(example.write_substeps_level!=-1) Write_Substep(LOG::sprintf("END Frame %d",current_frame),0,example.write_substeps_level);
        LOG::cout<<"TIME = "<<time<<std::endl;}
}
//#####################################################################
namespace PhysBAM{
template class SOLIDS_FLUIDS_DRIVER<VECTOR<float,1> >;
template class SOLIDS_FLUIDS_DRIVER<VECTOR<float,2> >;
template class SOLIDS_FLUIDS_DRIVER<VECTOR<float,3> >;
template class SOLIDS_FLUIDS_DRIVER<VECTOR<double,1> >;
template class SOLIDS_FLUIDS_DRIVER<VECTOR<double,2> >;
template class SOLIDS_FLUIDS_DRIVER<VECTOR<double,3> >;
}
