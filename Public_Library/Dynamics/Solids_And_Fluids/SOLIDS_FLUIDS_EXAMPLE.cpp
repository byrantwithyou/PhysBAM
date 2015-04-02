//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_EXAMPLE
//#####################################################################
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Log/SCOPE.h>
#include <Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_PARAMETERS.h>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_FLUIDS_EXAMPLE<TV>::
SOLIDS_FLUIDS_EXAMPLE(const STREAM_TYPE stream_type)
    :BASE(stream_type),use_melting(false),solids_parameters(*new SOLIDS_PARAMETERS<TV>),solids_fluids_parameters(*new SOLIDS_FLUIDS_PARAMETERS<TV>(this)),
    solid_body_collection(*new SOLID_BODY_COLLECTION<TV>),solids_evolution(new NEWMARK_EVOLUTION<TV>(solids_parameters,solid_body_collection,*this)),
    opt_solidssymmqmr(false),opt_solidscr(false),opt_solidscg(false)
{
    Set_Minimum_Collision_Thickness();
    Set_Write_Substeps_Level(-1);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_FLUIDS_EXAMPLE<TV>::
~SOLIDS_FLUIDS_EXAMPLE()
{
    delete solids_evolution;
    delete &solid_body_collection;
    delete &solids_parameters;
    delete &solids_fluids_parameters;
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Read_Output_Files_Solids(const int frame)
{
    solid_body_collection.Read(stream_type,output_directory,frame,frame,solids_parameters.write_static_variables_every_frame,solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,
        solids_parameters.write_deformable_body,solids_parameters.write_from_every_process);
    std::string f=LOG::sprintf("%d",frame);
    //if(NEWMARK_EVOLUTION<TV>* newmark=dynamic_cast<NEWMARK_EVOLUTION<TV>*>(solids_evolution))
    //    newmark->Read_Position_Update_Projection_Data(stream_type,output_directory+"/"+f+"/");
}
//#####################################################################
// Function Log_Parameters
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Log_Parameters() const
{
    LOG::SCOPE scope("SOLIDS_FLUIDS_EXAMPLE parameters");
    BASE::Log_Parameters();
    LOG::cout<<"minimum_collision_thickness="<<minimum_collision_thickness<<std::endl;
    LOG::cout<<"use_melting="<<use_melting<<std::endl;
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Register_Options()
{
    if(!parse_args) return;
    BASE::Register_Options();
    parse_args->Add("-solidscfl",&solids_parameters.cfl,"cfl","solids CFL");
    parse_args->Add("-solidscg",&opt_solidscg,"Use CG for time integration");
    parse_args->Add("-solidscr",&opt_solidscr,"Use CONJUGATE_RESIDUAL for time integration");
    parse_args->Add("-solidssymmqmr",&opt_solidssymmqmr,"Use SYMMQMR for time integration");
    parse_args->Add("-rigidcfl",&solids_parameters.rigid_body_evolution_parameters.rigid_cfl,"cfl","rigid CFL");
}
//#####################################################################
// Function Parse_Late_Options
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Parse_Late_Options()
{
    if(!parse_args) return;
    BASE::Parse_Late_Options();

    if(opt_solidscg) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cg;
    if(opt_solidscr) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cr;
    if(opt_solidssymmqmr) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_symmqmr;
}
//#####################################################################
// Function Adjust_Output_Directory_For_MPI
//#####################################################################
template<class TV> template<class T_MPI> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Adjust_Output_Directory_For_MPI(const T_MPI mpi)
{
    if(mpi && mpi->Number_Of_Processors()>1){
        output_directory+=LOG::sprintf("/%d",(mpi->rank+1));
        FILE_UTILITIES::Create_Directory(output_directory);
        FILE_UTILITIES::Create_Directory(output_directory+"/common");
        LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",restart);}
}
//#####################################################################
// Function Post_Initialization
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Post_Initialization()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Preprocess_Frame(const int frame)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Postprocess_Frame(const int frame)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Preprocess_Substep(const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Postprocess_Substep(const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Read_Output_Files_Fluids
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Read_Output_Files_Fluids(const int frame)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Initialize_Bodies() 
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Extrapolate_Phi_Into_Objects
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Extrapolate_Phi_Into_Objects(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Postprocess_Phi
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Postprocess_Phi(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
template<class TV> bool SOLIDS_FLUIDS_EXAMPLE<TV>::
Adjust_Phi_With_Sources(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return false;
}
//#####################################################################
// Function Adjust_Phi_With_Objects
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Adjust_Phi_With_Objects(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Add_SPH_Particles_For_Sources
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Add_SPH_Particles_For_Sources(const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Initialize_SPH_Particles
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Initialize_SPH_Particles()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Initialize_Velocities()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Initialize_Euler_State
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Initialize_Euler_State()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Setup_Initial_Refinement
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Setup_Initial_Refinement()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Initialize_Advection()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Clamp_Velocities
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Clamp_Velocities(const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Melting_Substep
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Melting_Substep(const T dt,const T time)
{
}
//#####################################################################
// Function Modify_Fluid_For_Melting
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Modify_Fluid_For_Melting(const T dt,const T time)
{
}
//#####################################################################
// Function Update_Melting_Substep_Parameters
//#####################################################################
template<class TV> void SOLIDS_FLUIDS_EXAMPLE<TV>::
Update_Melting_Substep_Parameters(const T dt,const T time)
{
}
//#####################################################################
namespace PhysBAM{
template class SOLIDS_FLUIDS_EXAMPLE<VECTOR<float,1> >;
template class SOLIDS_FLUIDS_EXAMPLE<VECTOR<float,2> >;
template class SOLIDS_FLUIDS_EXAMPLE<VECTOR<float,3> >;
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<float,1> >::Adjust_Output_Directory_For_MPI<MPI_SOLIDS<VECTOR<float,1> >*>(MPI_SOLIDS<VECTOR<float,1> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<float,2> >::Adjust_Output_Directory_For_MPI<MPI_SOLIDS<VECTOR<float,2> >*>(MPI_SOLIDS<VECTOR<float,2> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<float,3> >::Adjust_Output_Directory_For_MPI<MPI_SOLIDS<VECTOR<float,3> >*>(MPI_SOLIDS<VECTOR<float,3> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<float,1> >::Adjust_Output_Directory_For_MPI<MPI_UNIFORM_GRID<VECTOR<float,1> >*>(MPI_UNIFORM_GRID<VECTOR<float,1> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<float,2> >::Adjust_Output_Directory_For_MPI<MPI_UNIFORM_GRID<VECTOR<float,2> >*>(MPI_UNIFORM_GRID<VECTOR<float,2> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<float,3> >::Adjust_Output_Directory_For_MPI<MPI_UNIFORM_GRID<VECTOR<float,3> >*>(MPI_UNIFORM_GRID<VECTOR<float,3> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<float,1> >::Adjust_Output_Directory_For_MPI<MPI_SOLID_FLUID<VECTOR<float,1> >*>(MPI_SOLID_FLUID<VECTOR<float,1> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<float,2> >::Adjust_Output_Directory_For_MPI<MPI_SOLID_FLUID<VECTOR<float,2> >*>(MPI_SOLID_FLUID<VECTOR<float,2> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<float,3> >::Adjust_Output_Directory_For_MPI<MPI_SOLID_FLUID<VECTOR<float,3> >*>(MPI_SOLID_FLUID<VECTOR<float,3> >* const);
template class SOLIDS_FLUIDS_EXAMPLE<VECTOR<double,1> >;
template class SOLIDS_FLUIDS_EXAMPLE<VECTOR<double,2> >;
template class SOLIDS_FLUIDS_EXAMPLE<VECTOR<double,3> >;
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<double,1> >::Adjust_Output_Directory_For_MPI<MPI_SOLIDS<VECTOR<double,1> >*>(MPI_SOLIDS<VECTOR<double,1> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<double,2> >::Adjust_Output_Directory_For_MPI<MPI_SOLIDS<VECTOR<double,2> >*>(MPI_SOLIDS<VECTOR<double,2> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<double,3> >::Adjust_Output_Directory_For_MPI<MPI_SOLIDS<VECTOR<double,3> >*>(MPI_SOLIDS<VECTOR<double,3> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<double,1> >::Adjust_Output_Directory_For_MPI<MPI_UNIFORM_GRID<VECTOR<double,1> >*>(MPI_UNIFORM_GRID<VECTOR<double,1> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<double,2> >::Adjust_Output_Directory_For_MPI<MPI_UNIFORM_GRID<VECTOR<double,2> >*>(MPI_UNIFORM_GRID<VECTOR<double,2> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<double,3> >::Adjust_Output_Directory_For_MPI<MPI_UNIFORM_GRID<VECTOR<double,3> >*>(MPI_UNIFORM_GRID<VECTOR<double,3> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<double,1> >::Adjust_Output_Directory_For_MPI<MPI_SOLID_FLUID<VECTOR<double,1> >*>(MPI_SOLID_FLUID<VECTOR<double,1> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<double,2> >::Adjust_Output_Directory_For_MPI<MPI_SOLID_FLUID<VECTOR<double,2> >*>(MPI_SOLID_FLUID<VECTOR<double,2> >* const);
template void SOLIDS_FLUIDS_EXAMPLE<VECTOR<double,3> >::Adjust_Output_Directory_For_MPI<MPI_SOLID_FLUID<VECTOR<double,3> >*>(MPI_SOLID_FLUID<VECTOR<double,3> >* const);
}
