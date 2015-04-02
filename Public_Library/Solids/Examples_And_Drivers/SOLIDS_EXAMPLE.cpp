//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_EXAMPLE
//#####################################################################
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Log/SCOPE.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Matrices/MATRIX_4X4.h>
#include <Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Particles/PARTICLES.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_EXAMPLE<TV>::
SOLIDS_EXAMPLE(const STREAM_TYPE stream_type)
    :BASE(stream_type),use_melting(false),solids_parameters(*new SOLIDS_PARAMETERS<TV>),solid_body_collection(*new SOLID_BODY_COLLECTION<TV>),
    solids_evolution(new NEWMARK_EVOLUTION<TV>(solids_parameters,solid_body_collection,*this)),opt_solidssymmqmr(false),opt_solidscr(false),
    opt_solidscg(false),debug_particles(*new DEBUG_PARTICLES<TV>),opt_skip_debug_data(false)
{
    Set_Minimum_Collision_Thickness();
    Set_Write_Substeps_Level(-1);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_EXAMPLE<TV>::
~SOLIDS_EXAMPLE()
{
    delete solids_evolution;
    delete &solid_body_collection;
    delete &solids_parameters;
    delete &debug_particles;
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
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
template<class TV> void SOLIDS_EXAMPLE<TV>::
Log_Parameters() const
{
    LOG::SCOPE scope("SOLIDS_EXAMPLE parameters");
    BASE::Log_Parameters();
    LOG::cout<<"minimum_collision_thickness="<<minimum_collision_thickness<<std::endl;
    LOG::cout<<"use_melting="<<use_melting<<std::endl;
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Register_Options()
{
    if(!parse_args) return;
    BASE::Register_Options();
    parse_args->Add("-solidscfl",&solids_parameters.cfl,"cfl","solids CFL");
    parse_args->Add("-solidscg",&opt_solidscg,"Use CG for time integration");
    parse_args->Add("-solidscr",&opt_solidscr,"Use CONJUGATE_RESIDUAL for time integration");
    parse_args->Add("-solidssymmqmr",&opt_solidssymmqmr,"Use SYMMQMR for time integration");
    parse_args->Add("-rigidcfl",&solids_parameters.rigid_body_evolution_parameters.rigid_cfl,"cfl","rigid CFL");
    BASE::Register_Options();
    parse_args->Add("-skip_debug_data",&opt_skip_debug_data,"turn off file io for debug data");
}
//#####################################################################
// Function Parse_Late_Options
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
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
template<class TV> void SOLIDS_EXAMPLE<TV>::
Adjust_Output_Directory_For_MPI(const MPI_SOLIDS<TV>* mpi)
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
template<class TV> void SOLIDS_EXAMPLE<TV>::
Post_Initialization()
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Preprocess_Frame(const int frame)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Postprocess_Frame(const int frame)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Preprocess_Substep(const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Postprocess_Substep(const T dt,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Initialize_Bodies() 
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Write_Output_Files(const int frame) const
{
    if(this->use_test_output){
        std::string file=LOG::sprintf("%s/%s-%03d.txt",output_directory.c_str(),this->test_output_prefix.c_str(),frame);
        OCTAVE_OUTPUT<T> oo(file.c_str());
        if(solid_body_collection.deformable_body_collection.particles.X.m){
            oo.Write("db_X",solid_body_collection.deformable_body_collection.particles.X.Flattened());
            oo.Write("db_V",solid_body_collection.deformable_body_collection.particles.V.Flattened());}
        if(solid_body_collection.rigid_body_collection.rigid_body_particles.frame.m){
            RIGID_BODY_PARTICLES<TV>& particle=solid_body_collection.rigid_body_collection.rigid_body_particles;
            ARRAY_VIEW<T> f(particle.frame.m*(sizeof(FRAME<TV>)/sizeof(T)),(T*)particle.frame.Get_Array_Pointer());
            oo.Write("rb_frame",f);
            ARRAY_VIEW<T> t(particle.twist.m*TWIST<TV>::m,(T*)particle.twist.Get_Array_Pointer());
            oo.Write("rb_twist",t);}}

    FILE_UTILITIES::Create_Directory(output_directory);
    std::string f=LOG::sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    Write_Frame_Title(frame);
    debug_particles.Write_Debug_Particles(stream_type,output_directory,frame);
    solid_body_collection.Write(stream_type,output_directory,frame,first_frame,solids_parameters.write_static_variables_every_frame,
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,solids_parameters.write_deformable_body,solids_parameters.write_from_every_process,
        solids_parameters.triangle_collision_parameters.output_interaction_pairs);
    if(NEWMARK_EVOLUTION<TV>* newmark=dynamic_cast<NEWMARK_EVOLUTION<TV>*>(solids_evolution))
        newmark->Write_Position_Update_Projection_Data(stream_type,output_directory+"/"+f+"/");
}
//#####################################################################
namespace PhysBAM{
template class SOLIDS_EXAMPLE<VECTOR<float,1> >;
template class SOLIDS_EXAMPLE<VECTOR<float,2> >;
template class SOLIDS_EXAMPLE<VECTOR<float,3> >;
template class SOLIDS_EXAMPLE<VECTOR<double,1> >;
template class SOLIDS_EXAMPLE<VECTOR<double,2> >;
template class SOLIDS_EXAMPLE<VECTOR<double,3> >;
}
