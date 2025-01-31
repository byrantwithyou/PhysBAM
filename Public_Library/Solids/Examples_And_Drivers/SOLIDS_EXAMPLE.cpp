//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_EXAMPLE
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/MATRIX_4X4.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Particles/PARTICLES.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MPI.h>
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
SOLIDS_EXAMPLE(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
    :BASE(stream_type,parse_args),use_melting(false),solids_parameters(*new SOLIDS_PARAMETERS<TV>),solid_body_collection(*new SOLID_BODY_COLLECTION<TV>),
    solids_evolution(new NEWMARK_EVOLUTION<TV>(solids_parameters,solid_body_collection,*this)),
    debug_particles(*new DEBUG_PARTICLES<TV>),opt_skip_debug_data(false)
{
    Set_Minimum_Collision_Thickness();

    bool opt_solidssymmqmr=false,opt_solidscr=false,opt_solidscg=false;
    parse_args.Add("-solidscfl",&solids_parameters.cfl,"cfl","solids CFL");
    parse_args.Add("-solidscg",&opt_solidscg,"Use CG for time integration");
    parse_args.Add("-solidscr",&opt_solidscr,"Use CONJUGATE_RESIDUAL for time integration");
    parse_args.Add("-solidssymmqmr",&opt_solidssymmqmr,"Use SYMMQMR for time integration");
    parse_args.Add("-rigidcfl",&solids_parameters.rigid_body_evolution_parameters.rigid_cfl,"cfl","rigid CFL");
    parse_args.Add("-skip_debug_data",&opt_skip_debug_data,"turn off file io for debug data");

    if(opt_solidscg) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cg;
    if(opt_solidscr) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cr;
    if(opt_solidssymmqmr) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_symmqmr;
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
Read_Output_Files_Solids()
{
    solid_body_collection.Read(viewer_dir,solids_parameters.write_static_variables_every_frame);
    Read_From_File(viewer_dir.current_directory+"/triangle_collision_parameters",
        solids_parameters.triangle_collision_parameters.repulsion_pair_update_count,
        solids_parameters.triangle_collision_parameters.topological_hierarchy_build_count);
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
// Function Adjust_Output_Directory_For_MPI
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Adjust_Output_Directory_For_MPI(const MPI_SOLIDS<TV>* mpi)
{
    if(mpi && mpi->Number_Of_Processors()>1)
        viewer_dir.output_directory+=LOG::sprintf("/%d",(mpi->rank+1));
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
Write_Output_Files() const
{
    Create_Directory(viewer_dir.output_directory);
    if(this->use_test_output){
        std::string file=LOG::sprintf("%s/%s-%03d.txt",viewer_dir.output_directory.c_str(),this->test_output_prefix.c_str(),viewer_dir.frame_stack(0));
        OCTAVE_OUTPUT<T> oo(file.c_str());
        if(solid_body_collection.deformable_body_collection.particles.X.m){
            oo.Write("db_X",solid_body_collection.deformable_body_collection.particles.X.Flattened());
            oo.Write("db_V",solid_body_collection.deformable_body_collection.particles.V.Flattened());}
        if(solid_body_collection.rigid_body_collection.rigid_body_particles.frame.m){
            RIGID_BODY_PARTICLES<TV>& particle=solid_body_collection.rigid_body_collection.rigid_body_particles;
            ARRAY_VIEW<T> f((T*)particle.frame.Get_Array_Pointer(),particle.frame.m*(sizeof(FRAME<TV>)/sizeof(T)));
            oo.Write("rb_frame",f);
            ARRAY_VIEW<T> t((T*)particle.twist.Get_Array_Pointer(),particle.twist.m*TWIST<TV>::m);
            oo.Write("rb_twist",t);}}

    debug_particles.Write_Debug_Particles(stream_type,viewer_dir);
    solid_body_collection.Write(stream_type,viewer_dir,solids_parameters.write_static_variables_every_frame,
        solids_parameters.write_from_every_process);
    if(NEWMARK_EVOLUTION<TV>* newmark=dynamic_cast<NEWMARK_EVOLUTION<TV>*>(solids_evolution))
        newmark->Write_Position_Update_Projection_Data(stream_type,viewer_dir.current_directory+"/");

    Write_To_File(stream_type,viewer_dir.current_directory+"/triangle_collision_parameters",
        solids_parameters.triangle_collision_parameters.repulsion_pair_update_count,
        solids_parameters.triangle_collision_parameters.topological_hierarchy_build_count);
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
