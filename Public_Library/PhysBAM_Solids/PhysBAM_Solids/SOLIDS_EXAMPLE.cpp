//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_EXAMPLE
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/SOLIDS_EXAMPLE.h>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_EXAMPLE<TV>::
SOLIDS_EXAMPLE(const STREAM_TYPE stream_type)
    :BASE(stream_type),solids_parameters(*new SOLIDS_PARAMETERS<TV>),solid_body_collection(*new SOLID_BODY_COLLECTION<TV>(this)),
    solids_evolution(new NEWMARK_EVOLUTION<TV>(solids_parameters,solid_body_collection)),opt_solidscg(false),opt_solidscr(false),opt_solidssymmqmr(false)
{
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
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Read_Output_Files_Solids(const int frame)
{
    solid_body_collection.Read(stream_type,output_directory,frame,frame,solids_parameters.write_static_variables_every_frame,solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,
        solids_parameters.write_deformable_body,solids_parameters.write_from_every_process);
}
//#####################################################################
// Function Log_Parameters
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Log_Parameters() const
{
    LOG::SCOPE scope("SOLIDS_EXAMPLE parameters");
    BASE::Log_Parameters();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Write_Output_Files(const int frame) const
{
    if(this->use_test_output){
        std::string file=STRING_UTILITIES::string_sprintf("%s/%s-%03d.txt",output_directory.c_str(),this->test_output_prefix.c_str(),frame);
        OCTAVE_OUTPUT<T> oo(file.c_str());
        if(solid_body_collection.deformable_body_collection.particles.X.m){
            oo.Write("db_X",solid_body_collection.deformable_body_collection.particles.X.Flattened());
            oo.Write("db_V",solid_body_collection.deformable_body_collection.particles.V.Flattened());}
        if(solid_body_collection.rigid_body_collection.rigid_body_particle.frame.m){
            RIGID_BODY_PARTICLES<TV>& particle=solid_body_collection.rigid_body_collection.rigid_body_particle;
            ARRAY_VIEW<T> f(particle.frame.m*(sizeof(FRAME<TV>)/sizeof(T)),(T*)particle.frame.Get_Array_Pointer());
            oo.Write("rb_frame",f);
            ARRAY_VIEW<T> t(particle.twist.m*TWIST<TV>::m,(T*)particle.twist.Get_Array_Pointer());
            oo.Write("rb_twist",t);}}

    FILE_UTILITIES::Create_Directory(output_directory);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    Write_Frame_Title(frame);
    solid_body_collection.Write(stream_type,output_directory,frame,first_frame,solids_parameters.write_static_variables_every_frame,
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,solids_parameters.write_deformable_body,solids_parameters.write_from_every_process,
        solids_parameters.triangle_collision_parameters.output_interaction_pairs);
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Register_Options()
{
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
template<class TV> void SOLIDS_EXAMPLE<TV>::
Parse_Late_Options()
{
    BASE::Parse_Late_Options();
    if(opt_solidscg) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cg;
    if(opt_solidscr) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cr;
    if(opt_solidssymmqmr) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_symmqmr;
}
//#####################################################################
template class SOLIDS_EXAMPLE<VECTOR<float,1> >;
template class SOLIDS_EXAMPLE<VECTOR<float,2> >;
template class SOLIDS_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLIDS_EXAMPLE<VECTOR<double,1> >;
template class SOLIDS_EXAMPLE<VECTOR<double,2> >;
template class SOLIDS_EXAMPLE<VECTOR<double,3> >;
#endif
