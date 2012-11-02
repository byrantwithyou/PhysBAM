//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Fluids/Coupled_Evolution/COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/FLUID_COLLISION_BODY_INACCURATE_UNION.h>
#include <PhysBAM_Dynamics/Coupled_Driver/PLS_FSI_EXAMPLE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
#include <PhysBAM_Dynamics/Forces_And_Torques/EULER_FLUID_FORCES.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_PARAMETERS.h>
#include <stdexcept>
using namespace PhysBAM;
namespace PhysBAM{template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);}
//#####################################################################
// Constructor
//#####################################################################
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4355) // 'this' : used in base member initializer list
#endif
template<class TV> PLS_FSI_EXAMPLE<TV>::
PLS_FSI_EXAMPLE(const STREAM_TYPE stream_type,const int number_of_regions)
    :BASE(stream_type),solids_parameters(*new SOLIDS_PARAMETERS<TV>),solids_fluids_parameters(*new SOLIDS_FLUIDS_PARAMETERS<TV>(this)),
    solid_body_collection(*new SOLID_BODY_COLLECTION<TV>(this)),solids_evolution(new NEWMARK_EVOLUTION<TV>(solids_parameters,solid_body_collection)),
    fluids_parameters(number_of_regions,FLUIDS_PARAMETERS<GRID<TV> >::WATER),fluid_collection(*fluids_parameters.grid),resolution(8),convection_order(1),
    use_pls_evolution_for_structure(false),two_phase(false),use_kang(false),print_matrix(false),test_system(false),kang_poisson_viscosity(0),m(1),s(1),kg(1),
    opt_skip_debug_data(false),opt_solidscg(false),opt_solidscr(false),opt_solidssymmqmr(false)
{
    Set_Minimum_Collision_Thickness();
    Set_Write_Substeps_Level(-1);
}
#ifdef _MSC_VER
#pragma warning(pop)
#endif
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PLS_FSI_EXAMPLE<TV>::
~PLS_FSI_EXAMPLE()
{
    fluids_parameters.projection=0;
    delete solids_evolution;
    delete &solid_body_collection;
    delete &solids_parameters;
    delete &solids_fluids_parameters;
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Register_Options()
{
    if(!parse_args) return;
    BASE::Register_Options();
    parse_args->Add("-solidscfl",&solids_parameters.cfl,"cfl","solids CFL");
    parse_args->Add("-solidscg",&opt_solidscg,"Use CG for time integration");
    parse_args->Add("-solidscr",&opt_solidscr,"Use CONJUGATE_RESIDUAL for time integration");
    parse_args->Add("-solidssymmqmr",&opt_solidssymmqmr,"Use SYMMQMR for time integration");
    parse_args->Add("-rigidcfl",&solids_parameters.rigid_body_evolution_parameters.rigid_cfl,"cfl","rigid CFL");
    parse_args->Add("-skip_debug_data",&opt_skip_debug_data,"turn off file io for debug data");
    parse_args->Add("-resolution",&resolution,"resolution","simulation resolution");
    parse_args->Add("-m",&m,"scale","length unit");
    parse_args->Add("-s",&s,"scale","time unit");
    parse_args->Add("-kg",&kg,"scale","mass unit");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Parse_Options()
{
    BASE::Parse_Options();
    fluids_parameters.write_debug_data=!opt_skip_debug_data;
}
//#####################################################################
// Function Add_Volumetric_Body_To_Fluid_Simulation
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,bool add_collision,bool add_coupling)
{
    RIGID_COLLISION_GEOMETRY<TV>* collision_geometry=new RIGID_COLLISION_GEOMETRY<TV>(rigid_body);
    if(add_collision) fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(collision_geometry,rigid_body.particle_index,true);
    if(add_coupling)
        if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_evolution))
            coupled_evolution->iterator_info.coupling_bodies.Append(collision_geometry);
    if(rigid_body.simplicial_object) rigid_body.simplicial_object->Initialize_Hierarchy();
}
//#####################################################################
// Function Add_Thin_Shell_To_Fluid_Simulation
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,bool add_collision,bool add_coupling)
{
    RIGID_COLLISION_GEOMETRY<TV>* collision_geometry=new RIGID_COLLISION_GEOMETRY<TV>(rigid_body);
    if(add_collision) fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(collision_geometry,rigid_body.particle_index,true);
    if(add_coupling)
        if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_evolution))
            coupled_evolution->iterator_info.coupling_bodies.Append(collision_geometry);
    rigid_body.thin_shell=true;
    if(collision_geometry->collision_thickness<minimum_collision_thickness) collision_geometry->collision_thickness=minimum_collision_thickness;
    if(rigid_body.simplicial_object) rigid_body.simplicial_object->Initialize_Hierarchy();
}
//#####################################################################
// Function Add_To_Fluid_Simulation
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions,bool add_collision,bool add_coupling)
{
    if(add_collision) fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&deformable_collisions,0,false);
    if(add_coupling)
        if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_evolution))
            coupled_evolution->iterator_info.coupling_bodies.Append(&deformable_collisions);
    if(deformable_collisions.collision_thickness<minimum_collision_thickness) deformable_collisions.collision_thickness=minimum_collision_thickness;
    deformable_collisions.Initialize_For_Thin_Shells_Fluid_Coupling();
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class TV> template<class GEOMETRY> void PLS_FSI_EXAMPLE<TV>::
Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity)
{
    for(UNIFORM_GRID_ITERATOR_FACE<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
        int axis=iterator.Axis();fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(axis)(iterator.Face_Index())=true;
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(axis)(iterator.Face_Index())=constant_source_velocity[axis];}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class TV> template<class GEOMETRY> void PLS_FSI_EXAMPLE<TV>::
Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity,const ARRAY<bool,FACE_INDEX<TV::m> >& invalid_mask)
{
    for(UNIFORM_GRID_ITERATOR_FACE<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
        if(!invalid_mask(axis,face_index) && source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
            fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(axis)(face_index)=true;
            fluid_collection.incompressible_fluid_collection.face_velocities.Component(axis)(face_index)=constant_source_velocity[axis];}}
}
//#####################################################################
// Function Adjust_Phi_With_Source
//#####################################################################
template<class TV> template<class GEOMETRY> void PLS_FSI_EXAMPLE<TV>::
Adjust_Phi_With_Source(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source)
{
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){TV source_X=world_to_source.Homogeneous_Times(iterator.Location());
        if(source.Lazy_Inside(source_X)) 
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=min(fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index()),
                                                                                          source.Signed_Distance(source_X));}
}
//#####################################################################
// Function Adjust_Phi_With_Source
//#####################################################################
template<class TV> template<class GEOMETRY> void PLS_FSI_EXAMPLE<TV>::
Adjust_Phi_With_Source(const GEOMETRY& source,const int region,const T_TRANSFORMATION_MATRIX& world_to_source)
{
    T bandwidth=3*fluids_parameters.grid->Minimum_Edge_Length();
    ARRAY<ARRAY<T,TV_INT> >& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        TV source_X=world_to_source.Homogeneous_Times(iterator.Location());
        if(source.Inside(source_X,-bandwidth)){
            T source_signed_distance=source.Signed_Distance(source_X);
            for(int i=0;i<fluids_parameters.number_of_regions;i++){
                if(i==region) phis(i)(iterator.Cell_Index())=min(phis(i)(iterator.Cell_Index()),source_signed_distance);
                else phis(i)(iterator.Cell_Index())=max(phis(i)(iterator.Cell_Index()),-source_signed_distance);}}}
}
//#####################################################################
// Function Revalidate_Fluid_Scalars
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Revalidate_Fluid_Scalars()
{
    LEVELSET<TV>& levelset=fluids_parameters.particle_levelset_evolution->Levelset(0);
    FAST_LEVELSET_ADVECTION<GRID<TV> >& levelset_advection=fluids_parameters.particle_levelset_evolution->Levelset_Advection(0);
    if(levelset_advection.nested_semi_lagrangian_collidable)
        levelset_advection.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(*fluids_parameters.grid,fluids_parameters.collidable_phi_replacement_value,levelset.phi);
}
//#####################################################################
// Function Revalidate_Phi_After_Modify_Levelset
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Revalidate_Phi_After_Modify_Levelset()
{
    LEVELSET<TV>& levelset=fluids_parameters.particle_levelset_evolution->Levelset(0);
    FAST_LEVELSET_ADVECTION<GRID<TV> >& levelset_advection=fluids_parameters.particle_levelset_evolution->Levelset_Advection(0);
    if(levelset_advection.nested_semi_lagrangian_collidable){
        levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_current=levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_next;
        levelset_advection.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(*fluids_parameters.grid,fluids_parameters.collidable_phi_replacement_value,levelset.phi);}
}
//#####################################################################
// Function Revalidate_Fluid_Velocity
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Revalidate_Fluid_Velocity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities)
{
    if(fluids_parameters.incompressible->nested_semi_lagrangian_collidable) 
        fluids_parameters.incompressible->nested_semi_lagrangian_collidable->Average_To_Invalidated_Face(*fluids_parameters.grid,face_velocities);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Get_Object_Velocities(LAPLACE_UNIFORM<GRID<TV> >* elliptic_solver,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
{
    if(fluids_parameters.solid_affects_fluid && (fluids_parameters.fluid_affects_solid || fluids_parameters.use_slip)){
        if(!fluids_parameters.use_slip){
           SOLID_FLUID_COUPLED_EVOLUTION<TV>& coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION<TV>&>(*solids_evolution);
           coupled_evolution.Apply_Solid_Boundary_Conditions(time,false,face_velocities);}
        else{
            SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>& coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution);
            coupled_evolution.Get_Coupled_Faces_And_Interpolated_Solid_Velocities(time,elliptic_solver->psi_N,face_velocities);}}
    else fluids_parameters.collision_bodies_affecting_fluid->Compute_Psi_N(elliptic_solver->psi_N,&face_velocities);
}
//#####################################################################
// Function Get_Levelset_Velocity
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET_MULTIPLE<GRID<TV> >& levelset_multiple,ARRAY<T,FACE_INDEX<TV::m> >& V_levelset,const T time) const
{
    ARRAY<T,FACE_INDEX<TV::m> >::Put(fluid_collection.incompressible_fluid_collection.face_velocities,V_levelset);
}
//#####################################################################
// Function Initialize_Swept_Occupied_Blocks_For_Advection
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Initialize_Swept_Occupied_Blocks_For_Advection(const T dt,const T time,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities)
{
    GRID<TV>& grid=*fluids_parameters.grid;
    T maximum_fluid_speed=face_velocities.Max_Abs().Max(),maximum_particle_speed=0;
    PARTICLE_LEVELSET_UNIFORM<GRID<TV> >& particle_levelset=fluids_parameters.particle_levelset_evolution->Particle_Levelset(0);
    if(particle_levelset.use_removed_negative_particles) for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(particle_levelset.levelset.grid);iterator.Valid();iterator.Next()){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=particle_levelset.removed_negative_particles(iterator.Cell_Index());
            if(particles) maximum_particle_speed=max(maximum_particle_speed,particles->V.Maximum_Magnitude());}
    if(particle_levelset.use_removed_positive_particles) for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(particle_levelset.levelset.grid);iterator.Valid();iterator.Next()){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=particle_levelset.removed_positive_particles(iterator.Cell_Index());
            if(particles) maximum_particle_speed=max(maximum_particle_speed,particles->V.Maximum_Magnitude());}
    T max_particle_collision_distance=0;
    max_particle_collision_distance=max(max_particle_collision_distance,fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).max_collision_distance_factor*grid.dX.Max());
    fluids_parameters.collision_bodies_affecting_fluid->Compute_Occupied_Blocks(true,
        dt*max(maximum_fluid_speed,maximum_particle_speed)+2*max_particle_collision_distance+(T).5*fluids_parameters.p_grid.dX.Max(),10);
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{}
//#####################################################################
// Function Read_Output_Files_Fluids
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Read_Output_Files_Fluids(const int frame)
{
    fluids_parameters.Read_Output_Files(stream_type,output_directory,frame);
    fluid_collection.Read_Output_Files(stream_type,output_directory,frame);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
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
            oo.Write("rb_twist",t);}
        if(fluids_parameters.incompressible){
            const ARRAY<T,FACE_INDEX<TV::m> >& u=fluid_collection.incompressible_fluid_collection.face_velocities;
            ARRAY_VIEW<T> a(u.buffer_size,u.base_pointer);
            oo.Write("if_u",a);}
        if(fluids_parameters.euler) oo.Write("cf_U",fluids_parameters.euler->U.array.Flattened());}

    FILE_UTILITIES::Create_Directory(output_directory);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    Write_Frame_Title(frame);
    solid_body_collection.Write(stream_type,output_directory,frame,first_frame,solids_parameters.write_static_variables_every_frame,
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,solids_parameters.write_deformable_body,solids_parameters.write_from_every_process,
        solids_parameters.triangle_collision_parameters.output_interaction_pairs);
    if(NEWMARK_EVOLUTION<TV>* newmark=dynamic_cast<NEWMARK_EVOLUTION<TV>*>(solids_evolution))
        newmark->Write_Position_Update_Projection_Data(stream_type,output_directory+"/"+f+"/");
    
    fluids_parameters.Write_Output_Files(stream_type,output_directory,first_frame,frame);
    fluid_collection.Write_Output_Files(stream_type,output_directory,frame);
}
//#####################################################################
// Function Log_Parameters 
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Log_Parameters() const
{       
    LOG::SCOPE scope("PLS_FSI_EXAMPLE parameters");
    BASE::Log_Parameters();
    LOG::cout<<"minimum_collision_thickness="<<minimum_collision_thickness<<std::endl;
    fluids_parameters.Log_Parameters();
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Read_Output_Files_Solids(const int frame)
{
    solid_body_collection.Read(stream_type,output_directory,frame,frame,solids_parameters.write_static_variables_every_frame,solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,
        solids_parameters.write_deformable_body,solids_parameters.write_from_every_process);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    //if(NEWMARK_EVOLUTION<TV>* newmark=dynamic_cast<NEWMARK_EVOLUTION<TV>*>(solids_evolution))
    //    newmark->Read_Position_Update_Projection_Data(stream_type,output_directory+"/"+f+"/");
}
//#####################################################################
// Function Parse_Late_Options
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Parse_Late_Options()
{
    if(!parse_args) return;
    BASE::Parse_Late_Options();
    if(opt_solidscg) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cg;
    if(opt_solidscr) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cr;
    if(opt_solidssymmqmr) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_symmqmr;
}
//#####################################################################
// Function Adjust_Particle_For_Domain_Boundaries
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    // remove this for speed - don't call the function for the other particles
    if(particle_type==PARTICLE_LEVELSET_POSITIVE || particle_type==PARTICLE_LEVELSET_REMOVED_POSITIVE) return;

    TV& X=particles.X(index);TV X_new=X+dt*V;
    T max_collision_distance=fluids_parameters.particle_levelset_evolution->particle_levelset.Particle_Collision_Distance(particles.quantized_collision_distance(index));
    T min_collision_distance=fluids_parameters.particle_levelset_evolution->particle_levelset.min_collision_distance_factor*max_collision_distance;
    TV min_corner=fluids_parameters.grid->domain.Minimum_Corner(),max_corner=fluids_parameters.grid->domain.Maximum_Corner();
    for(int axis=0;axis<TV::m;axis++){
        if(fluids_parameters.domain_walls[axis][0] && X_new[axis]<min_corner[axis]+max_collision_distance){
            T collision_distance=X[axis]-min_corner[axis];
            if(collision_distance>max_collision_distance)collision_distance=X_new[axis]-min_corner[axis];
            collision_distance=max(min_collision_distance,collision_distance);
            X_new[axis]+=max((T)0,min_corner[axis]-X_new[axis]+collision_distance);
            V[axis]=max((T)0,V[axis]);X=X_new-dt*V;}
        if(fluids_parameters.domain_walls[axis][1] && X_new[axis]>max_corner[axis]-max_collision_distance){
            T collision_distance=max_corner[axis]-X[axis];
            if(collision_distance>max_collision_distance) collision_distance=max_corner[axis]-X_new[axis];
            collision_distance=max(min_collision_distance,collision_distance);
            X_new[axis]-=max((T)0,X_new[axis]-max_corner[axis]+collision_distance);
            V[axis]=min((T)0,V[axis]);X=X_new-dt*V;}}
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Update_Fluid_Parameters(const T dt,const T time)
{
    fluids_parameters.incompressible->Set_Gravity(fluids_parameters.gravity,fluids_parameters.gravity_direction);
    fluids_parameters.incompressible->Set_Body_Force(fluids_parameters.use_body_force);
    fluids_parameters.incompressible->projection.Use_Non_Zero_Divergence(fluids_parameters.use_non_zero_divergence);
    fluids_parameters.incompressible->projection.elliptic_solver->Solve_Neumann_Regions(fluids_parameters.solve_neumann_regions);
    fluids_parameters.incompressible->projection.elliptic_solver->solve_single_cell_neumann_regions=fluids_parameters.solve_single_cell_neumann_regions;
    fluids_parameters.incompressible->Use_Explicit_Part_Of_Implicit_Viscosity(fluids_parameters.use_explicit_part_of_implicit_viscosity);
    if(fluids_parameters.implicit_viscosity && fluids_parameters.implicit_viscosity_iterations)
        fluids_parameters.incompressible->Set_Maximum_Implicit_Viscosity_Iterations(fluids_parameters.implicit_viscosity_iterations);
    fluids_parameters.incompressible->Use_Variable_Vorticity_Confinement(fluids_parameters.use_variable_vorticity_confinement);
    fluids_parameters.incompressible->Set_Surface_Tension(fluids_parameters.surface_tension);
    fluids_parameters.incompressible->Set_Variable_Surface_Tension(fluids_parameters.variable_surface_tension);
    fluids_parameters.incompressible->Set_Viscosity(fluids_parameters.viscosity);
    fluids_parameters.incompressible->Set_Variable_Viscosity(fluids_parameters.variable_viscosity);
    fluids_parameters.incompressible->projection.Set_Density(fluids_parameters.density);
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void PLS_FSI_EXAMPLE<TV>::
Set_Boundary_Conditions(ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& psi_D_value,
    ARRAY<T,FACE_INDEX<TV::dimension> >& psi_N_value) const
{
    void Resize(const RANGE<TV_INT>& box,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T());

    psi_D.Resize(fluids_parameters.grid->Domain_Indices(3),false,false);
    psi_N.Resize(fluids_parameters.grid->Domain_Indices(3),false,false);
    psi_D_value.Resize(fluids_parameters.grid->Domain_Indices(3),false,false);
    psi_N_value.Resize(fluids_parameters.grid->Domain_Indices(3),false,false);
    psi_D.Fill(false);
    psi_N.Fill(false);
    psi_D_value.Fill(0);
    psi_N_value.Fill(0);

    GRID<TV>& grid=*fluids_parameters.grid;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,3,GRID<TV>::GHOST_REGION);it.Valid();it.Next()){
        psi_D(it.index)=true;
        /*Add_Debug_Particle(it.Location(), VECTOR<T,3>(1,0,0));*/}
    for(int d=0;d<TV::m;d++)
        for(int i=0;i<2;i++)
            if(fluids_parameters.domain_walls(d)(i))
                for(UNIFORM_GRID_ITERATOR_FACE<TV> it(grid,0,GRID<TV>::BOUNDARY_REGION,i+2*(d-1),-1);it.Valid();it.Next()){
                    psi_N(it.Full_Index())=true;
                    /*Add_Debug_Particle(it.Location(),VECTOR<T,3>(0,1,0));*/}

    Set_Boundary_Conditions_Callback(psi_D,psi_N,psi_D_value,psi_N_value);
}
//#####################################################################
template class PLS_FSI_EXAMPLE<VECTOR<float,1> >;
template class PLS_FSI_EXAMPLE<VECTOR<float,2> >;
template class PLS_FSI_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PLS_FSI_EXAMPLE<VECTOR<double,1> >;
template class PLS_FSI_EXAMPLE<VECTOR<double,2> >;
template class PLS_FSI_EXAMPLE<VECTOR<double,3> >;
#endif
