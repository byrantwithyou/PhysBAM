//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Particles/PARTICLES.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Fluids/Coupled_Evolution/COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/FLUID_COLLISION_BODY_INACCURATE_UNION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
#include <PhysBAM_Dynamics/Forces_And_Torques/EULER_FLUID_FORCES.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <PhysBAM_Dynamics/Interpolation/FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
SOLIDS_FLUIDS_EXAMPLE_UNIFORM(const STREAM_TYPE stream_type,const int number_of_regions,const typename FLUIDS_PARAMETERS<T_GRID>::TYPE type)
    :SOLIDS_FLUIDS_EXAMPLE<TV>(stream_type),fluids_parameters(number_of_regions,type),fluid_collection(*fluids_parameters.grid),resolution(8),
    debug_particles(*new DEBUG_PARTICLES<TV>),opt_skip_debug_data(false)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
~SOLIDS_FLUIDS_EXAMPLE_UNIFORM()
{
    if(dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_evolution)) fluids_parameters.projection=0;
    delete &debug_particles;
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Register_Options()
{
    BASE::Register_Options();
    this->parse_args->Add("-skip_debug_data",&opt_skip_debug_data,"turn off file io for debug data");
    this->parse_args->Add("-use_fmm_extrapolation",&fluids_parameters.euler_solid_fluid_coupling_utilities->use_fast_marching,
        "use fast marching to extrapolate compressible flow data into solid state.");
    this->parse_args->Add("-resolution",&resolution,"resolution","simulation resolution");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Parse_Options()
{
    BASE::Parse_Options();
    fluids_parameters.write_debug_data=!opt_skip_debug_data;
}
//#####################################################################
// Function Add_Volumetric_Body_To_Fluid_Simulation
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
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
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
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
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
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
// Function Initialize_MPI
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Initialize_MPI()
{
    fluids_parameters.mpi_grid->Initialize(fluids_parameters.domain_walls);
    if(fluids_parameters.number_of_regions==1)fluids_parameters.phi_boundary=new BOUNDARY_MPI<TV>(fluids_parameters.mpi_grid,*fluids_parameters.phi_boundary);
    else if(fluids_parameters.number_of_regions>=2)for(int i=0;i<fluids_parameters.number_of_regions;i++)
        fluids_parameters.phi_boundary_multiphase(i)=new BOUNDARY_MPI<TV>(fluids_parameters.mpi_grid,*fluids_parameters.phi_boundary_multiphase(i));
    fluids_parameters.fluid_boundary=new BOUNDARY_MPI<TV>(fluids_parameters.mpi_grid,*fluids_parameters.fluid_boundary);
    if(fluids_parameters.use_soot){
        if(fluids_parameters.soot_boundary) fluids_parameters.soot_boundary=new BOUNDARY_MPI<TV>(fluids_parameters.mpi_grid,*fluids_parameters.soot_boundary);
        fluids_parameters.soot_container.Set_Custom_Boundary(*new BOUNDARY_MPI<TV>(fluids_parameters.mpi_grid,*fluids_parameters.soot_container.boundary));
        fluids_parameters.soot_fuel_container.Set_Custom_Boundary(*new BOUNDARY_MPI<TV>(fluids_parameters.mpi_grid,*fluids_parameters.soot_fuel_container.boundary));}
    if(fluids_parameters.use_density){
        if(fluids_parameters.density_boundary) fluids_parameters.density_boundary=new BOUNDARY_MPI<TV>(fluids_parameters.mpi_grid,*fluids_parameters.density_boundary);
        fluids_parameters.density_container.Set_Custom_Boundary(*new BOUNDARY_MPI<TV>(fluids_parameters.mpi_grid,*fluids_parameters.density_container.boundary));}
    if(fluids_parameters.use_temperature){
        if(fluids_parameters.temperature_boundary) fluids_parameters.temperature_boundary=new BOUNDARY_MPI<TV>(fluids_parameters.mpi_grid,*fluids_parameters.temperature_boundary);
        fluids_parameters.temperature_container.Set_Custom_Boundary(*new BOUNDARY_MPI<TV>(fluids_parameters.mpi_grid,*fluids_parameters.temperature_container.boundary));}
    for(int i=0;i<fluids_parameters.number_of_regions;i++) fluids_parameters.particle_levelset_evolution->Particle_Levelset(i).mpi_grid=fluids_parameters.mpi_grid;
    if(fluids_parameters.use_strain||fluids_parameters.use_multiphase_strain.Count_Matches(0)<fluids_parameters.number_of_regions)
        fluids_parameters.strain_boundary=new BOUNDARY_MPI<TV,SYMMETRIC_MATRIX<T,TV::m> >(fluids_parameters.mpi_grid,*fluids_parameters.strain_boundary);
    if(fluids_parameters.incompressible){
        fluids_parameters.incompressible->mpi_grid=fluids_parameters.mpi_grid;
        fluids_parameters.incompressible->projection.elliptic_solver->mpi_grid=fluids_parameters.mpi_grid;}
    if(fluids_parameters.compressible){
        fluids_parameters.compressible_boundary=new BOUNDARY_MPI<TV,VECTOR<typename T_GRID::SCALAR,T_GRID::dimension+2> >(fluids_parameters.mpi_grid,*fluids_parameters.compressible_boundary);
        fluids_parameters.euler->euler_projection.pressure_boundary=new BOUNDARY_MPI<TV>(fluids_parameters.mpi_grid,*fluids_parameters.euler->euler_projection.pressure_boundary);
        fluids_parameters.euler->mpi_grid=fluids_parameters.mpi_grid;
        fluids_parameters.euler->euler_projection.elliptic_solver->mpi_grid=fluids_parameters.mpi_grid;}
    if(!restart && fluids_parameters.store_particle_ids){ // TODO: this should be fixed so that id's are scalable
        for(int i=0;i<fluids_parameters.number_of_regions;i++)
            fluids_parameters.particle_levelset_evolution->Particle_Levelset(i).last_unique_particle_id=fluids_parameters.mpi_grid->rank*30000000;}       
}
//#####################################################################
// Function Initialize_Solid_Fluid_Coupling
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization()
{
    fluids_parameters.use_poisson=true;
    if(fluids_parameters.use_slip){
        SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=new SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>(solids_parameters,solid_body_collection,*this,
            fluids_parameters,solids_fluids_parameters,fluid_collection);
        delete solids_evolution;
        solids_evolution=coupled_evolution;
        fluids_parameters.Set_Projection(coupled_evolution);}
    else{
        delete solids_evolution;
        solids_evolution=new SOLID_FLUID_COUPLED_EVOLUTION<TV>(solids_parameters,solid_body_collection,*this,fluids_parameters,fluid_collection,solids_fluids_parameters);}
    // TODO: set up anything that needs to be set up for solid_affects_fluid only case
}
//#####################################################################
// Function Initialize_Solid_Fluid_Coupling
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Initialize_Solid_Fluid_Coupling_After_Grid_Initialization()
{
    fluids_parameters.collision_bodies_affecting_fluid->Initialize_Grids();

#if 0
    bool initialize_valid_value_mask=!restart; // If we're restarting we'll be reading in the values anyway
    // TODO: update for new advections!
    thin_shells_semi_lagrangian_density.Initialize(fluids_parameters.grid,initialize_valid_value_mask);
    thin_shells_semi_lagrangian_density.interpolation->Set_Default_Replacement_Value(fluids_parameters.density_container.ambient_density);
    thin_shells_semi_lagrangian_temperature.Initialize(fluids_parameters.grid,initialize_valid_value_mask);
    thin_shells_semi_lagrangian_temperature.interpolation->Set_Default_Replacement_Value(fluids_parameters.temperature_container.ambient_temperature);
    thin_shells_semi_lagrangian_velocity.Initialize(fluids_parameters.grid,initialize_valid_value_mask);
    if(fluids_parameters.water||fluids_parameters.fire) thin_shells_semi_lagrangian_phi.Initialize(fluids_parameters.grid,initialize_valid_value_mask);
    if(fluids_parameters.water) thin_shells_semi_lagrangian_phi.interpolation->Set_Default_Replacement_Value((T)1e-5);
    if(fluids_parameters.fire){ // special treatment for fire because we advect different quantities using different velocity fields
        thin_shells_semi_lagrangian_phi.interpolation->Set_Default_Replacement_Value((T)1e-5);} // bias to negative phi for fire
#endif

    if(fluids_parameters.compressible){
        EULER_UNIFORM<T_GRID>& euler=*fluids_parameters.euler;
        SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>& euler_solid_fluid_coupling_utilities=*fluids_parameters.euler_solid_fluid_coupling_utilities;
        fluids_parameters.euler_solid_fluid_coupling_utilities->Initialize_Solid_Fluid_Coupling(fluids_parameters.collision_bodies_affecting_fluid);
        if(fluids_parameters.fluid_affects_solid && !euler.timesplit){
            EULER_FLUID_FORCES<T_GRID>* euler_fluid_forces=new EULER_FLUID_FORCES<T_GRID>(euler.grid,euler_solid_fluid_coupling_utilities.pressure_at_faces,euler_solid_fluid_coupling_utilities.solid_fluid_face_time_n,
                euler_solid_fluid_coupling_utilities.cells_inside_fluid_time_n,euler_solid_fluid_coupling_utilities.collision_bodies_affecting_fluid,
                solid_body_collection.deformable_body_collection.particles,solid_body_collection.rigid_body_collection);
            solid_body_collection.Add_Force(euler_fluid_forces);}}

    if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_evolution))
        coupled_evolution->Setup_Boundary_Condition_Collection();
}
//#####################################################################
// Function Initialize_Compressible_Incompressible_Coupling
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Initialize_Compressible_Incompressible_Coupling()
{
    if(fluids_parameters.compressible && fluids_parameters.number_of_regions){
        fluids_parameters.euler->euler_projection.Set_Incompressible_Coupling_Callbacks(fluids_parameters.compressible_incompressible_coupling_utilities);}
}
//#####################################################################
// Function Set_Ghost_Density_And_Temperature_Inside_Flame_Core
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Set_Ghost_Density_And_Temperature_Inside_Flame_Core()
{
    for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        int region=fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Inside_Region(iterator.Cell_Index());
        if(fluids_parameters.fuel_region(region)){
            fluids_parameters.density_container.density(iterator.Cell_Index())=fluids_parameters.density;
            fluids_parameters.temperature_container.temperature(iterator.Cell_Index())=fluids_parameters.temperature_products;}}
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Set_Dirichlet_Boundary_Conditions(const T time)
{
    if(fluids_parameters.smoke||fluids_parameters.fire||fluids_parameters.water){
        if(fluids_parameters.number_of_regions>=2)
            fluids_parameters.incompressible_multiphase->Set_Dirichlet_Boundary_Conditions(fluids_parameters.particle_levelset_evolution_multiple->phis,fluids_parameters.dirichlet_regions);
        else if(fluids_parameters.number_of_regions==1) fluids_parameters.incompressible->Set_Dirichlet_Boundary_Conditions(&fluids_parameters.particle_levelset_evolution->phi,0);
        else fluids_parameters.incompressible->Set_Dirichlet_Boundary_Conditions(0,0);} // 1 phase
    else if(fluids_parameters.sph) for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        fluids_parameters.incompressible->projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;fluids_parameters.incompressible->projection.p(iterator.Cell_Index())=0;}
    else if(fluids_parameters.compressible){
        fluids_parameters.euler->euler_projection.Set_Dirichlet_Boundary_Conditions(time);
        if(fluids_parameters.number_of_regions) fluids_parameters.incompressible->Set_Dirichlet_Boundary_Conditions(&fluids_parameters.particle_levelset_evolution->phi,
            fluids_parameters.compressible_incompressible_coupling_utilities->p_dirichlet_incompressible);}

    if(fluids_parameters.solid_affects_fluid && (fluids_parameters.use_slip || fluids_parameters.fluid_affects_solid)){
        if(fluids_parameters.use_slip){
            // slip sets up its own boundary conditions
            //SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_parameters.solids_evolution);
            //coupled_evolution->Set_Dirichlet_Boundary_Conditions(fluid_collection.incompressible_fluid_collection.face_velocities,time);
        }
        else{
            SOLID_FLUID_COUPLED_EVOLUTION<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION<TV>*>(solids_evolution);
            coupled_evolution->Set_Dirichlet_Boundary_Conditions(time);}}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity)
{
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
        int axis=iterator.Axis();fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(axis)(iterator.Face_Index())=true;
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(axis)(iterator.Face_Index())=constant_source_velocity[axis];}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity,const T_FACE_ARRAYS_BOOL& invalid_mask)
{
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
        if(!invalid_mask(axis,face_index) && source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
            fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(axis)(face_index)=true;
            fluid_collection.incompressible_fluid_collection.face_velocities.Component(axis)(face_index)=constant_source_velocity[axis];}}
}
//#####################################################################
// Function Adjust_Phi_With_Source
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Adjust_Phi_With_Source(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source)
{
    for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){TV source_X=world_to_source.Homogeneous_Times(iterator.Location());
        if(source.Lazy_Inside(source_X)) 
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=min(fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index()),
                                                                                          source.Signed_Distance(source_X));}
}
//#####################################################################
// Function Adjust_Phi_With_Source
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Adjust_Phi_With_Source(const GEOMETRY& source,const int region,const T_TRANSFORMATION_MATRIX& world_to_source)
{
    T bandwidth=3*fluids_parameters.grid->dX.Min();
    ARRAY<T_ARRAYS_SCALAR>& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
    for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        TV source_X=world_to_source.Homogeneous_Times(iterator.Location());
        if(source.Inside(source_X,-bandwidth)){
            T source_signed_distance=source.Signed_Distance(source_X);
            for(int i=0;i<fluids_parameters.number_of_regions;i++){
                if(i==region) phis(i)(iterator.Cell_Index())=min(phis(i)(iterator.Cell_Index()),source_signed_distance);
                else phis(i)(iterator.Cell_Index())=max(phis(i)(iterator.Cell_Index()),-source_signed_distance);}}}
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Get_Source_Reseed_Mask(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,ARRAY<bool,TV_INT>*& cell_centered_mask,const bool reset_mask)
{
    if(reset_mask){if(cell_centered_mask) delete cell_centered_mask;cell_centered_mask=new ARRAY<bool,TV_INT>(fluids_parameters.grid->Domain_Indices(1));}
    T padding=3*fluids_parameters.grid->dX.Max();
    for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        if(!source.Outside(world_to_source.Homogeneous_Times(iterator.Location()),padding)) (*cell_centered_mask)(iterator.Cell_Index())=true;
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Adjust_Density_And_Temperature_With_Sources(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const T source_density,const T source_temperature)
{
    for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
        if(fluids_parameters.use_density) fluids_parameters.density_container.density(iterator.Cell_Index())=source_density;
        if(fluids_parameters.use_temperature) fluids_parameters.temperature_container.temperature(iterator.Cell_Index())=source_temperature;}
}
//#####################################################################
// Function Revalidate_Fluid_Scalars
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Revalidate_Fluid_Scalars()
{
    for(int i=0;i<fluids_parameters.number_of_regions;i++){
        LEVELSET<TV>& levelset=fluids_parameters.particle_levelset_evolution->Levelset(i);
        LEVELSET_ADVECTION<TV>& levelset_advection=fluids_parameters.particle_levelset_evolution->Levelset_Advection(i);
        int sign=1;if(fluids_parameters.number_of_regions>=2&&fluids_parameters.dirichlet_regions(i))sign=-1;
        if(levelset_advection.nested_semi_lagrangian_collidable)
            levelset_advection.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(*fluids_parameters.grid,sign*fluids_parameters.collidable_phi_replacement_value,levelset.phi);}
    if(fluids_parameters.use_density){
        if(fluids_parameters.density_container.nested_semi_lagrangian_collidable)
            fluids_parameters.density_container.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(fluids_parameters.density_container.grid,fluids_parameters.density_container.ambient_density,
            fluids_parameters.density_container.density);}
    if(fluids_parameters.use_temperature){
        if(fluids_parameters.temperature_container.nested_semi_lagrangian_collidable)
            fluids_parameters.temperature_container.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(fluids_parameters.temperature_container.grid,
                fluids_parameters.temperature_container.ambient_temperature,fluids_parameters.temperature_container.temperature);}
}
//#####################################################################
// Function Revalidate_Phi_After_Modify_Levelset
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Revalidate_Phi_After_Modify_Levelset()
{
    for(int i=0;i<fluids_parameters.number_of_regions;i++){
        LEVELSET<TV>& levelset=fluids_parameters.particle_levelset_evolution->Levelset(i);
        LEVELSET_ADVECTION<TV>& levelset_advection=fluids_parameters.particle_levelset_evolution->Levelset_Advection(i);
        int sign=1;if(fluids_parameters.number_of_regions>=2&&fluids_parameters.dirichlet_regions(i))sign=-1;
        if(levelset_advection.nested_semi_lagrangian_collidable){
            levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_current=levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_next;
            levelset_advection.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(*fluids_parameters.grid,sign*fluids_parameters.collidable_phi_replacement_value,levelset.phi);}}
}
//#####################################################################
// Function Revalidate_Fluid_Velocity
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Revalidate_Fluid_Velocity(T_FACE_ARRAYS_SCALAR& face_velocities)
{
    if(fluids_parameters.incompressible->nested_nested_semi_lagrangian_fire_multiphase_collidable) 
        fluids_parameters.incompressible->nested_nested_semi_lagrangian_fire_multiphase_collidable->Average_To_Invalidated_Face(*fluids_parameters.grid,face_velocities);
    if(fluids_parameters.incompressible->nested_semi_lagrangian_collidable) 
        fluids_parameters.incompressible->nested_semi_lagrangian_collidable->Average_To_Invalidated_Face(*fluids_parameters.grid,face_velocities);
    //if(fluids_parameters.incompressible->nested_semi_lagrangian_collidable_slip)
    //    fluids_parameters.incompressible->nested_semi_lagrangian_collidable_slip->Average_To_Invalidated_Face(*fluids_parameters.grid,face_velocities);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Get_Object_Velocities(LAPLACE_UNIFORM<T_GRID>* elliptic_solver,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
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
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Get_Levelset_Velocity(const T_GRID& grid,LEVELSET_MULTIPLE<T_GRID>& levelset_multiple,T_FACE_ARRAYS_SCALAR& V_levelset,const T time) const
{
    if(!fluids_parameters.use_reacting_flow) T_FACE_ARRAYS_SCALAR::Put(fluid_collection.incompressible_fluid_collection.face_velocities,V_levelset);
    else{
        const PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection=fluids_parameters.incompressible_multiphase->projection;
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
            TV_INT cell1,cell2;grid.Cells_Touching_Face(axis,face,cell1,cell2);
            int region1,other_region1,region2,other_region2;T phi1,other_phi1,phi2,other_phi2;
            levelset_multiple.Two_Minimum_Regions(cell1,region1,other_region1,phi1,other_phi1);
            levelset_multiple.Two_Minimum_Regions(cell2,region2,other_region2,phi2,other_phi2);
            // TODO: break out of this if phi1 and ph2 are far from the interface
            int face_region=region1;
            if(region1==region2){ // not close enough to the interface
                if(phi1>phi2){region2=other_region1;phi2=other_phi1;} // cell1 is close to the interface
                else{region1=other_region2;phi1=other_phi2;}} // cell2 is closer to the interface
            else if(phi1>phi2) face_region=region2; // see LEVELSET_MULTIPLE::Inside_Region_Face()
            if(projection.flame_speed_constants(region1,region2).z==0) V_levelset.Component(axis)(face)=fluid_collection.incompressible_fluid_collection.face_velocities.Component(axis)(face);
            else{
                int fuel_region=region1,product_region=region2;if(fluids_parameters.densities(fuel_region)<fluids_parameters.densities(product_region)){fuel_region=region2;product_region=region1;}
                V_levelset.Component(axis)(face)=projection.Face_Velocity_With_Ghost_Value_Multiphase(fluid_collection.incompressible_fluid_collection.face_velocities,axis,face,fuel_region,face_region)-
                    projection.Flame_Speed_Face_Multiphase(axis,face,fuel_region,product_region)*
                    (levelset_multiple.Phi(fuel_region,cell2)-levelset_multiple.Phi(fuel_region,cell1))*grid.one_over_dX[axis];}}}
    if(fluids_parameters.pseudo_dirichlet_regions.Number_True()>0){
        T_ARRAYS_SCALAR phi_for_pseudo_dirichlet_regions;T_GRID grid_temp(grid);LEVELSET<TV> levelset_for_pseudo_dirichlet_regions(grid_temp,phi_for_pseudo_dirichlet_regions);
        fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Get_Single_Levelset(fluids_parameters.pseudo_dirichlet_regions,levelset_for_pseudo_dirichlet_regions,false);
        fluids_parameters.incompressible->Extrapolate_Velocity_Across_Interface(V_levelset,phi_for_pseudo_dirichlet_regions);}
}
//#####################################################################
// Function Initialize_Swept_Occupied_Blocks_For_Advection
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Initialize_Swept_Occupied_Blocks_For_Advection(const T dt,const T time,T maximum_fluid_speed,const bool include_removed_particle_velocities)
{
    T_GRID& grid=*fluids_parameters.grid;
    T maximum_particle_speed=0;
    if(fluids_parameters.fire){
        if(fluids_parameters.number_of_regions==1){
            Get_Levelset_Velocity(*fluids_parameters.grid,fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).levelset,
                fluids_parameters.particle_levelset_evolution->V,time);
            maximum_fluid_speed=max(maximum_fluid_speed,fluids_parameters.particle_levelset_evolution->V.Max_Abs().Max());}
        else{
            Get_Levelset_Velocity(*fluids_parameters.grid,fluids_parameters.particle_levelset_evolution_multiple->Levelset_Multiple(),
                fluids_parameters.particle_levelset_evolution_multiple->V,time);
            maximum_fluid_speed=max(maximum_fluid_speed,fluids_parameters.particle_levelset_evolution_multiple->V.Max_Abs().Max());}}
    if(include_removed_particle_velocities){
        for(int i=0;i<fluids_parameters.number_of_regions;i++){
            PARTICLE_LEVELSET_UNIFORM<T_GRID>& particle_levelset=fluids_parameters.particle_levelset_evolution->Particle_Levelset(i);
            if(particle_levelset.use_removed_negative_particles) for(CELL_ITERATOR<TV> iterator(particle_levelset.levelset.grid);iterator.Valid();iterator.Next()){
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=particle_levelset.removed_negative_particles(iterator.Cell_Index());
                if(particles) maximum_particle_speed=max(maximum_particle_speed,particles->V.Maximum_Magnitude());}
            if(particle_levelset.use_removed_positive_particles) for(CELL_ITERATOR<TV> iterator(particle_levelset.levelset.grid);iterator.Valid();iterator.Next()){
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=particle_levelset.removed_positive_particles(iterator.Cell_Index());
                if(particles) maximum_particle_speed=max(maximum_particle_speed,particles->V.Maximum_Magnitude());}}}
    T max_particle_collision_distance=0;
    for(int i=0;i<fluids_parameters.number_of_regions;i++)
        max_particle_collision_distance=max(max_particle_collision_distance,fluids_parameters.particle_levelset_evolution->Particle_Levelset(i).max_collision_distance_factor*grid.dX.Max());
    fluids_parameters.collision_bodies_affecting_fluid->Compute_Occupied_Blocks(true,
        dt*max(maximum_fluid_speed,maximum_particle_speed)+2*max_particle_collision_distance+(T).5*fluids_parameters.p_grid.dX.Max(),10);

    if(fluids_parameters.use_maccormack_semi_lagrangian_advection && fluids_parameters.use_maccormack_compute_mask){
        typedef typename T_GRID::VECTOR_INT TV_INT;
        fluids_parameters.maccormack_cell_mask.Resize(grid.Domain_Indices(fluids_parameters.number_of_ghost_cells),false,false);
        fluids_parameters.maccormack_face_mask.Resize(grid,fluids_parameters.number_of_ghost_cells);
        // don't use maccormack near domain boundary conditions
        for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
            fluids_parameters.maccormack_cell_mask(iterator.Cell_Index())=TV_INT::Componentwise_Min(iterator.Cell_Index()-grid.Domain_Indices().Minimum_Corner(),
                grid.Domain_Indices().Maximum_Corner()-iterator.Cell_Index()).Min()>fluids_parameters.cfl;
        VECTOR<T_GRID,T_GRID::dimension> grids;
        for(int i=0;i<T_GRID::dimension;i++) grids[i]=grid.Get_Face_Grid(i);
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT index=iterator.Face_Index();
            fluids_parameters.maccormack_face_mask(axis,index)
                =TV_INT::Componentwise_Min(index-grids[axis].Domain_Indices().Minimum_Corner(),grids[axis].Domain_Indices().Maximum_Corner()-index).Min()>fluids_parameters.cfl;}
        if(fluids_parameters.mpi_grid){ // if mpi check turned off domain boundary regions to make sure they are global domain boundaries
            RANGE<TV> global_domain=fluids_parameters.mpi_grid->global_grid.domain;T double_cfl_dx=(T)fluids_parameters.cfl*fluids_parameters.p_grid.dX.Max();
            for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
                if(!fluids_parameters.maccormack_cell_mask(iterator.Cell_Index()) && global_domain.Inside(iterator.Location(),(T)double_cfl_dx))
                    fluids_parameters.maccormack_cell_mask(iterator.Cell_Index())=true;
            for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
                int axis=iterator.Axis();TV_INT index=iterator.Face_Index();
                if(!fluids_parameters.maccormack_face_mask(axis,index) && global_domain.Inside(iterator.Location(),(T)double_cfl_dx))
                    fluids_parameters.maccormack_face_mask(axis,index)=true;}}
        // don't use maccormack near the interface
        if(fluids_parameters.bandwidth_without_maccormack_near_interface){
            T dx_times_band=grid.dX.Max()*fluids_parameters.bandwidth_without_maccormack_near_interface;
            for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())if(fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())>-dx_times_band){
                fluids_parameters.maccormack_cell_mask(iterator.Cell_Index())=false;
                for(int axis=0;axis<T_GRID::dimension;axis++){
                    fluids_parameters.maccormack_face_mask(axis,iterator.First_Face_Index(axis))=false;fluids_parameters.maccormack_face_mask(axis,iterator.Second_Face_Index(axis))=false;}}}
        // turn off maccormack near objects
        for(NODE_ITERATOR<TV> node_iterator(grid,2);node_iterator.Valid();node_iterator.Next()) 
            if(fluids_parameters.collision_bodies_affecting_fluid->occupied_blocks(node_iterator.Node_Index())){
                TV_INT block_index=node_iterator.Node_Index();BLOCK_UNIFORM<T_GRID> block(grid,block_index);
                for(int cell_index=0;cell_index<T_GRID::number_of_cells_per_block;cell_index++) fluids_parameters.maccormack_cell_mask(block.Cell(cell_index))=false;
                for(int axis=0;axis<T_GRID::dimension;axis++) for(int face=0;face<T_GRID::number_of_incident_faces_per_block;face++)
                    fluids_parameters.maccormack_face_mask(axis,block.Incident_Face(axis,face))=false;}}
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{}
//#####################################################################
// Function Read_Output_Files_Fluids
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Read_Output_Files_Fluids(const int frame)
{
    fluids_parameters.Read_Output_Files(stream_type,output_directory,frame);
    fluid_collection.Read_Output_Files(stream_type,output_directory,frame);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    if(fluids_parameters.smoke||fluids_parameters.fire||fluids_parameters.water){
        if(fluids_parameters.solid_affects_fluid && fluids_parameters.fluid_affects_solid){std::string filename;
            /*
            if(fluids_parameters.smoke && fluids_parameters.use_density && fluids_parameters.semi_lagrangian_collidable_density.use_valid_mask){
                filename=output_directory+"/density_valid_mask."+f;
                if(FILE_UTILITIES::File_Exists(filename)){
                    LOG::cout<<"Reading "<<filename<<std::endl;
                    FILE_UTILITIES::Read_From_File(stream_type,filename,fluids_parameters.semi_lagrangian_collidable_density.valid_points_current);}}
            if(fluids_parameters.smoke && fluids_parameters.use_temperature && fluids_parameters.semi_lagrangian_collidable_temperature.use_valid_mask){
                filename=output_directory+"/temperature_valid_mask."+f;
                if(FILE_UTILITIES::File_Exists(filename)){
                    LOG::cout<<"Reading "<<filename<<std::endl;
                    FILE_UTILITIES::Read_From_File(stream_type,filename,fluids_parameters.semi_lagrangian_collidable_temperature.valid_points_current);}}
            if(fluids_parameters.semi_lagrangian_collidable_velocity.use_valid_mask){
                filename=output_directory+"/velocity_valid_mask."+f;
                if(FILE_UTILITIES::File_Exists(filename)){
                    LOG::cout<<"Reading "<<filename<<std::endl;
                    FILE_UTILITIES::Read_From_File(stream_type,filename,fluids_parameters.semi_lagrangian_collidable_velocity.valid_points_current);}}
            if((fluids_parameters.water || fluids_parameters.fire) && fluids_parameters.semi_lagrangian_collidable_phi.use_valid_mask){
                filename=output_directory+"/phi_valid_mask."+f;
                if(FILE_UTILITIES::File_Exists(filename)){
                    LOG::cout<<"Reading "<<filename<<std::endl;
                    FILE_UTILITIES::Read_From_File(stream_type,filename,fluids_parameters.semi_lagrangian_collidable_phi.valid_points_current);}}*/}}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Write_Output_Files(const int frame) const
{
    if(this->use_test_output){
        std::string file=STRING_UTILITIES::string_sprintf("%s/%s-%03d.txt",output_directory.c_str(),this->test_output_prefix.c_str(),frame);
        OCTAVE_OUTPUT<T> oo(file.c_str());
        if(solid_body_collection.deformable_body_collection.particles.X.m){
            oo.Write("db_X",solid_body_collection.deformable_body_collection.particles.X.Flattened());
            oo.Write("db_V",solid_body_collection.deformable_body_collection.particles.V.Flattened());}
        if(solid_body_collection.rigid_body_collection.rigid_body_particles.frame.m){
            RIGID_BODY_PARTICLES<TV>& particle=solid_body_collection.rigid_body_collection.rigid_body_particles;
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
    debug_particles.Write_Debug_Particles(stream_type,output_directory,frame);
    solid_body_collection.Write(stream_type,output_directory,frame,first_frame,solids_parameters.write_static_variables_every_frame,
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,solids_parameters.write_deformable_body,solids_parameters.write_from_every_process,
        solids_parameters.triangle_collision_parameters.output_interaction_pairs);
    if(NEWMARK_EVOLUTION<TV>* newmark=dynamic_cast<NEWMARK_EVOLUTION<TV>*>(solids_evolution))
        newmark->Write_Position_Update_Projection_Data(stream_type,output_directory+"/"+f+"/");
    
    fluids_parameters.Write_Output_Files(stream_type,output_directory,first_frame,frame);
    if(fluids_parameters.incompressible) fluid_collection.Write_Output_Files(stream_type,output_directory,frame);
    if(fluids_parameters.smoke||fluids_parameters.fire||fluids_parameters.water){
        if(fluids_parameters.solid_affects_fluid && fluids_parameters.fluid_affects_solid){
            if(fluids_parameters.write_debug_data){
                /*
                FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/thin_shells_grid_visibility."+f,fluids_parameters.collision_bodies_affecting_fluid->cell_neighbors_visible,
                    fluids_parameters.collision_bodies_affecting_fluid->face_corners_visible_from_face_center);
                if(fluids_parameters.smoke && fluids_parameters.use_density && fluids_parameters.semi_lagrangian_collidable_density.use_valid_mask)
                    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/density_valid_mask."+f,fluids_parameters.semi_lagrangian_collidable_density.valid_points_current);
                if(fluids_parameters.smoke && fluids_parameters.use_temperature && fluids_parameters.semi_lagrangian_collidable_temperature.use_valid_mask)
                    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/temperature_valid_mask."+f,fluids_parameters.semi_lagrangian_collidable_temperature.valid_points_current);
                if(fluids_parameters.semi_lagrangian_collidable_velocity.use_valid_mask) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/velocity_valid_mask."+f,
                    fluids_parameters.semi_lagrangian_collidable_velocity.valid_points_current);
                if((fluids_parameters.water || fluids_parameters.fire) && fluids_parameters.semi_lagrangian_collidable_phi.use_valid_mask)
                    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/phi_valid_mask."+f,fluids_parameters.semi_lagrangian_collidable_phi.valid_points_current);*/}}}
}
//#####################################################################
// Function Log_Parameters 
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Log_Parameters() const
{       
    LOG::SCOPE scope("SOLIDS_FLUIDS_EXAMPLE_UNIFORM parameters");
    BASE::Log_Parameters();
    fluids_parameters.Log_Parameters();
}
//#####################################################################
template SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::~SOLIDS_FLUIDS_EXAMPLE_UNIFORM();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Initialize_Swept_Occupied_Blocks_For_Advection(float,float,float,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Revalidate_Fluid_Scalars();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Revalidate_Fluid_Velocity(ARRAY<float,FACE_INDEX<1> >&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Revalidate_Phi_After_Modify_Levelset();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> >&,PARTICLE_LEVELSET_PARTICLE_TYPE,float);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Get_Object_Velocities(LAPLACE_UNIFORM<GRID<VECTOR<float,2> > >*,ARRAY<float,FACE_INDEX<2> >&,float,float);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Initialize_Swept_Occupied_Blocks_For_Advection(float,float,float,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Log_Parameters() const;
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Parse_Options();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Read_Output_Files_Fluids(int);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Register_Options();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Revalidate_Fluid_Scalars();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Revalidate_Fluid_Velocity(ARRAY<float,FACE_INDEX<2> >&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Revalidate_Phi_After_Modify_Levelset();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Set_Dirichlet_Boundary_Conditions(float);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Set_Ghost_Density_And_Temperature_Inside_Flame_Core();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Write_Output_Files(int) const;
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<VECTOR<float,3> >&,PARTICLE_LEVELSET_PARTICLE_TYPE,float);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Get_Levelset_Velocity(GRID<VECTOR<float,3> > const&,LEVELSET_MULTIPLE<GRID<VECTOR<float,3> > >&,ARRAY<float,FACE_INDEX<3> >&,float) const;
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Get_Object_Velocities(LAPLACE_UNIFORM<GRID<VECTOR<float,3> > >*,ARRAY<float,FACE_INDEX<3> >&,float,float);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Initialize_Compressible_Incompressible_Coupling();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Initialize_MPI();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Initialize_Solid_Fluid_Coupling_After_Grid_Initialization();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Initialize_Swept_Occupied_Blocks_For_Advection(float,float,float,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Log_Parameters() const;
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Parse_Options();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Read_Output_Files_Fluids(int);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Register_Options();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Revalidate_Fluid_Scalars();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Revalidate_Fluid_Velocity(ARRAY<float,FACE_INDEX<3> >&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Revalidate_Phi_After_Modify_Levelset();
template SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::SOLIDS_FLUIDS_EXAMPLE_UNIFORM(STREAM_TYPE,int,FLUIDS_PARAMETERS<GRID<VECTOR<float,3> > >::TYPE);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Set_Dirichlet_Boundary_Conditions(float);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Set_Ghost_Density_And_Temperature_Inside_Flame_Core();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Write_Output_Files(int) const;
template SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::~SOLIDS_FLUIDS_EXAMPLE_UNIFORM();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<VECTOR<float,2> >&,bool,bool);
template SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::SOLIDS_FLUIDS_EXAMPLE_UNIFORM(STREAM_TYPE,int,FLUIDS_PARAMETERS<GRID<VECTOR<float,2> > >::TYPE);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Adjust_Density_And_Temperature_With_Sources<RANGE<VECTOR<float,2> > >(RANGE<VECTOR<float,2> > const&,MATRIX<float,3,3> const&,float,float);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Adjust_Phi_With_Source<RANGE<VECTOR<float,2> > >(RANGE<VECTOR<float,2> > const&,MATRIX<float,3,3> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Adjust_Phi_With_Source<RANGE<VECTOR<float,2> > >(RANGE<VECTOR<float,2> > const&,int,MATRIX<float,3,3> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Get_Source_Reseed_Mask<RANGE<VECTOR<float,2> > >(RANGE<VECTOR<float,2> > const&,MATRIX<float,3,3> const&,ARRAY<bool,VECTOR<int,2> >*&,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Get_Source_Velocities<RANGE<VECTOR<float,2> > >(RANGE<VECTOR<float,2> > const&,MATRIX<float,3,3> const&,VECTOR<float,2> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<VECTOR<float,3> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<VECTOR<float,1> >&,PARTICLE_LEVELSET_PARTICLE_TYPE,float);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Get_Levelset_Velocity(GRID<VECTOR<float,1> > const&,LEVELSET_MULTIPLE<GRID<VECTOR<float,1> > >&,ARRAY<float,FACE_INDEX<1> >&,float) const;
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Get_Object_Velocities(LAPLACE_UNIFORM<GRID<VECTOR<float,1> > >*,ARRAY<float,FACE_INDEX<1> >&,float,float);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Initialize_Compressible_Incompressible_Coupling();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Initialize_MPI();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Initialize_Solid_Fluid_Coupling_After_Grid_Initialization();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Log_Parameters() const;
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Parse_Options();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Read_Output_Files_Fluids(int);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Register_Options();
template SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::SOLIDS_FLUIDS_EXAMPLE_UNIFORM(STREAM_TYPE,int,FLUIDS_PARAMETERS<GRID<VECTOR<float,1> > >::TYPE);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Set_Dirichlet_Boundary_Conditions(float);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Set_Ghost_Density_And_Temperature_Inside_Flame_Core();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Write_Output_Files(int) const;
template SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::~SOLIDS_FLUIDS_EXAMPLE_UNIFORM();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<VECTOR<float,2> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,2> > >::Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<VECTOR<float,2> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<VECTOR<float,3> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<VECTOR<float,3> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Adjust_Density_And_Temperature_With_Sources<CYLINDER<float> >(CYLINDER<float> const&,MATRIX<float,4,4> const&,float,float);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Adjust_Phi_With_Source<CYLINDER<float> >(CYLINDER<float> const&,MATRIX<float,4,4> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Adjust_Phi_With_Source<RANGE<VECTOR<float,3> > >(RANGE<VECTOR<float,3> > const&,MATRIX<float,4,4> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Get_Source_Velocities<CYLINDER<float> >(CYLINDER<float> const&,MATRIX<float,4,4> const&,VECTOR<float,3> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Get_Source_Velocities<RANGE<VECTOR<float,3> > >(RANGE<VECTOR<float,3> > const&,MATRIX<float,4,4> const&,VECTOR<float,3> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<VECTOR<float,1> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,1> > >::Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<VECTOR<float,1> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Adjust_Density_And_Temperature_With_Sources<RANGE<VECTOR<float,3> > >(RANGE<VECTOR<float,3> > const&,MATRIX<float,4,4> const&,float,float);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Adjust_Phi_With_Source<CYLINDER<float> >(CYLINDER<float> const&,int,MATRIX<float,4,4> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Adjust_Phi_With_Source<SPHERE<VECTOR<float,3> > >(SPHERE<VECTOR<float,3> > const&,int,MATRIX<float,4,4> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Get_Source_Reseed_Mask<CYLINDER<float> >(CYLINDER<float> const&,MATRIX<float,4,4> const&,ARRAY<bool,VECTOR<int,3> >*&,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<float,3> > >::Get_Source_Reseed_Mask<SPHERE<VECTOR<float,3> > >(SPHERE<VECTOR<float,3> > const&,MATRIX<float,4,4> const&,ARRAY<bool,VECTOR<int,3> >*&,bool);
template SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::~SOLIDS_FLUIDS_EXAMPLE_UNIFORM();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Initialize_Swept_Occupied_Blocks_For_Advection(double,double,double,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Revalidate_Fluid_Scalars();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Revalidate_Fluid_Velocity(ARRAY<double,FACE_INDEX<1> >&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Revalidate_Phi_After_Modify_Levelset();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<VECTOR<double,2> >&,PARTICLE_LEVELSET_PARTICLE_TYPE,double);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Get_Object_Velocities(LAPLACE_UNIFORM<GRID<VECTOR<double,2> > >*,ARRAY<double,FACE_INDEX<2> >&,double,double);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Initialize_Swept_Occupied_Blocks_For_Advection(double,double,double,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Log_Parameters() const;
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Parse_Options();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Read_Output_Files_Fluids(int);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Register_Options();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Revalidate_Fluid_Scalars();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Revalidate_Fluid_Velocity(ARRAY<double,FACE_INDEX<2> >&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Revalidate_Phi_After_Modify_Levelset();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Set_Dirichlet_Boundary_Conditions(double);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Set_Ghost_Density_And_Temperature_Inside_Flame_Core();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Write_Output_Files(int) const;
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<VECTOR<double,3> >&,PARTICLE_LEVELSET_PARTICLE_TYPE,double);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Get_Levelset_Velocity(GRID<VECTOR<double,3> > const&,LEVELSET_MULTIPLE<GRID<VECTOR<double,3> > >&,ARRAY<double,FACE_INDEX<3> >&,double) const;
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Get_Object_Velocities(LAPLACE_UNIFORM<GRID<VECTOR<double,3> > >*,ARRAY<double,FACE_INDEX<3> >&,double,double);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Initialize_Compressible_Incompressible_Coupling();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Initialize_MPI();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Initialize_Solid_Fluid_Coupling_After_Grid_Initialization();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Initialize_Swept_Occupied_Blocks_For_Advection(double,double,double,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Log_Parameters() const;
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Parse_Options();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Read_Output_Files_Fluids(int);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Register_Options();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Revalidate_Fluid_Scalars();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Revalidate_Fluid_Velocity(ARRAY<double,FACE_INDEX<3> >&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Revalidate_Phi_After_Modify_Levelset();
template SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::SOLIDS_FLUIDS_EXAMPLE_UNIFORM(STREAM_TYPE,int,FLUIDS_PARAMETERS<GRID<VECTOR<double,3> > >::TYPE);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Set_Dirichlet_Boundary_Conditions(double);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Set_Ghost_Density_And_Temperature_Inside_Flame_Core();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Write_Output_Files(int) const;
template SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::~SOLIDS_FLUIDS_EXAMPLE_UNIFORM();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<VECTOR<double,2> >&,bool,bool);
template SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::SOLIDS_FLUIDS_EXAMPLE_UNIFORM(STREAM_TYPE,int,FLUIDS_PARAMETERS<GRID<VECTOR<double,2> > >::TYPE);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Adjust_Density_And_Temperature_With_Sources<RANGE<VECTOR<double,2> > >(RANGE<VECTOR<double,2> > const&,MATRIX<double,3,3> const&,double,double);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Adjust_Phi_With_Source<RANGE<VECTOR<double,2> > >(RANGE<VECTOR<double,2> > const&,MATRIX<double,3,3> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Adjust_Phi_With_Source<RANGE<VECTOR<double,2> > >(RANGE<VECTOR<double,2> > const&,int,MATRIX<double,3,3> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Get_Source_Reseed_Mask<RANGE<VECTOR<double,2> > >(RANGE<VECTOR<double,2> > const&,MATRIX<double,3,3> const&,ARRAY<bool,VECTOR<int,2> >*&,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Get_Source_Velocities<RANGE<VECTOR<double,2> > >(RANGE<VECTOR<double,2> > const&,MATRIX<double,3,3> const&,VECTOR<double,2> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<VECTOR<double,3> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<VECTOR<double,1> >&,PARTICLE_LEVELSET_PARTICLE_TYPE,double);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Get_Levelset_Velocity(GRID<VECTOR<double,1> > const&,LEVELSET_MULTIPLE<GRID<VECTOR<double,1> > >&,ARRAY<double,FACE_INDEX<1> >&,double) const;
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Get_Object_Velocities(LAPLACE_UNIFORM<GRID<VECTOR<double,1> > >*,ARRAY<double,FACE_INDEX<1> >&,double,double);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Initialize_Compressible_Incompressible_Coupling();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Initialize_MPI();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Initialize_Solid_Fluid_Coupling_After_Grid_Initialization();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Log_Parameters() const;
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Parse_Options();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Read_Output_Files_Fluids(int);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Register_Options();
template SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::SOLIDS_FLUIDS_EXAMPLE_UNIFORM(STREAM_TYPE,int,FLUIDS_PARAMETERS<GRID<VECTOR<double,1> > >::TYPE);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Set_Dirichlet_Boundary_Conditions(double);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Set_Ghost_Density_And_Temperature_Inside_Flame_Core();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Write_Output_Files(int) const;
template SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::~SOLIDS_FLUIDS_EXAMPLE_UNIFORM();
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<VECTOR<double,2> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,2> > >::Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<VECTOR<double,2> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<VECTOR<double,3> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<VECTOR<double,3> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Adjust_Density_And_Temperature_With_Sources<CYLINDER<double> >(CYLINDER<double> const&,MATRIX<double,4,4> const&,double,double);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Adjust_Phi_With_Source<CYLINDER<double> >(CYLINDER<double> const&,MATRIX<double,4,4> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Adjust_Phi_With_Source<RANGE<VECTOR<double,3> > >(RANGE<VECTOR<double,3> > const&,MATRIX<double,4,4> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Get_Source_Velocities<CYLINDER<double> >(CYLINDER<double> const&,MATRIX<double,4,4> const&,VECTOR<double,3> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Get_Source_Velocities<RANGE<VECTOR<double,3> > >(RANGE<VECTOR<double,3> > const&,MATRIX<double,4,4> const&,VECTOR<double,3> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<VECTOR<double,1> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,1> > >::Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<VECTOR<double,1> >&,bool,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Adjust_Density_And_Temperature_With_Sources<RANGE<VECTOR<double,3> > >(RANGE<VECTOR<double,3> > const&,MATRIX<double,4,4> const&,double,double);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Adjust_Phi_With_Source<CYLINDER<double> >(CYLINDER<double> const&,int,MATRIX<double,4,4> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Adjust_Phi_With_Source<SPHERE<VECTOR<double,3> > >(SPHERE<VECTOR<double,3> > const&,int,MATRIX<double,4,4> const&);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Get_Source_Reseed_Mask<CYLINDER<double> >(CYLINDER<double> const&,MATRIX<double,4,4> const&,ARRAY<bool,VECTOR<int,3> >*&,bool);
template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<double,3> > >::Get_Source_Reseed_Mask<SPHERE<VECTOR<double,3> > >(SPHERE<VECTOR<double,3> > const&,MATRIX<double,4,4> const&,ARRAY<bool,VECTOR<int,3> >*&,bool);
