//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Jon Gretarsson, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUIDS_PARAMETERS_UNIFORM
//#####################################################################
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/SCOPE.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/Find_Type.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Computations/GRADIENT_UNIFORM.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_PDE/Advection/ADVECTION_HAMILTON_JACOBI_WENO.h>
#include <Grid_PDE/Advection/ADVECTION_MACCORMACK_UNIFORM.h>
#include <Grid_PDE/Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Grid_PDE/Boundaries/BOUNDARY_REFLECTION_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <Geometry/Projection/BOUNDARY_CONDITION_DOUBLE_FINE.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
#include <Incompressible/Collisions_And_Interactions/FLUID_COLLISION_BODY_INACCURATE_UNION.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <Compressible/Compressible_Fluids/COMPRESSIBLE_AUXILIARY_DATA.h>
#include <Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <Fluids/Coupled_Evolution/COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES.h>
#include <Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <Dynamics/Incompressible_Flows/DETONATION_SHOCK_DYNAMICS.h>
#include <Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/TURBULENCE.h>
#include <Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
#include <Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUIDS_PARAMETERS_UNIFORM<TV>::
FLUIDS_PARAMETERS_UNIFORM(const int number_of_regions_input,const typename FLUIDS_PARAMETERS<TV>::TYPE type)
    :FLUIDS_PARAMETERS<TV>(type),mpi_grid(0),particle_levelset_evolution(0),incompressible(0),particle_levelset_evolution_multiple(0),incompressible_multiphase(0),sph_evolution(0),
    maccormack_node_mask(*new ARRAY<bool,TV_INT>),maccormack_cell_mask(*new ARRAY<bool,TV_INT>),maccormack_face_mask(*new ARRAY<bool,FACE_INDEX<TV::m> >),
    maccormack_semi_lagrangian(*new ADVECTION_MACCORMACK_UNIFORM<TV,T,T_ADVECTION_SEMI_LAGRANGIAN_SCALAR>(semi_lagrangian,&maccormack_node_mask,&maccormack_cell_mask,
        &maccormack_face_mask)),euler(0),euler_solid_fluid_coupling_utilities(0),compressible_incompressible_coupling_utilities(0),projection(0),use_reacting_flow(false),
    use_flame_speed_multiplier(false),use_dsd(false),use_psi_R(false),use_levelset_viscosity(false),print_viscosity_matrix(false),use_second_order_pressure(false),
    use_surface_solve(true),projection_scale(1)
{
    Initialize_Number_Of_Regions(number_of_regions_input);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUIDS_PARAMETERS_UNIFORM<TV>::
~FLUIDS_PARAMETERS_UNIFORM()
{
    delete mpi_grid;
    delete &maccormack_semi_lagrangian;
    delete &maccormack_node_mask;
    delete &maccormack_cell_mask;
    delete &maccormack_face_mask;
    delete particle_levelset_evolution;
    delete incompressible;
    // if incompressible_multiphase!=0 it is the sasme as incompressible
    incompressible_multiphase=0;
    delete sph_evolution;
    delete projection;
    delete euler;
    delete euler_solid_fluid_coupling_utilities;
    delete compressible_incompressible_coupling_utilities;
    delete bc_fine;
}
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Initialize_Number_Of_Regions(const int number_of_regions_input)
{
    number_of_regions=number_of_regions_input;
    masses.Resize(number_of_regions);
    densities.Resize(number_of_regions);
    viscosities.Resize(number_of_regions);
    surface_tensions.Resize(VECTOR<int,2>()+number_of_regions);
    dirichlet_regions.Resize(number_of_regions);
    pseudo_dirichlet_regions.Resize(number_of_regions);
    fuel_region.Resize(number_of_regions);
    normal_flame_speeds.Resize(VECTOR<int,2>()+number_of_regions);
    confinement_parameters.Resize(number_of_regions);
    curvature_flame_speeds.Resize(VECTOR<int,2>()+number_of_regions);
    use_multiphase_strain.Resize(number_of_regions);
    elastic_moduli.Resize(number_of_regions);
    plasticity_alphas.Resize(number_of_regions);
    plasticity_gammas.Resize(number_of_regions);
}
//#####################################################################
// Function Initialize_Fluid_Evolution
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Initialize_Fluid_Evolution(ARRAY<T,FACE_INDEX<TV::m> >& incompressible_face_velocities)
{
    if(number_of_regions>=2){ // multiphase
        particle_levelset_evolution_multiple=new PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<TV>(grid->Get_MAC_Grid(),*collision_bodies_affecting_fluid,number_of_ghost_cells);
        if(!projection) projection=new PROJECTION_DYNAMICS_UNIFORM<TV>(*grid,fire,true,false);
        incompressible_multiphase=new INCOMPRESSIBLE_MULTIPHASE_UNIFORM<TV>(grid->Get_MAC_Grid(),*projection);
        phi_boundary_multiphase.Resize(number_of_regions);for(int i=0;i<number_of_regions;i++) phi_boundary_multiphase(i)=&phi_boundary_reflection;
        phi_boundary=0;
        particle_levelset_evolution=particle_levelset_evolution_multiple;
        incompressible=incompressible_multiphase;}
    else if(number_of_regions==1){ // free surface_flow
        particle_levelset_evolution=new PARTICLE_LEVELSET_EVOLUTION_UNIFORM<TV>(*grid,*collision_bodies_affecting_fluid,number_of_ghost_cells,false);
        if(!projection)
            projection=new PROJECTION_DYNAMICS_UNIFORM<TV>(*grid,fire,false,false,use_poisson);
        incompressible=new INCOMPRESSIBLE_UNIFORM<TV>(*grid,*projection);
        phi_boundary=&phi_boundary_water; // override default
        phi_boundary_water.Set_Velocity_Pointer(incompressible_face_velocities);
        boundary_mac_slip.Set_Phi(particle_levelset_evolution->phi);}
    else if(!compressible){ // smoke or sph
        if(!projection)
            projection=new PROJECTION_DYNAMICS_UNIFORM<TV>(*grid,fire,false,false,use_poisson);
        incompressible=new INCOMPRESSIBLE_UNIFORM<TV>(*grid,*projection);}

    if(sph){
        sph_evolution=new SPH_EVOLUTION_UNIFORM<TV>(*grid,*incompressible,*this);
        fluid_boundary=&fluid_boundary_water;}
    if(compressible){
        euler=new EULER_UNIFORM<TV>(*grid);
        euler_solid_fluid_coupling_utilities=new SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>(*euler,mpi_grid);
        if(number_of_regions){
            compressible_incompressible_coupling_utilities=new COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>(*grid,&incompressible_face_velocities,
                &(incompressible->projection.density),&(particle_levelset_evolution->phi));}}

    if(compressible){
        soot_container.Set_Custom_Advection(semi_lagrangian);
        soot_container.Set_Velocity(&euler->euler_projection.face_velocities);
    
        soot_fuel_container.Set_Custom_Advection(semi_lagrangian);
        soot_fuel_container.Set_Velocity(&euler->euler_projection.face_velocities);
    }
    else{
        soot_container.Set_Velocity(&incompressible_face_velocities);
        soot_fuel_container.Set_Velocity(&incompressible_face_velocities);}
    density_container.Set_Velocity(&incompressible_face_velocities);
    density_container.Set_Velocity(&incompressible_face_velocities);
    temperature_container.Set_Velocity(&incompressible_face_velocities);
}
//#####################################################################
// Function Use_Fluid_Coupling_Defaults
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Use_Fluid_Coupling_Defaults()
{
    PHYSBAM_ASSERT(solid_affects_fluid,"solid_affects_fluid should be true when using coupling"); // This might break out-of-date examples
    if(!compressible) soot_container.Use_Semi_Lagrangian_Collidable_Advection(*collision_bodies_affecting_fluid,incompressible->valid_mask);
    if(!compressible) soot_fuel_container.Use_Semi_Lagrangian_Collidable_Advection(*collision_bodies_affecting_fluid,incompressible->valid_mask);
    density_container.Use_Semi_Lagrangian_Collidable_Advection(*collision_bodies_affecting_fluid,incompressible->valid_mask);
    temperature_container.Use_Semi_Lagrangian_Collidable_Advection(*collision_bodies_affecting_fluid,incompressible->valid_mask);
    for(int i=0;i<number_of_regions;i++){
        particle_levelset_evolution->Levelset_Advection(i).Use_Semi_Lagrangian_Collidable_Advection(*collision_bodies_affecting_fluid,collidable_phi_replacement_value,incompressible->valid_mask);}
    /*if(use_slip) incompressible->Use_Semi_Lagrangian_Collidable_Slip_Advection(*collision_bodies_affecting_fluid);
      else */if(!use_reacting_flow) incompressible->Use_Semi_Lagrangian_Collidable_Advection(*collision_bodies_affecting_fluid);
    else incompressible->Use_Semi_Lagrangian_Fire_Multiphase_Collidable_Advection(*collision_bodies_affecting_fluid,incompressible->projection,
        particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple); // TODO: fire does not work with thin shells, need to test for this
    if(use_maccormack_semi_lagrangian_advection){
        if(use_maccormack_for_level_set) for(int i=0;i<number_of_regions;i++) particle_levelset_evolution->Levelset_Advection(i).Use_Maccormack_Advection(maccormack_cell_mask);
        soot_container.Use_Maccormack_Advection(maccormack_cell_mask);
        soot_fuel_container.Use_Maccormack_Advection(maccormack_cell_mask);
        density_container.Use_Maccormack_Advection(maccormack_cell_mask);
        temperature_container.Use_Maccormack_Advection(maccormack_cell_mask);
        if(use_maccormack_for_incompressible) incompressible->Use_Maccormack_Advection(0,0,&maccormack_face_mask);}
    incompressible->collision_body_list=collision_bodies_affecting_fluid;
}
//#####################################################################
// Function Use_No_Fluid_Coupling_Defaults
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Use_No_Fluid_Coupling_Defaults()
{
    if(!compressible) soot_container.Set_Custom_Advection(semi_lagrangian);
    if(!compressible) soot_fuel_container.Set_Custom_Advection(semi_lagrangian);
    density_container.Set_Custom_Advection(semi_lagrangian);
    temperature_container.Set_Custom_Advection(semi_lagrangian);
    for(int i=0;i<number_of_regions;i++){
        if(false && analytic_test){assert(particle_levelset_evolution->runge_kutta_order_levelset==3);particle_levelset_evolution->Levelset_Advection(i).Set_Custom_Advection(hamilton_jacobi_weno);}
        else particle_levelset_evolution->Levelset_Advection(i).Set_Custom_Advection(semi_lagrangian);}
    if(!use_reacting_flow){
        incompressible->Set_Custom_Advection(semi_lagrangian);}
    else incompressible->Use_Semi_Lagrangian_Fire_Multiphase_Advection(incompressible_multiphase->projection,particle_levelset_evolution_multiple->Levelset_Multiple());
    if(use_maccormack_semi_lagrangian_advection){
        if(use_maccormack_for_level_set) for(int i=0;i<number_of_regions;i++) particle_levelset_evolution->Levelset_Advection(i).Set_Custom_Advection(maccormack_semi_lagrangian);
        soot_container.Set_Custom_Advection(maccormack_semi_lagrangian);
        soot_fuel_container.Set_Custom_Advection(maccormack_semi_lagrangian);
        density_container.Set_Custom_Advection(maccormack_semi_lagrangian);
        temperature_container.Set_Custom_Advection(maccormack_semi_lagrangian);
        if(use_maccormack_for_incompressible && !use_reacting_flow) incompressible->Set_Custom_Advection(maccormack_semi_lagrangian);}
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Initialize_Grids()
{
    if(!compressible)
        *grid=grid->Get_MAC_Grid();
    p_grid=*grid;callbacks->Initialize_Fluids_Grids();
}
//#####################################################################
// Function Adjust_Particle_For_Domain_Boundaries
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    // remove this for speed - don't call the function for the other particles
    if(particle_type==PARTICLE_LEVELSET_POSITIVE || particle_type==PARTICLE_LEVELSET_REMOVED_POSITIVE) return;

    TV& X=particles.X(index);TV X_new=X+dt*V;
    T max_collision_distance=0;
    if(number_of_regions==1) particle_levelset_evolution->Particle_Levelset(0).Particle_Collision_Distance(particles.quantized_collision_distance(index));
    else particle_levelset_evolution_multiple->particle_levelset_multiple.Particle_Collision_Distance(particles.quantized_collision_distance(index));
    T min_collision_distance=particle_levelset_evolution->Particle_Levelset(0).min_collision_distance_factor*max_collision_distance;
    TV min_corner=grid->domain.Minimum_Corner(),max_corner=grid->domain.Maximum_Corner();
    for(int axis=0;axis<TV::m;axis++){
        if(domain_walls[axis][0] && X_new[axis]<min_corner[axis]+max_collision_distance){
            T collision_distance=X[axis]-min_corner[axis];
            if(collision_distance>max_collision_distance)collision_distance=X_new[axis]-min_corner[axis];
            collision_distance=max(min_collision_distance,collision_distance);
            X_new[axis]+=max((T)0,min_corner[axis]-X_new[axis]+collision_distance);
            V[axis]=max((T)0,V[axis]);X=X_new-dt*V;}
        if(domain_walls[axis][1] && X_new[axis]>max_corner[axis]-max_collision_distance){
            T collision_distance=max_corner[axis]-X[axis];
            if(collision_distance>max_collision_distance) collision_distance=max_corner[axis]-X_new[axis];
            collision_distance=max(min_collision_distance,collision_distance);
            X_new[axis]-=max((T)0,X_new[axis]-max_corner[axis]+collision_distance);
            V[axis]=min((T)0,V[axis]);X=X_new-dt*V;}}
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Delete_Particles_Inside_Objects(const T time)
{
    for(int i=0;i<number_of_regions;i++){
        PARTICLE_LEVELSET_UNIFORM<TV>* particle_levelset;
        if(number_of_regions==1) particle_levelset=&particle_levelset_evolution->Particle_Levelset(0);
        else particle_levelset=particle_levelset_evolution_multiple->particle_levelset_multiple.particle_levelsets(i);
        Delete_Particles_Inside_Objects<PARTICLE_LEVELSET_PARTICLES<TV> >(particle_levelset->positive_particles,PARTICLE_LEVELSET_POSITIVE,time);
        Delete_Particles_Inside_Objects<PARTICLE_LEVELSET_PARTICLES<TV> >(particle_levelset->negative_particles,PARTICLE_LEVELSET_NEGATIVE,time);
        if(particle_levelset->use_removed_positive_particles)
            Delete_Particles_Inside_Objects<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(particle_levelset->removed_positive_particles,PARTICLE_LEVELSET_REMOVED_POSITIVE,time);
        if(particle_levelset->use_removed_negative_particles)
            Delete_Particles_Inside_Objects<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(particle_levelset->removed_negative_particles,PARTICLE_LEVELSET_REMOVED_NEGATIVE,time);}
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class TV> template<class T_PARTICLES> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Delete_Particles_Inside_Objects(ARRAY<T_PARTICLES*,TV_INT>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{
    for(NODE_ITERATOR<TV> iterator(*grid);iterator.Valid();iterator.Next()){TV_INT block_index=iterator.Node_Index();if(particles(block_index)){
        BLOCK_UNIFORM<TV> block(*grid,block_index);
        COLLISION_GEOMETRY_ID body_id;int aggregate_id;
        T_PARTICLES& block_particles=*particles(block_index);
        if(collision_bodies_affecting_fluid->Occupied_Block(block)){
            // TODO: add back contour value??
            // for(int k=particles.Size()-1;k>=0;k--) if(collision_bodies_affecting_fluid->Inside_Any_Simplex_Of_Any_Body(particles.X(k),contour_value,body_id,aggregate_id)) particles.Delete_Particle(k);
            // TODO(jontg): Shouldn't this delete particles inside the solid, not just on the surface?
            for(int k=block_particles.Size()-1;k>=0;k--) if(collision_bodies_affecting_fluid->Inside_Any_Simplex_Of_Any_Body(block_particles.X(k),body_id,aggregate_id)) block_particles.Delete_Element(k);}
        callbacks->Delete_Particles_Inside_Objects(block_particles,particle_type,time);
        if(block_particles.Size()==0) for(int i=0;i<number_of_regions;i++) particle_levelset_evolution->Particle_Levelset(i).Free_Particle_And_Clear_Pointer(particles(block_index));}}
}
//#####################################################################
// Function Set_Projection
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Set_Projection(PROJECTION_DYNAMICS_UNIFORM<TV>* projection_input)
{
    if(incompressible) PHYSBAM_FATAL_ERROR("projection cannot be set after initialization");
    projection=projection_input;
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Update_Fluid_Parameters(const T dt,const T time)
{
    if(!incompressible) return;
    incompressible->gravity=gravity;
    incompressible->Set_Body_Force(use_body_force);
    incompressible->projection.Use_Non_Zero_Divergence(use_non_zero_divergence);
    incompressible->projection.elliptic_solver->Solve_Neumann_Regions(solve_neumann_regions);
    incompressible->projection.elliptic_solver->solve_single_cell_neumann_regions=solve_single_cell_neumann_regions;
    incompressible->Use_Explicit_Part_Of_Implicit_Viscosity(use_explicit_part_of_implicit_viscosity);
    if(implicit_viscosity && implicit_viscosity_iterations) incompressible->Set_Maximum_Implicit_Viscosity_Iterations(implicit_viscosity_iterations);
    if(use_vorticity_confinement) incompressible->Set_Vorticity_Confinement(confinement_parameter);
    incompressible->Use_Variable_Vorticity_Confinement(use_variable_vorticity_confinement);
    if(use_variable_vorticity_confinement) callbacks->Get_Variable_Vorticity_Confinement(incompressible->variable_vorticity_confinement,time);

    if(number_of_regions>=2){
        for(int i=0;i<number_of_regions;i++) if(use_multiphase_strain(i)){
            incompressible_multiphase->strains(i)->Set_Elastic_Modulus(elastic_moduli(i));
            incompressible_multiphase->strains(i)->Set_Plasticity_Components(plasticity_alphas(i),plasticity_gammas(i));}
        if(use_reacting_flow) incompressible_multiphase->projection.Set_Flame_Speed_Constants(densities,normal_flame_speeds,curvature_flame_speeds);
        // set all the densities, vorticity confinenements etc here
        incompressible_multiphase->projection.Set_Densities(densities);
        incompressible_multiphase->Set_Viscosity(viscosities);
        incompressible_multiphase->Set_Surface_Tension(surface_tensions);
        incompressible_multiphase->Set_Vorticity_Confinement(confinement_parameters);}
    else{
        if(use_strain){
            incompressible->strain->Set_Elastic_Modulus(elastic_modulus);
            incompressible->strain->Set_Plasticity_Components(plasticity_alpha,plasticity_gamma);}
        incompressible->Set_Surface_Tension(surface_tension);
        incompressible->Set_Variable_Surface_Tension(variable_surface_tension);
        if(variable_surface_tension) callbacks->Get_Variable_Surface_Tension(incompressible->variable_surface_tension,time);
        incompressible->Set_Viscosity(viscosity);
        incompressible->Set_Variable_Viscosity(variable_viscosity);
        if(variable_viscosity) callbacks->Get_Variable_Viscosity(incompressible->variable_viscosity,time);
        incompressible->projection.Set_Density(density);}
}
//#####################################################################
// Function Get_Neumann_And_Dirichlet_Boundary_Conditions
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Get_Neumann_And_Dirichlet_Boundary_Conditions(LAPLACE_UNIFORM<TV>* elliptic_solver,
        ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
{
    if(bc_fine)
        return Get_Unified_Boundary_Conditions(elliptic_solver,face_velocities,time);
    elliptic_solver->psi_N.Fill(false);
    elliptic_solver->psi_D.Fill(false);
    if(incompressible){
        POISSON_COLLIDABLE_UNIFORM<TV>* poisson=incompressible->projection.poisson_collidable;
        if(poisson) poisson->beta_face.Fill(1/density);}
    Set_Domain_Boundary_Conditions(*elliptic_solver,face_velocities,time);
    callbacks->Set_Dirichlet_Boundary_Conditions(time);
    callbacks->Get_Source_Velocities(face_velocities,elliptic_solver->psi_N,time);
    callbacks->Get_Object_Velocities(elliptic_solver,face_velocities,dt,time);

    if(solid_affects_fluid && !fluid_affects_solid) for(CELL_ITERATOR<TV> iterator(*grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(elliptic_solver->All_Cell_Faces_Neumann(cell_index)){
            elliptic_solver->psi_D(cell_index)=true;
            elliptic_solver->u(cell_index)=(T)0;}}  // Never flag a solid cell as Dirichlet unless it's surrounded by Neumann BC
}
//#####################################################################
// Function Get_Unified_Boundary_Conditions
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Get_Unified_Boundary_Conditions(LAPLACE_UNIFORM<TV>* elliptic_solver,
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T time)
{
    if(incompressible){
        POISSON_COLLIDABLE_UNIFORM<TV>* poisson=incompressible->projection.poisson_collidable;
        if(poisson) poisson->beta_face.Fill(1/density);}

    // for(int axis=0;axis<TV::m;axis++)
    //     for(int axis_side=0;axis_side<2;axis_side++){
    //         char bc_type=domain_walls(axis)(axis_side)?bc_fine->bc_slip:bc_fine->bc_free;
    //         int side=2*axis+axis_side;
    //         bc_fine->Set_Domain_Walls(1<<side,bc_type,0);}

    callbacks->Get_Unified_Boundary_Conditions(bc_fine,time);
    bc_fine->Get_Pressure_Boundary_Conditions(elliptic_solver->psi_D,elliptic_solver->psi_N,elliptic_solver->u,face_velocities);
}
//#####################################################################
// Function Set_Domain_Boundary_Conditions
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Set_Domain_Boundary_Conditions(LAPLACE_UNIFORM<TV>& elliptic_solver,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T time)
{
    ARRAY<bool,TV_INT>& psi_D=elliptic_solver.psi_D;ARRAY<bool,FACE_INDEX<TV::m> >& psi_N=elliptic_solver.psi_N;

    for(int axis=0;axis<TV::m;axis++) for(int axis_side=0;axis_side<2;axis_side++){
        int side=2*axis+axis_side;
        TV_INT interior_cell_offset=axis_side==0?TV_INT():-TV_INT::Axis_Vector(axis);
        TV_INT exterior_cell_offset=axis_side==0?-TV_INT::Axis_Vector(axis):TV_INT();
        TV_INT boundary_face_offset=axis_side==0?TV_INT::Axis_Vector(axis):-TV_INT::Axis_Vector(axis);
        if(domain_walls(axis)(axis_side)){
            if(number_of_regions>=2){
                // TODO: clean this up. currently iterating over faces 1 grid cell out from boundary in order to get the corners as well. change iterators to give corners instead.
                for(FACE_ITERATOR<TV> iterator(elliptic_solver.grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                    TV_INT boundary_face=iterator.Face_Index()+boundary_face_offset;
                    int region=particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Inside_Region_Face(iterator.Axis(),boundary_face);
                    if(!dirichlet_regions(region) || (flood_fill_for_bubbles && incompressible_multiphase->levelset_for_dirichlet_regions->phi(boundary_face+interior_cell_offset)<=0)){
                        if(face_velocities.Component(axis).Valid_Index(boundary_face)){
                            psi_N.Component(axis)(boundary_face)=true;face_velocities.Component(axis)(boundary_face)=0;}}
                    else{TV_INT cell=boundary_face+exterior_cell_offset;
                        psi_D(cell)=true;elliptic_solver.u(cell)=0;}}}
            else{
                // TODO: clean this up. currently iterating over faces 1 grid cell out from boundary in order to get the corners as well. change iterators to give corners instead.
                for(FACE_ITERATOR<TV> iterator(elliptic_solver.grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                    TV_INT face=iterator.Face_Index()+boundary_face_offset;
                    if(!water || particle_levelset_evolution->phi(face+interior_cell_offset)<=0 || use_sph_for_removed_negative_particles){
                        if(face_velocities.Component(axis).Valid_Index(face)){
                            psi_N.Component(axis)(face)=true;face_velocities.Component(axis)(face)=0;}}
                    else{TV_INT cell=face+exterior_cell_offset;
                        psi_D(cell)=true;elliptic_solver.u(cell)=0;}}}}
        else
            // TODO: clean this up. currently iterating over faces 1 grid cell out from boundary in order to get the corners as well. change iterators to give corners instead.
            for(FACE_ITERATOR<TV> iterator(elliptic_solver.grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
                psi_D(cell)=true;elliptic_solver.u(cell)=0;}}
}
//#####################################################################
// Function Get_Body_Force
//#####################################################################
// flames and smoke have a buoyancy force
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Get_Body_Force(ARRAY<T,FACE_INDEX<TV::m> >& force,const T dt,const T time)
{
    if(smoke){
        ARRAY<T,TV_INT> density_ghost(grid->Domain_Indices(number_of_ghost_cells),no_init),temperature_ghost(grid->Domain_Indices(number_of_ghost_cells),no_init);
        density_container.boundary->Fill_Ghost_Cells_Cell(*grid,density_container.density,density_ghost,time,number_of_ghost_cells);
        temperature_container.boundary->Fill_Ghost_Cells_Cell(*grid,temperature_container.temperature,temperature_ghost,time,number_of_ghost_cells);
        for(FACE_ITERATOR<TV> iterator(*grid,0,GRID<TV>::WHOLE_REGION,-1,1);iterator.Valid();iterator.Next()){ // y-direction forces only
            T rho_atm=rho_bottom+(rho_top-rho_bottom)*(iterator.Location().y-grid->domain.min_corner.y)/(grid->domain.max_corner.y-grid->domain.min_corner.y);
            T face_density=(density_ghost(iterator.First_Cell_Index())+density_ghost(iterator.Second_Cell_Index()))/(T)2;
            T face_temperature=(temperature_ghost(iterator.First_Cell_Index())+temperature_ghost(iterator.Second_Cell_Index()))/(T)2;
            if(face_density>density_buoyancy_threshold){
                T density_difference=face_density-rho_atm,temperature_difference=face_temperature-temperature_container.ambient_temperature;
                if(density_difference>0 || temperature_difference>0)
                    force.Component(1)(iterator.Face_Index())=temperature_buoyancy_constant*temperature_difference-density_buoyancy_constant*density_difference;}}}
    else if(use_reacting_flow){
        ARRAY<T,TV_INT> temperature_ghost(grid->Domain_Indices(number_of_ghost_cells),no_init);
        temperature_container.boundary->Fill_Ghost_Cells_Cell(*grid,temperature_container.temperature,temperature_ghost,time,number_of_ghost_cells);
        LEVELSET_MULTIPLE<TV>& levelset_multiple=particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple;
        ARRAY<T> one_over_densities(number_of_regions);for(int i=0;i<number_of_regions;i++) one_over_densities(i)=(T)1/densities(i);
        for(FACE_ITERATOR<TV> iterator(*grid,0,GRID<TV>::WHOLE_REGION,-1,1);iterator.Valid();iterator.Next()){ // y-direction forces only
            T temperature;
            int region=levelset_multiple.Inside_Region_Face(iterator.Axis(),iterator.Face_Index());
            if(fuel_region(region)) temperature=temperature_fuel;
            else temperature=(T).5*(temperature_ghost(iterator.First_Cell_Index())+temperature_ghost(iterator.Second_Cell_Index()));
            force.Component(1)(iterator.Face_Index())=one_over_densities(region)*temperature_buoyancy_constant*(temperature-temperature_container.ambient_temperature);}}
}
template<> void FLUIDS_PARAMETERS_UNIFORM<VECTOR<double,1> >::
Get_Body_Force(ARRAY<T,FACE_INDEX<1> >& force,const T dt,const T time)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<> void FLUIDS_PARAMETERS_UNIFORM<VECTOR<float,1> >::
Get_Body_Force(ARRAY<T,FACE_INDEX<1> >& force,const T dt,const T time)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Apply_Isobaric_Fix
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Apply_Isobaric_Fix(const T dt,const T time)
{
    if(compressible_apply_isobaric_fix) euler_solid_fluid_coupling_utilities->Apply_Isobaric_Fix(dt,time);
}
//#####################################################################
// Function Blend_In_External_Velocity
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Blend_In_External_Velocity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
{
    if(use_external_velocity){
        ARRAY<TV,TV_INT> V_blend(grid->Domain_Indices(1));ARRAY<T,TV_INT> blend(grid->Domain_Indices(1));callbacks->Get_External_Velocity(V_blend,blend,time);
        for(FACE_ITERATOR<TV> iterator(*grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();T& face_velocity=face_velocities.Component(axis)(iterator.Face_Index());
            T b=(T).5*(blend(iterator.First_Cell_Index())+blend(iterator.Second_Cell_Index()));b=(T)1-pow(max((T)1-b,(T)0),dt);
            face_velocity=(1-b)*face_velocity+b*(T).5*(V_blend(iterator.First_Cell_Index())[axis]+V_blend(iterator.Second_Cell_Index())[axis]);}}
    if(kolmogorov){
        while(time>turbulence.time_end) turbulence.Advance_Turbulence();
        T b=(T)1-pow(max((T)1-kolmogorov,(T)0),dt),fraction=(time-turbulence.time_start)/(turbulence.time_end-turbulence.time_start);
        for(FACE_ITERATOR<TV> iterator(*grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();T& face_velocity=face_velocities.Component(axis)(iterator.Face_Index());
            face_velocity=(1-b)*face_velocity+b*turbulence.Turbulent_Face_Velocity(axis,iterator.Location(),fraction);}}
}
//#####################################################################
// Function Move_Grid
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Move_Grid(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const TV_INT& shift_domain,const T time)
{
    TV_INT max_shift=2*TV_INT::All_Ones_Vector();
    TV_INT remaining_shift=shift_domain;
    while(remaining_shift!=TV_INT()){
        TV_INT temp_shift=clamp(remaining_shift,-max_shift,max_shift);remaining_shift-=temp_shift;
        TV temp_shift_t(temp_shift);

        RANGE<TV> new_domain(grid->domain.Minimum_Corner()+temp_shift_t*grid->dX,grid->domain.Maximum_Corner()+temp_shift_t*grid->dX);
        grid->Initialize(grid->counts,new_domain,grid->Is_MAC_Grid());

        {GRID<TV>& pls_grid(particle_levelset_evolution->grid);
        RANGE<TV> new_domain(pls_grid.domain.Minimum_Corner()+temp_shift_t*pls_grid.dX,pls_grid.domain.Maximum_Corner()+temp_shift_t*pls_grid.dX);
        pls_grid.Initialize(pls_grid.counts,new_domain,pls_grid.Is_MAC_Grid());}

        if(mpi_grid){
            GRID<TV>& global_grid(mpi_grid->global_grid);
            RANGE<TV> new_domain(global_grid.domain.Minimum_Corner()+temp_shift_t*global_grid.dX,global_grid.domain.Maximum_Corner()+temp_shift_t*global_grid.dX);
            global_grid.Initialize(global_grid.counts,new_domain,global_grid.Is_MAC_Grid());}


        {ARRAY<T,TV_INT> phi_ghost(grid->Domain_Indices(number_of_ghost_cells),no_init);
            phi_boundary->Fill_Ghost_Cells(*grid,particle_levelset_evolution->phi,phi_ghost,0,time,number_of_ghost_cells);
        ARRAY<T,TV_INT>::Limited_Shifted_Get(particle_levelset_evolution->phi,phi_ghost,temp_shift);}
        {ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost(*grid,number_of_ghost_cells,no_init);
            fluid_boundary->Fill_Ghost_Faces(*grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);
        for(int axis=0;axis<TV::m;axis++)
            ARRAY<T,TV_INT>::Limited_Shifted_Get(face_velocities.Component(axis),face_velocities_ghost.Component(axis),temp_shift);}
        {ARRAY<T,TV_INT> p_ghost(p_grid.Domain_Indices(number_of_ghost_cells),no_init);
            BOUNDARY<TV,T>().Fill_Ghost_Cells(p_grid,incompressible->projection.p,p_ghost,0,time,number_of_ghost_cells);
        ARRAY<T,TV_INT>::Limited_Shifted_Get(incompressible->projection.p,p_ghost,temp_shift);}

        for(int i=0;i<number_of_regions;i++){
            PARTICLE_LEVELSET_UNIFORM<TV>& particle_levelset=particle_levelset_evolution->Particle_Levelset(i);
            particle_levelset.Update_Particle_Cells(particle_levelset.positive_particles);
            particle_levelset.Update_Particle_Cells(particle_levelset.negative_particles);
            if(particle_levelset.use_removed_positive_particles) particle_levelset.Update_Particle_Cells(particle_levelset.removed_positive_particles);
            if(particle_levelset.use_removed_negative_particles) particle_levelset.Update_Particle_Cells(particle_levelset.removed_negative_particles);
            particle_levelset.Delete_Particles_Outside_Grid();}}
}
//#####################################################################
// Function Move_Grid
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Move_Grid(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T time)
{
    assert(!smoke && !fire);

    if(move_grid_explicitly){callbacks->Move_Grid_Explicitly(time);return;} // otherwise move automatically

    TV_INT shift;
    for(int axis=0;axis<TV::m;axis++) for(int axis_side=0;axis_side<2;axis_side++){
        int side=2*axis+axis_side;TV_INT offset=(axis_side==0?-1:1)*moving_grid_number_of_cells*TV_INT::Axis_Vector(axis);
        // loop over boundary region by looping over ghost region and shifting inwards
        for(CELL_ITERATOR<TV> iterator(*grid,moving_grid_number_of_cells,GRID<TV>::GHOST_REGION,side);iterator.Valid();iterator.Next())
            if(particle_levelset_evolution->phi(iterator.Cell_Index()-offset)<=0){shift+=offset;break;}}

    Move_Grid(face_velocities,shift,time);
    *grid=GRID<TV>(grid->Domain_Indices().Maximum_Corner(),grid->domain+(TV)shift*grid->dX);
}
//#####################################################################
// Function Adjust_Strain_For_Object
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Adjust_Strain_For_Object(LEVELSET<TV>& levelset_object,ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT>& e_ghost,const T time)
{
    assert(!smoke && !fire);
    if(adhesion_coefficient==1 && !adhesion_normal_strain) return;
    T epsilon=adhesion_half_bandwidth*grid->dX.Min();
    for(CELL_ITERATOR<TV> iterator(*grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        T heaviside=LEVELSET_UTILITIES<T>::Heaviside(levelset_object.phi(cell),epsilon);
        e_ghost(cell)*=adhesion_coefficient+(1-adhesion_coefficient)*heaviside;
        if(adhesion_normal_strain && heaviside!=1) e_ghost(cell)+=(1-heaviside)*adhesion_normal_strain*SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(levelset_object.Normal(iterator.Location()));}
}
//#####################################################################
// Function Combustion
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Combustion(const T dt,const T time)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Before Combustion",1);

    for(CELL_ITERATOR<TV> iterator(*grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        T fuel_fraction=soot_fuel_container.density(cell_index);
        T cell_temperature=euler->Get_Temperature(cell_index);
        T rho=euler->U(cell_index)(0);
        T fuel_burnt=0;
        if(cell_temperature>burn_temperature_threshold && fuel_fraction>0){
            if(rho*fuel_fraction>burn_rate*dt) fuel_burnt=burn_rate*dt;
            else fuel_burnt=rho*fuel_fraction;
            T fuel_fraction_burnt=fuel_burnt/rho;

            T energy_generated=fuel_burnt*soot_fuel_calorific_value;

            soot_fuel_container.density(cell_index)-=fuel_fraction_burnt;
            soot_container.density(cell_index)+=fuel_fraction_burnt;
            euler->U(cell_index)(TV::m+1)+=energy_generated;}}
    euler->Invalidate_Ghost_Cells();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After Combustion",1);
}
//#####################################################################
// Function Evolve_Soot
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Evolve_Soot(const T dt,const T time)
{
    if(compressible){
        if(use_soot_fuel_combustion) Combustion(dt,time);
        euler->euler_projection.Compute_Density_Weighted_Face_Velocities(dt,time,euler->euler_projection.elliptic_solver->psi_N);}
    BASE::Evolve_Soot(dt,time);
}
//#####################################################################
// Function Total_Number_Of_Particles
//#####################################################################
template<class TV> template<class T_ARRAYS_PARTICLES> int FLUIDS_PARAMETERS_UNIFORM<TV>::
Total_Number_Of_Particles(const T_ARRAYS_PARTICLES& particles) const
{
    int total=0;
    for(CELL_ITERATOR<TV> iterator(*grid,3);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())) total+=particles(iterator.Cell_Index())->Size();
    return total;
}
//#####################################################################
// Function Write_Particles
//#####################################################################
template<class TV> template<class T_ARRAYS_PARTICLES> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Write_Particles(const STREAM_TYPE stream_type,const PARTICLES<TV>& template_particles,const T_ARRAYS_PARTICLES& particles,const VIEWER_DIR& viewer_dir,const std::string& prefix) const
{
    Write_To_File(stream_type,viewer_dir.current_directory+"/"+prefix,particles);
    int total_number_of_particles=Total_Number_Of_Particles(particles);
    LOG::cout<<"Writing "<<total_number_of_particles<<" "<<prefix<<std::endl;

    // write out one T_PARTICLES that contains all the particles
    if(write_flattened_particles){
        PARTICLES<TV>* all_particles=template_particles.Clone();
        all_particles->Initialize(particles.array);
        Write_To_File(stream_type,viewer_dir.current_directory+"/"+prefix+"_all",*all_particles);}
}
//#####################################################################
// Function Read_Particles
//#####################################################################
template<class TV> template<class T_PARTICLES,class T_ARRAYS_PARTICLES> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Read_Particles(const T_PARTICLES& template_particles,T_ARRAYS_PARTICLES& particles,const VIEWER_DIR& viewer_dir,const std::string& prefix)
{
    STATIC_ASSERT((is_same<T_PARTICLES,typename remove_pointer<typename T_ARRAYS_PARTICLES::ELEMENT>::type>::value)); 
    Read_From_File(viewer_dir.current_directory+"/"+prefix,particles);
    if(typeid(template_particles)!=typeid(T_PARTICLES)) // swap in clones of template_particle for pure T_PARTICLESs
        for(int i=0;i<particles.array.Size();i++) if(particles.array(i)){
            T_PARTICLES* replacement=template_particles.Clone();
            replacement->Initialize(*particles.array(i));
            delete particles.array(i);
            particles.array(i)=replacement;}
    LOG::cout<<"Reading "<<Total_Number_Of_Particles(particles)<<" "<<prefix<<std::endl;
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Read_Output_Files(const VIEWER_DIR& viewer_dir)
{
    Read_From_File(viewer_dir.current_directory+"/grid",*grid);
    if(mpi_grid) Read_From_File(viewer_dir.current_directory+"/global_grid",mpi_grid->global_grid);

    if(use_soot) Read_From_File(viewer_dir.current_directory+"/soot",soot_container.density);
    if(use_soot && use_soot_fuel_combustion) Read_From_File(viewer_dir.current_directory+"/soot_fuel",soot_fuel_container.density);
    if(smoke || fire || water){
        // scalar fields
        if(smoke || fire){
            if(use_density) Read_From_File(viewer_dir.current_directory+"/density",density_container.density);
            if(use_temperature) Read_From_File(viewer_dir.current_directory+"/temperature",temperature_container.temperature);}
        // particle levelset
        if(write_levelset){
            if(number_of_regions==1){
                PARTICLE_LEVELSET_UNIFORM<TV>& particle_levelset=particle_levelset_evolution->Particle_Levelset(0);
                Read_From_File(viewer_dir.current_directory+"/levelset",particle_levelset.levelset);
                if(write_particles && viewer_dir.frame_stack(0)%restart_data_write_rate==0){
                    Read_Particles(particle_levelset.template_particles,particle_levelset.positive_particles,viewer_dir,"positive_particles");
                    Read_Particles(particle_levelset.template_particles,particle_levelset.negative_particles,viewer_dir,"negative_particles");}
                if(write_removed_positive_particles)
                    Read_Particles(particle_levelset.template_removed_particles,particle_levelset.removed_positive_particles,viewer_dir,
                        "removed_positive_particles");
                if(write_removed_negative_particles)
                    Read_Particles(particle_levelset.template_removed_particles,particle_levelset.removed_negative_particles,viewer_dir,
                        "removed_negative_particles");
                if(store_particle_ids)
                    Read_From_Text_File(viewer_dir.current_directory+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
                if(use_strain && write_strain) Read_From_File(viewer_dir.current_directory+"/strain",incompressible->strain->e);}
            else if(number_of_regions>=2){
                for(int i=0;i<number_of_regions;i++){
                    std::string ii=LOG::sprintf("%d",i); // TODO(jontg): This still does .%d.gz
                    PARTICLE_LEVELSET_UNIFORM<TV>& particle_levelset=*particle_levelset_evolution_multiple->particle_levelset_multiple.particle_levelsets(i);
                    Read_From_File(viewer_dir.current_directory+"/levelset_"+ii,particle_levelset.levelset);
                    if(write_particles && viewer_dir.frame_stack(0)%restart_data_write_rate==0){
                        Read_Particles(particle_levelset.template_particles,particle_levelset.positive_particles,viewer_dir,"positive_particles_"+ii);
                        Read_Particles(particle_levelset.template_particles,particle_levelset.negative_particles,viewer_dir,"negative_particles_"+ii);}
                    if(write_removed_positive_particles)
                        Read_Particles(particle_levelset.template_removed_particles,particle_levelset.removed_positive_particles,viewer_dir,
                            "removed_positive_particles_"+ii);
                    if(write_removed_negative_particles)
                        Read_Particles(particle_levelset.template_removed_particles,particle_levelset.removed_negative_particles,viewer_dir,
                            "removed_negative_particles_"+ii);
                    //if(store_particle_ids)
                    //    Read_From_Text_File(viewer_dir.current_directory+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
                    if(write_strain && use_multiphase_strain.Count_Matches(0)<number_of_regions) if(incompressible_multiphase->strains(i)){
                        if(incompressible_multiphase->strains(i))
                            Read_From_File(viewer_dir.current_directory+"/strain_"+ii,incompressible_multiphase->strains(i)->e);}}}}

        // pressure and velocities
        std::string filename;
        filename=viewer_dir.current_directory+"/pressure";
        if(File_Exists(filename)){LOG::cout<<"Reading pressure "<<filename<<std::endl;
            Read_From_File(filename,incompressible->projection.p);}}

    else if(compressible){
        Read_From_File(viewer_dir.current_directory+"/euler_U",euler->U);
        Read_From_File(viewer_dir.current_directory+"/euler_psi",euler->psi);}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Write_Output_Files(const STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir) const
{
    if(!simulate) return;
    if(use_soot){
        ARRAY<T,TV_INT> soot_ghost(grid->Cell_Indices(number_of_ghost_cells),no_init);
        ARRAY<T,TV_INT> soot_fuel_ghost(grid->Cell_Indices(number_of_ghost_cells),no_init);
        soot_container.boundary->Fill_Ghost_Cells(*grid,soot_container.density,soot_ghost,0,0,number_of_ghost_cells);
        soot_fuel_container.boundary->Fill_Ghost_Cells(*grid,soot_fuel_container.density,soot_fuel_ghost,0,0,number_of_ghost_cells);
        Write_To_File(stream_type,viewer_dir.current_directory+"/soot",soot_ghost);
        if(use_soot_fuel_combustion) Write_To_File(stream_type,viewer_dir.current_directory+"/soot_fuel",soot_fuel_ghost);}
    if(smoke || fire || water || sph){
        if(viewer_dir.First_Frame()){Write_To_File(stream_type,viewer_dir.output_directory+"/common/grid",*grid);
        if(mpi_grid) Write_To_File(stream_type,viewer_dir.output_directory+"/common/global_grid",mpi_grid->global_grid);}
        Write_To_File(stream_type,viewer_dir.current_directory+"/grid",*grid);
        if(mpi_grid) Write_To_File(stream_type,viewer_dir.current_directory+"/global_grid",mpi_grid->global_grid);
        // object levelset
        if(FLUID_COLLISION_BODY_INACCURATE_UNION<TV>* innacurate_union=
            Find_Type<FLUID_COLLISION_BODY_INACCURATE_UNION<TV>*>(collision_bodies_affecting_fluid->collision_geometry_collection.bodies))
            Write_To_File(stream_type,viewer_dir.current_directory+"/object_levelset",innacurate_union->levelset);
        // scalar fields
        if(smoke || fire){
            if(use_density){
                ARRAY<T,TV_INT> density_ghost(grid->Cell_Indices(number_of_ghost_cells),no_init);
                density_container.boundary->Fill_Ghost_Cells(*grid,density_container.density,density_ghost,0,0,number_of_ghost_cells);
                Write_To_File(stream_type,viewer_dir.current_directory+"/density",density_ghost);}
            if(use_temperature){
                ARRAY<T,TV_INT> temperature_ghost(grid->Cell_Indices(number_of_ghost_cells),no_init);
                temperature_container.boundary->Fill_Ghost_Cells(*grid,temperature_container.temperature,temperature_ghost,0,0,number_of_ghost_cells);
                Write_To_File(stream_type,viewer_dir.current_directory+"/temperature",temperature_ghost);}}
        // sph particles
        if(sph) Write_To_File(stream_type,viewer_dir.current_directory+"/sph_particles",sph_evolution->sph_particles);
        // velocities
        if(write_velocity && viewer_dir.frame_stack(0)%restart_data_write_rate==0){
            if(fire) Write_To_File(stream_type,viewer_dir.current_directory+"/levelset_velocities",particle_levelset_evolution->V);}
        // particle levelset
        if(write_levelset){
            if(number_of_regions==1){
                PARTICLE_LEVELSET_UNIFORM<TV>& particle_levelset=particle_levelset_evolution->Particle_Levelset(0);
                Write_To_File(stream_type,viewer_dir.current_directory+"/levelset",particle_levelset.levelset);
                if(write_particles && viewer_dir.frame_stack(0)%restart_data_write_rate==0){
                    Write_Particles(stream_type,particle_levelset.template_particles,particle_levelset.positive_particles,viewer_dir,"positive_particles");
                    Write_Particles(stream_type,particle_levelset.template_particles,particle_levelset.negative_particles,viewer_dir,"negative_particles");}
                if(write_removed_positive_particles)
                    Write_Particles(stream_type,particle_levelset.template_removed_particles,particle_levelset.removed_positive_particles,viewer_dir,
                        "removed_positive_particles");
                if(write_removed_negative_particles)
                    Write_Particles(stream_type,particle_levelset.template_removed_particles,particle_levelset.removed_negative_particles,viewer_dir,
                        "removed_negative_particles");
                if(store_particle_ids)
                    Write_To_Text_File(viewer_dir.current_directory+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);}
            else if(number_of_regions>=2){
                for(int i=0;i<number_of_regions;i++){
                    std::string ii=LOG::sprintf("%d",i);
                    PARTICLE_LEVELSET_UNIFORM<TV>& particle_levelset=*particle_levelset_evolution_multiple->particle_levelset_multiple.particle_levelsets(i);
                    Write_To_File(stream_type,viewer_dir.current_directory+"/levelset_"+ii,particle_levelset.levelset);
                    if(write_particles && viewer_dir.frame_stack(0)%restart_data_write_rate==0){
                        Write_Particles(stream_type,particle_levelset.template_particles,particle_levelset.positive_particles,viewer_dir,"positive_particles_"+ii);
                        Write_Particles(stream_type,particle_levelset.template_particles,particle_levelset.negative_particles,viewer_dir,"negative_particles_"+ii);}
                    if(write_removed_positive_particles)
                        Write_Particles(stream_type,particle_levelset.template_removed_particles,particle_levelset.removed_positive_particles,viewer_dir,
                            "removed_positive_particles_"+ii);
                    if(write_removed_negative_particles)
                        Write_Particles(stream_type,particle_levelset.template_removed_particles,particle_levelset.removed_negative_particles,viewer_dir,
                            "removed_negative_particles_"+ii);
                    //if(store_particle_ids)
                    //    Write_To_Text_File(viewer_dir.current_directory+"/last_unique_particle_id."+f,particle_levelset.last_unique_particle_id);
                }}}
        if(water){
            if(write_strain){
                if(use_strain && number_of_regions==1) Write_To_File(stream_type,viewer_dir.current_directory+"/strain",incompressible->strain->e);
                if(use_multiphase_strain.Count_Matches(0)<number_of_regions && number_of_regions>=2)
                    for(int i=0;i<number_of_regions;i++){
                        if(incompressible_multiphase->strains(i)){
                            std::string ii=LOG::sprintf("%d",i); // TODO(jontg): ...
                            ARRAY<SYMMETRIC_MATRIX<T,TV::m>,TV_INT> e_ghost(grid->Domain_Indices(number_of_ghost_cells),no_init);
                            incompressible_multiphase->strains(i)->e_boundary->Fill_Ghost_Cells(*grid,incompressible_multiphase->strains(i)->e,e_ghost,0,0,number_of_ghost_cells); // TODO: use real dt/time
                            Write_To_File(stream_type,viewer_dir.current_directory+"/strain_"+ii,e_ghost);}}}}

        // sph
        if(sph_evolution) Write_To_File(stream_type,viewer_dir.current_directory+"/sph_cell_weights",sph_evolution->cell_weight);

        // dsd
        if(use_dsd){
            ARRAY<T,TV_INT> reaction_speed_ghost(grid->Domain_Indices(number_of_ghost_cells));incompressible->boundary->Fill_Ghost_Cells(*grid,incompressible->projection.dsd->Dn.array,reaction_speed_ghost,0,0,number_of_ghost_cells); // TODO: use real dt/time
            Write_To_File(stream_type,viewer_dir.current_directory+"/reaction_speed",reaction_speed_ghost);}

        if(write_debug_data || (write_restart_data && viewer_dir.frame_stack(0)%restart_data_write_rate==0)){
            Write_To_File(stream_type,viewer_dir.current_directory+"/pressure",incompressible->projection.p);}
        // debugging
        if(write_debug_data){
            if(incompressible->projection.poisson) Write_To_File(stream_type,viewer_dir.current_directory+"/beta_face",incompressible->projection.poisson->beta_face);
            if(use_body_force) Write_To_File(stream_type,viewer_dir.current_directory+"/forces",incompressible->force);
            if(use_maccormack_semi_lagrangian_advection){
                Write_To_File(stream_type,viewer_dir.current_directory+"/maccormack_cell_mask",maccormack_cell_mask);
                Write_To_File(stream_type,viewer_dir.current_directory+"/maccormack_face_mask",maccormack_face_mask);}
            Write_To_File(stream_type,viewer_dir.current_directory+"/psi_N",incompressible->projection.elliptic_solver->psi_N);
            Write_To_File(stream_type,viewer_dir.current_directory+"/psi_D",incompressible->projection.elliptic_solver->psi_D);
            Write_To_File(stream_type,viewer_dir.current_directory+"/colors",incompressible->projection.elliptic_solver->filled_region_colors);}}
    else if(compressible){
        if(viewer_dir.First_Frame()){
            Write_To_File(stream_type,viewer_dir.output_directory+"/common/grid",euler->grid);
            if(mpi_grid) Write_To_File(stream_type,viewer_dir.output_directory+"/common/global_grid",mpi_grid->global_grid);}
        Write_To_File(stream_type,viewer_dir.current_directory+"/grid",euler->grid);
        if(mpi_grid) Write_To_File(stream_type,viewer_dir.current_directory+"/global_grid",mpi_grid->global_grid);

        Write_To_File(stream_type,viewer_dir.current_directory+"/euler_U",euler->U);
        Write_To_File(stream_type,viewer_dir.current_directory+"/euler_psi",euler->psi);
        euler->Fill_Ghost_Cells(0,0,number_of_ghost_cells);
        COMPRESSIBLE_AUXILIARY_DATA::Write_Auxiliary_Files(stream_type,viewer_dir,euler->grid,3,euler->U_ghost,euler->psi,*euler->eos,write_debug_data,&(euler->conservation->fluxes));
        if(write_debug_data && euler->timesplit) Write_To_File(stream_type,viewer_dir.current_directory+"/compressible_implicit_pressure",euler->euler_projection.p);
        if(write_debug_data){
            Write_To_File(stream_type,viewer_dir.current_directory+"/mac_velocities",euler->euler_projection.face_velocities);
            Write_To_File(stream_type,viewer_dir.current_directory+"/psi_D",euler->euler_projection.elliptic_solver->psi_D);
            Write_To_File(stream_type,viewer_dir.current_directory+"/psi_N",euler->euler_projection.elliptic_solver->psi_N);
            if(euler->apply_cavitation_correction){
                Write_To_File(stream_type,viewer_dir.current_directory+"/p_cavitation",euler->euler_cavitation_density.p_cavitation);
                Write_To_File(stream_type,viewer_dir.current_directory+"/p_internal_energy",euler->euler_cavitation_internal_energy.p_cavitation);}}

        if(number_of_regions && write_levelset){
            Write_To_File(stream_type,viewer_dir.current_directory+"/levelset",particle_levelset_evolution->Particle_Levelset(0).levelset);}}
}
//#####################################################################
// Function Log_Parameters 
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Log_Parameters() const
{
    LOG::SCOPE scope("FLUIDS_PARAMETERS_UNIFORM parameters");
    BASE::Log_Parameters();
    if(euler) euler->Log_Parameters();
}
//#####################################################################
// Function Use_Unified_Boundary_Conditions
//#####################################################################
template<class TV> void FLUIDS_PARAMETERS_UNIFORM<TV>::
Use_Unified_Boundary_Conditions()
{
    bc_fine=new BOUNDARY_CONDITION_DOUBLE_FINE<TV>(*grid,number_of_ghost_cells);
}
//#####################################################################
template class FLUIDS_PARAMETERS_UNIFORM<VECTOR<float,1> >;
template class FLUIDS_PARAMETERS_UNIFORM<VECTOR<float,2> >;
template class FLUIDS_PARAMETERS_UNIFORM<VECTOR<float,3> >;
template class FLUIDS_PARAMETERS_UNIFORM<VECTOR<double,1> >;
template class FLUIDS_PARAMETERS_UNIFORM<VECTOR<double,2> >;
template class FLUIDS_PARAMETERS_UNIFORM<VECTOR<double,3> >;
}
