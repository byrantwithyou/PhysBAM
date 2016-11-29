//#####################################################################
// Copyright 2002-2010, Mridul Aanjaneya, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Duc Nguyen, Nick Rasmussen, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Computations/VORTICITY_UNIFORM.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Grid_PDE/Advection/ADVECTION_MACCORMACK_UNIFORM.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <Incompressible/Advection_Collidable/ADVECTION_WRAPPER_COLLIDABLE_CELL.h>
#include <Incompressible/Advection_Collidable/ADVECTION_WRAPPER_COLLIDABLE_FACE.h>
#include <Incompressible/Advection_Collidable/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <Incompressible/Boundaries/BOUNDARY_MAC_GRID_SOLID_WALL_SLIP.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Forces/VISCOSITY.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/IMPLICIT_VISCOSITY_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/LEVELSET_VISCOSITY_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INCOMPRESSIBLE_UNIFORM<TV>::
INCOMPRESSIBLE_UNIFORM(const GRID<TV>& grid_input,PROJECTION_DYNAMICS_UNIFORM<TV>& projection_input,THREAD_QUEUE* thread_queue_input)
    :grid(grid_input.Get_MAC_Grid()),mpi_grid(0),projection(projection_input),strain(0),collision_body_list(0),momentum_conserving_vorticity(false),use_vorticity_weights(false),energy_clamp(0),vc_projection_direction(0),buoyancy_constant(0),thread_queue(thread_queue_input),
    boundary_default(*new BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV>),advection_maccormack(0)
{ 
    boundary=&boundary_default;
    Initialize_Grids(grid);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INCOMPRESSIBLE_UNIFORM<TV>::
~INCOMPRESSIBLE_UNIFORM()
{
    delete strain;delete &boundary_default;
    delete advection_maccormack;
}
//#####################################################################
// Function Advance_One_Time_Step_Convection
//#####################################################################
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Advance_One_Time_Step_Convection(const T dt,const T time,const ARRAY<T,FACE_INDEX<TV::m> >& advecting_face_velocities,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_to_advect,const int number_of_ghost_cells)
{
    // TODO: make efficient if advection velocities are same as advected velocities
    // find ghost cells
    ARRAY<T,FACE_INDEX<TV::m> > advection_face_velocities_ghost(grid,number_of_ghost_cells,false);
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_to_advect_ghost(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Faces(grid,face_velocities_to_advect,face_velocities_to_advect_ghost,time,number_of_ghost_cells);
    boundary->Fill_Ghost_Faces(grid,advecting_face_velocities,advection_face_velocities_ghost,time,number_of_ghost_cells);

    // update convection
    advection->Update_Advection_Equation_Face(grid,face_velocities_to_advect,face_velocities_to_advect_ghost,advection_face_velocities_ghost,*boundary,dt,time);
}
//#####################################################################
// Function Advance_One_Time_Step_Forces
//#####################################################################
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Advance_One_Time_Step_Forces(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time,const bool implicit_viscosity,const ARRAY<T,TV_INT>* phi_ghost,const int number_of_ghost_cells)
{
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost;face_velocities_ghost.Resize(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Faces(grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);

    // update strain and apply elastic forces
    if(strain){
        assert(!projection.flame);assert(phi_ghost);strain->Update_Strain_Equation(dt,time,projection.density,face_velocities,face_velocities_ghost,*phi_ghost,number_of_ghost_cells);}

    // update gravity
    for(int axis=0;axis<TV::m;axis++)
        if(gravity(axis))
            DOMAIN_ITERATOR_THREADED_ALPHA<INCOMPRESSIBLE_UNIFORM<TV>,TV>(grid.Face_Indices()[axis],thread_queue).template Run<ARRAY<T,FACE_INDEX<TV::m> >&,const T,int>(*this,&INCOMPRESSIBLE_UNIFORM<TV>::Add_Gravity_Threaded,face_velocities,dt,axis);

    // update body force
    if(use_force){
        for(int axis=0;axis<TV::m;axis++)
            DOMAIN_ITERATOR_THREADED_ALPHA<INCOMPRESSIBLE_UNIFORM<TV>,TV>(grid.Face_Indices()[axis],thread_queue).template Run<ARRAY<T,FACE_INDEX<TV::m> >&,const T,int>(*this,&INCOMPRESSIBLE_UNIFORM<TV>::Add_Body_Force_Threaded,face_velocities,dt,axis);
        boundary->Apply_Boundary_Condition_Face(grid,face_velocities,time+dt);
        boundary->Fill_Ghost_Faces(grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);}

    // update viscosity explicitly
    //if(dt && (viscosity || use_variable_viscosity) && (!implicit_viscosity || use_explicit_part_of_implicit_viscosity)){
    //    if(!implicit_viscosity) PHYSBAM_NOT_IMPLEMENTED();
    //    IMPLICIT_VISCOSITY_UNIFORM<TV>::Variable_Viscosity_Explicit_Part(projection.density,variable_viscosity,grid,face_velocities,face_velocities_ghost,dt,time);}
    
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_old=face_velocities;
    if(vorticity_confinement || use_variable_vorticity_confinement){
        T tolerance=0; //TODO (mlentine): Look into what are good values here
        ARRAY<TV,TV_INT> F(grid.Cell_Indices(1),false);
        Compute_Vorticity_Confinement_Force(grid,face_velocities_ghost,F);
        if(collision_body_list){
            if(use_variable_vorticity_confinement){F*=dt;F*=variable_vorticity_confinement;}else F*=dt*vorticity_confinement;}
        else{
            if(use_variable_vorticity_confinement){F*=dt*(T).5;F*=variable_vorticity_confinement;}else F*=dt*vorticity_confinement*(T).5;}
        for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            for(int i=0;i<TV::m;i++) if(abs(F(cell)(i))<tolerance) F(cell)(i)=0;}
        Apply_Vorticity_Confinement_Force(face_velocities,F);}

    boundary->Apply_Boundary_Condition_Face(grid,face_velocities,time+dt);
}
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Add_Gravity_Threaded(RANGE<TV_INT>& domain,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,int axis)
{
    for(FACE_ITERATOR<TV> iterator(grid,domain,axis);iterator.Valid();iterator.Next())
        if(gravity[axis])
            face_velocities.Component(axis)(iterator.Face_Index())+=dt*gravity[axis];
}
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Add_Body_Force_Threaded(RANGE<TV_INT>& domain,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,int axis)
{
    for(FACE_ITERATOR<TV> iterator(grid,domain,axis);iterator.Valid();iterator.Next()){
        face_velocities.Component(iterator.Axis())(iterator.Face_Index())+=dt*force.Component(iterator.Axis())(iterator.Face_Index());}
}
//#####################################################################
// Function Update_Kinetic_Energy
//#####################################################################
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Apply_Pressure_Kinetic_Energy(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_new,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_old,const T dt,const T time)
{
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities(grid);
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::m> face=iterator.Full_Index();face_velocities(face)=(T)((face_velocities_new(face)+face_velocities_old(face))/2.);}
}
//#####################################################################
// Function Correct_Kinetic_Energy
//#####################################################################
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Add_Energy_With_Vorticity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const VECTOR<VECTOR<bool,2>,TV::m>& domain_boundary,const T dt,const T time,const int number_of_ghost_cells,LEVELSET<TV>* lsv,ARRAY<T,TV_INT>* density)
{
    use_vorticity_weights=true;
    vorticity_weights.Resize(grid);
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::m> index=iterator.Full_Index();
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner(iterator.Axis())++;
        for(int i=0;i<TV::m;i++){if(domain_boundary(i)(0)) domain.min_corner(i)++;if(domain_boundary(i)(1)) domain.max_corner(i)--;}
        vorticity_weights(index)=1;
        if(projection.elliptic_solver->psi_N(index) || !domain.Lazy_Inside_Half_Open(index.index)) vorticity_weights(index)=0;}
}
//#####################################################################
// Function Advance_One_Time_Step_Implicit_Part
//#####################################################################
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Advance_One_Time_Step_Implicit_Part(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time,const bool implicit_viscosity,BOUNDARY<TV,T>* projection_boundary,
    bool use_levelset_viscosity,BOUNDARY_CONDITIONS_CALLBACKS<TV>* bc_callbacks,bool print_viscosity_matrix)
{
    int ghost_cells=3;
    // boundary conditions
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Before apply boundary",0,0);
    if(!projection_boundary) projection_boundary=boundary;
    projection_boundary->Apply_Boundary_Condition_Face(projection.p_grid,face_velocities,time+dt);
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost(projection.p_grid,ghost_cells);
    projection_boundary->Fill_Ghost_Faces(projection.p_grid,face_velocities,face_velocities_ghost,time,ghost_cells);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After apply boundary",0,0);

    assert(Consistent_Boundary_Conditions(face_velocities));

    // viscosity
    if(use_levelset_viscosity && implicit_viscosity && dt && viscosity!=0){
        LEVELSET_VISCOSITY_UNIFORM<TV> levelset_viscosity(bc_callbacks,grid,dt,projection.density,viscosity);
        levelset_viscosity.print_matrix=print_viscosity_matrix;
        levelset_viscosity.Apply_Viscosity(face_velocities,false,true,false);}
    else if(implicit_viscosity && dt && (use_variable_viscosity || viscosity!=0)){
        projection.Make_Divergence_Free(face_velocities,dt,time);
        Implicit_Viscous_Update(face_velocities,dt,time);}
        //VISCOSITY<TV> viscosity_helper(*projection.elliptic_solver,variable_viscosity,projection.density,viscosity,implicit_viscosity,use_explicit_part_of_implicit_viscosity,use_variable_viscosity,maximum_implicit_viscosity_iterations);
        //viscosity_helper.Add_Implicit_Forces_Before_Projection(grid,face_velocities,face_velocities,dt,time);}

    projection.Make_Divergence_Free(face_velocities,dt,time);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After final projection",0,0);
}
//#####################################################################
// Function Implicit_Viscous_Update
//#####################################################################
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Implicit_Viscous_Update(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
{ 
    int number_of_ghost_cells=3;
    ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost;face_velocities_ghost.Resize(grid,number_of_ghost_cells,false);
    boundary->Fill_Ghost_Faces(grid,face_velocities,face_velocities_ghost,time,number_of_ghost_cells);

    for(int axis=0;axis<TV::m;axis++){
        IMPLICIT_VISCOSITY_UNIFORM<TV> implicit_viscosity(*projection.elliptic_solver,variable_viscosity,projection.density,viscosity,0,axis,false,false);
        implicit_viscosity.Viscous_Update(grid,face_velocities,face_velocities_ghost,dt,time,maximum_implicit_viscosity_iterations);}
    if(mpi_grid) mpi_grid->Copy_Common_Face_Data(face_velocities);
}
//#####################################################################
// Function Real_CFL
//#####################################################################
template<class TV> int INCOMPRESSIBLE_UNIFORM<TV>::
Real_CFL(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const bool inviscid,const bool viscous_only,T input_dt) const
{
    T dt=CFL(face_velocities,inviscid,viscous_only);
    return (int) (input_dt/dt + 1);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR INCOMPRESSIBLE_UNIFORM<TV>::
CFL(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const bool inviscid,const bool viscous_only) const
{
    TV DX=grid.dX,sqr_DX=DX*DX,max_abs_V,one_over_DX=grid.one_over_dX;
    int ghost_cells=3;
    // convection
    T dt_convection=0;
    if(!viscous_only){
        for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            T local_V_norm=0;for(int axis=0;axis<TV::m;axis++)
                local_V_norm+=grid.one_over_dX[axis]*maxabs(face_velocities(axis,grid.First_Face_Index_In_Cell(axis,cell)),
                    face_velocities(axis,grid.Second_Face_Index_In_Cell(axis,cell)));
            dt_convection=max(dt_convection,local_V_norm);}}
    // surface tension
    T dt_surface_tension=0;
    if(nonzero_surface_tension){
        LEVELSET<TV>& levelset=*projection.collidable_solver->levelset;levelset.Compute_Curvature();
        ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells));levelset.boundary->Fill_Ghost_Cells(grid,levelset.phi,phi_ghost,0,0,ghost_cells); // TODO: use real dt, time
        T kappa_cfl=0;
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();
            T phi_1=phi_ghost(cell_1);T phi_2=phi_ghost(cell_2);
            if(LEVELSET_UTILITIES<T>::Interface(phi_1,phi_2)){
                T surface_tension_coefficient=surface_tension;
                if(use_variable_surface_tension) surface_tension_coefficient=(T).5*(variable_surface_tension(cell_1)+variable_surface_tension(cell_2));
                if(surface_tension_coefficient){
                    T curvature=LEVELSET_UTILITIES<T>::Average(phi_1,(*levelset.curvature)(cell_1),phi_2,(*levelset.curvature)(cell_2));
                    kappa_cfl=max(kappa_cfl,abs(curvature*surface_tension_coefficient/projection.density));}}}
        dt_surface_tension=sqrt(kappa_cfl)/grid.dX.Min();}
    T dt_viscosity=0;
    /*if(!inviscid && nonzero_viscosity){
        T norm_2_over_sqr_DX=2*Inverse(sqr_DX).Sum_Abs();
        if(viscosity) dt_viscosity=viscosity/projection.density*norm_2_over_sqr_DX;
        if(use_variable_viscosity){assert(!projection.flame);dt_viscosity=variable_viscosity.Max_Abs()/projection.density*norm_2_over_sqr_DX;}}*/
    if(viscous_only) return 1/max(dt_viscosity,1/max_time_step);
    TV max_force;
    if(use_force) max_force=force.Max_Abs();
    T dt_force=0;
    if(use_force) dt_force=(max_force*one_over_DX).Sum_Abs();
    dt_force+=(gravity*one_over_DX).Sum_Abs();
    if(strain) dt_force+=1/strain->CFL(projection.density);
    T dt_overall=(dt_convection+dt_viscosity+sqrt(sqr(dt_convection+dt_viscosity)+4*dt_force+4*sqr(dt_surface_tension)))/2; 
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function Extrapolate_Velocity_Across_Interface
//#####################################################################
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Extrapolate_Velocity_Across_Interface(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<T,TV_INT>& phi_ghost,const bool enforce_divergence_free,const T band_width,
    const T damping,const TV& air_speed,const ARRAY<VECTOR<bool,TV::m>,FACE_INDEX<TV::m> >* face_neighbors_visible,const ARRAY<bool,FACE_INDEX<TV::m> >* fixed_faces_input)
{
    T delta=band_width*grid.dX.Max();
    const int ghost_cells=2*(int)std::ceil(band_width)+1;
    if(mpi_grid){
        ARRAY<T,FACE_INDEX<TV::m> > phi_faces(grid,ghost_cells);
        ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost(grid,ghost_cells,false);boundary->Fill_Ghost_Faces(grid,face_velocities,face_velocities_ghost,0,ghost_cells); // TODO: use real time
        ARRAY<bool,FACE_INDEX<TV::m> > fixed_faces=fixed_faces_input?*fixed_faces_input:ARRAY<bool,FACE_INDEX<TV::m> >(grid,ghost_cells);
        for(int axis=0;axis<TV::m;axis++){
            T_ARRAYS_BASE &phi_face=phi_faces.Component(axis),&face_velocity=face_velocities_ghost.Component(axis);ARRAYS_ND_BASE<bool,TV_INT>& fixed_face=fixed_faces.Component(axis);
            for(FACE_ITERATOR<TV> iterator(grid,ghost_cells,GRID<TV>::INTERIOR_REGION,-1,axis);iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Face_Index();phi_face(index)=(T).5*(phi_ghost(iterator.First_Cell_Index())+phi_ghost(iterator.Second_Cell_Index()));
                if(phi_face(index)<=0) fixed_face(index)=true;if(phi_face(index) >= delta && !fixed_face(index)) face_velocity(index)=(T)0;}}
        //mpi_grid->Exchange_Boundary_Face_Data(fixed_faces);
        for(int axis=0;axis<TV::m;axis++){
            T_ARRAYS_BASE &phi_face=phi_faces.Component(axis),&face_velocity=face_velocities_ghost.Component(axis);ARRAYS_ND_BASE<bool,TV_INT>& fixed_face=fixed_faces.Component(axis);
            GRID<TV> face_grid=grid.Get_Face_Grid(axis);T_EXTRAPOLATION_SCALAR extrapolate(face_grid,phi_face,face_velocity,ghost_cells);
            extrapolate.Set_Band_Width(band_width);extrapolate.Set_Custom_Seed_Done(&fixed_face);
            if(face_neighbors_visible) extrapolate.Set_Collision_Aware_Extrapolation(face_neighbors_visible->Component(axis));
            extrapolate.Extrapolate(0,false);
            ARRAY<T,TV_INT>::Get(face_velocities.Component(axis),face_velocity);
            if(damping) for(FACE_ITERATOR<TV> iterator(grid,0,GRID<TV>::INTERIOR_REGION,ghost_cells,axis);iterator.Valid();iterator.Next()){TV_INT index=iterator.Face_Index();
                if(!fixed_face(index) && phi_face(index)<delta) face_velocity(index)=(1-damping)*face_velocity(index)+damping*air_speed[axis];}}}
    else{
        for(int axis=0;axis<TV::m;axis++){
            GRID<TV> face_grid=grid.Get_Face_Grid(axis);ARRAY<T,TV_INT> phi_face(face_grid.Domain_Indices(),false);T_ARRAYS_BASE& face_velocity=face_velocities.Component(axis);
            ARRAY<bool,TV_INT> fixed_face;
            if(fixed_faces_input) fixed_face=fixed_faces_input->Component(axis); else fixed_face=ARRAY<bool,TV_INT>(face_grid.Domain_Indices());
            for(FACE_ITERATOR<TV> iterator(grid,0,GRID<TV>::WHOLE_REGION,-1,axis);iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Face_Index();phi_face(index)=(T).5*(phi_ghost(iterator.First_Cell_Index())+phi_ghost(iterator.Second_Cell_Index()));
                if(phi_face(index)<=0) fixed_face(index)=true;if(phi_face(index) >= delta && !fixed_face(index)) face_velocity(index)=(T)0;}
            LOG::cout<<"something..."<<std::endl;  // TODO(jontg): If this log statement doesn't appear, the code crashes in release mode...
            T_EXTRAPOLATION_SCALAR extrapolate(face_grid,phi_face,face_velocity,ghost_cells);extrapolate.Set_Band_Width(band_width);extrapolate.Set_Custom_Seed_Done(&fixed_face);
            if(face_neighbors_visible) extrapolate.Set_Collision_Aware_Extrapolation(face_neighbors_visible->Component(axis));
            extrapolate.Extrapolate();
            if(damping) for(FACE_ITERATOR<TV> iterator(grid,0,GRID<TV>::WHOLE_REGION,-1,axis);iterator.Valid();iterator.Next()){TV_INT index=iterator.Face_Index();
                if(!fixed_face(index) && phi_face(index)<delta) face_velocity(index)=(1-damping)*face_velocity(index)+damping*air_speed[axis];}}}

    // make extrapolated velocity divergence free
    if(enforce_divergence_free){
        ARRAY<T,TV_INT> p_new(grid.Domain_Indices(1));ARRAY<bool,TV_INT> psi_D_new(grid.Domain_Indices(1));ARRAY<bool,FACE_INDEX<TV::m> > psi_N_new(grid);
        ARRAY<T,TV_INT>::Exchange(p_new,projection.p);ARRAY<bool,TV_INT>::Exchange(psi_D_new,projection.elliptic_solver->psi_D);
        ARRAY<bool,FACE_INDEX<TV::m> >::Exchange(psi_N_new,projection.elliptic_solver->psi_N);
        projection.elliptic_solver->Set_Dirichlet_Outer_Boundaries();
        for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
            if(phi_ghost(index) >= delta) projection.elliptic_solver->psi_D(index)=true;
            else if(phi_ghost(index) <= 0) projection.elliptic_solver->psi_N.Set_All_Faces(true,index);
            else{
                bool local_maximum=true;
                for(int i=0;i<GRID<TV>::number_of_neighbors_per_cell;i++) if(phi_ghost(index)<phi_ghost(iterator.Cell_Neighbor(i))){local_maximum=false;break;}
                if(local_maximum)projection.elliptic_solver->psi_D(index)=true;}}
        projection.Make_Divergence_Free(face_velocities,0,0); // TODO: use real dt/time
        ARRAY<T,TV_INT>::Exchange(p_new,projection.p);ARRAY<bool,TV_INT>::Exchange(psi_D_new,projection.elliptic_solver->psi_D);
        ARRAY<bool,FACE_INDEX<TV::m> >::Exchange(psi_N_new,projection.elliptic_solver->psi_N);} // restore pressure for use as initial guess for incompressible projection
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Set_Dirichlet_Boundary_Conditions(const ARRAY<T,TV_INT>* phi,const T pressure)
{
    ARRAY<bool,TV_INT>& psi_D=projection.elliptic_solver->psi_D;
    if(phi) for(CELL_ITERATOR<TV> iterator(projection.p_grid);iterator.Valid();iterator.Next()) if((*phi)(iterator.Cell_Index())>0){
        psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=pressure;}
    if(projection.elliptic_solver->mpi_grid){
        projection.elliptic_solver->mpi_grid->Exchange_Boundary_Cell_Data(psi_D,1,false);
        projection.elliptic_solver->mpi_grid->Exchange_Boundary_Cell_Data(projection.p,1,false);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After exchange dirichlet",0,0);
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Set_Dirichlet_Boundary_Conditions(const ARRAY<T,TV_INT>* phi,const ARRAY<T,TV_INT>& pressure)
{
    ARRAY<bool,TV_INT>& psi_D=projection.elliptic_solver->psi_D;
    if(phi) for(CELL_ITERATOR<TV> iterator(projection.p_grid);iterator.Valid();iterator.Next()) if((*phi)(iterator.Cell_Index())>0){
        psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=pressure(iterator.Cell_Index());}
    if(mpi_grid){
        mpi_grid->Exchange_Boundary_Cell_Data(psi_D,1,false);
        mpi_grid->Exchange_Boundary_Cell_Data(projection.p,1,false);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After exchange dirichlet",0,0);
}
//#####################################################################
// Function Add_Surface_Tension
//#####################################################################
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Add_Surface_Tension(LEVELSET<TV>& levelset,const T time)
{
    
    LAPLACE_UNIFORM<TV>& elliptic_solver=*projection.elliptic_solver;GRID<TV>& p_grid=elliptic_solver.grid;
    ARRAY<T,TV_INT>& phi=levelset.phi;
    LINEAR_INTERPOLATION_UNIFORM<TV,T> interpolation;

    if(projection.collidable_solver->second_order_cut_cell_method) for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
        if(!projection.elliptic_solver->psi_N.Component(iterator.Axis())(iterator.Face_Index()) && LEVELSET_UTILITIES<T>::Interface(phi(iterator.First_Cell_Index()),phi(iterator.Second_Cell_Index()))){
            T theta=LEVELSET_UTILITIES<T>::Theta(phi(iterator.First_Cell_Index()),phi(iterator.Second_Cell_Index()));
            TV location=theta*(grid.Center(iterator.Second_Cell_Index())-grid.Center(iterator.First_Cell_Index()))+grid.Center(iterator.First_Cell_Index());
            T curvature_at_interface=levelset.Compute_Curvature(location);
            T surface_tension_coefficient=surface_tension;
            if(use_variable_surface_tension)surface_tension_coefficient=interpolation.Clamped_To_Array(grid,variable_surface_tension,location);
            projection.collidable_solver->u_interface.Component(iterator.Axis())(iterator.Face_Index())=-surface_tension_coefficient*curvature_at_interface;}}
    else{
        for(CELL_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()) if(elliptic_solver.psi_D(iterator.Cell_Index()) && phi(iterator.Cell_Index()) < 5*grid.dX.Max()){
            T surface_tension_coefficient=surface_tension;
            if(use_variable_surface_tension) surface_tension_coefficient=variable_surface_tension(iterator.Cell_Index());
            projection.p(iterator.Cell_Index())=-surface_tension_coefficient*levelset.Compute_Curvature(phi,iterator.Cell_Index());}}
}
//#####################################################################
// Function Apply_Vorticity_Confinement_Force
//#####################################################################
template<class TV,class T_ARRAYS_TV,class T> static void
Apply_Vorticity_Confinement_Force_Helper(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,T_ARRAYS_TV& F,const INCOMPRESSIBLE_UNIFORM<TV>& incompressible)
{
    typedef VECTOR<int,TV::m> TV_INT;
    // want cells to face averaging here
    if(incompressible.collision_body_list){
        AVERAGING_COLLIDABLE_UNIFORM<TV,FACE_LOOKUP_COLLIDABLE_UNIFORM<TV> > vorticity_averaging_collidable(*incompressible.collision_body_list,T());
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            face_velocities.Component(axis)(iterator.Face_Index())+=vorticity_averaging_collidable.Cell_To_Face(grid,axis,iterator.Face_Index(),F);}}
    else
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
            face_velocities.Component(axis)(iterator.Face_Index())+=(incompressible.use_vorticity_weights?incompressible.vorticity_weights(iterator.Full_Index()):1)*(F(iterator.First_Cell_Index())[axis]+F(iterator.Second_Cell_Index())[axis]);}
}
static void
Apply_Vorticity_Confinement_Force_Helper(const GRID<VECTOR<float,1> >&,ARRAY<float,FACE_INDEX<1> >&,ARRAY<VECTOR<float,1> ,VECTOR<int,1> >&,const INCOMPRESSIBLE_UNIFORM<VECTOR<float,1> >&)
{PHYSBAM_NOT_IMPLEMENTED();}
static void
Apply_Vorticity_Confinement_Force_Helper(const GRID<VECTOR<double,1> >&,ARRAY<double,FACE_INDEX<1> >&,ARRAY<VECTOR<double,1> ,VECTOR<int,1> >&,const INCOMPRESSIBLE_UNIFORM<VECTOR<double,1> >&)
{PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Apply_Vorticity_Confinement_Force(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<TV,TV_INT>& F)
{
    Apply_Vorticity_Confinement_Force_Helper(grid,face_velocities,F,*this);
}
//#####################################################################
// Function Compute_Vorticity_Confinement_Force_Helper
//#####################################################################
template<class TV,class T,class TV_INT> static void 
Compute_Vorticity_Confinement_Force_Helper(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<TV,TV_INT>& F,const INCOMPRESSIBLE_UNIFORM<TV>& incompressible)
{
    typedef ARRAY<typename TV::SPIN,TV_INT> T_ARRAYS_SPIN;typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;
    T_ARRAYS_SPIN vorticity(grid.Cell_Indices(2),false);
    ARRAY<T,TV_INT> vorticity_magnitude(grid.Cell_Indices(2));
    if(incompressible.collision_body_list){
        FACE_LOOKUP_UNIFORM<TV> face_velocities_lookup_uniform(face_velocities_ghost);
        FACE_LOOKUP_COLLIDABLE_UNIFORM<TV> face_velocities_lookup(face_velocities_lookup_uniform,*incompressible.collision_body_list,&incompressible.valid_mask);
        VORTICITY_UNIFORM<TV>::Vorticity(grid,face_velocities_lookup,vorticity,vorticity_magnitude);}
    else VORTICITY_UNIFORM<TV>::Vorticity(grid,T_FACE_LOOKUP(face_velocities_ghost),vorticity,vorticity_magnitude);
    for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){ // do collision awareness when these are averaged to faces
        TV vortex_normal_vector=LEVELSET<TV>::Normal_At_Node(grid,vorticity_magnitude,iterator.Cell_Index());
        F(iterator.Cell_Index())=TV::Cross_Product(vortex_normal_vector,vorticity(iterator.Cell_Index()));}

    if(incompressible.vc_projection_direction){
        for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()) F(iterator.Cell_Index())(incompressible.vc_projection_direction)=0;}
    if(incompressible.momentum_conserving_vorticity){
        for(int axis=0;axis<TV::m;axis++){T sum=0; int count=0;
            for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
                if(incompressible.variable_vorticity_confinement(cell)!=(T)0){sum+=F(cell)[axis];count++;}}
            PHYSBAM_ASSERT(count>0);
            sum/=count;T final_sum=0;
            for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index(); 
                if(incompressible.variable_vorticity_confinement(cell)!=(T)0){F(cell)[axis]-=sum;final_sum+=F(cell)[axis];}}}}
}
//#####################################################################
// Function Compute_Vorticity_Confinement_Force_Helper
//#####################################################################
template<class T,class TV_INT> static void
Compute_Vorticity_Confinement_Force_Helper(const GRID<VECTOR<T,1> >&,const ARRAY<T,FACE_INDEX<1> >&,ARRAY<VECTOR<T,1>,TV_INT>&,const INCOMPRESSIBLE_UNIFORM<VECTOR<T,1> >&)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Compute_Vorticity_Confinement_Force(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<TV,TV_INT>& F)
{
    Compute_Vorticity_Confinement_Force_Helper(grid,face_velocities_ghost,F,*this);
}
//#####################################################################
// Function Use_Maccormack_Advection
//#####################################################################
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Use_Maccormack_Advection(const ARRAY<bool,TV_INT>* node_mask, const ARRAY<bool,TV_INT>* cell_mask, const ARRAY<bool,FACE_INDEX<TV::m> >* face_mask)
{
    advection_maccormack=new ADVECTION_MACCORMACK_UNIFORM<TV,T,ADVECTION<TV,T> >(*advection,node_mask,cell_mask,face_mask,thread_queue);
    Set_Custom_Advection(*advection_maccormack);
}
//#####################################################################
// Function Consistent_Boundary_Conditions
//#####################################################################
template<class TV> bool INCOMPRESSIBLE_UNIFORM<TV>::
Consistent_Boundary_Conditions(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities) const
{
    if(projection.elliptic_solver->mpi_grid){
        LOG::cout<<"checking for consistent mpi boundaries"<<std::endl;
        ARRAY<bool,TV_INT> psi_D_ghost(projection.elliptic_solver->psi_D);
        ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost(face_velocities);ARRAY<T,FACE_INDEX<TV::m> > psi_N_ghost(projection.p_grid);
        projection.elliptic_solver->mpi_grid->Exchange_Boundary_Cell_Data(psi_D_ghost,1);
        for(int axis=0;axis<TV::m;axis++)for(int axis_side=0;axis_side<2;axis_side++){int side=2*axis+axis_side;
            if(projection.elliptic_solver->mpi_grid->Neighbor(axis,axis_side)){TV_INT exterior_cell_offset=axis_side==0?-TV_INT::Axis_Vector(axis):TV_INT();
                for(FACE_ITERATOR<TV> iterator(projection.p_grid,0,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                    TV_INT face=iterator.Face_Index(),cell=face+exterior_cell_offset;int axis=iterator.Axis();
                    psi_N_ghost(axis,face)=(T)projection.elliptic_solver->psi_N(axis,face);
                    face_velocities_ghost(axis,face)=face_velocities(axis,face);
                    assert(projection.elliptic_solver->psi_D(cell)==psi_D_ghost(cell));}}}
        projection.elliptic_solver->mpi_grid->Assert_Common_Face_Data(face_velocities_ghost);projection.elliptic_solver->mpi_grid->Assert_Common_Face_Data(psi_N_ghost);}
    return true;
}
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Initialize_Grids(const GRID<TV>& grid_input)
{
    INCOMPRESSIBLE<TV>::Initialize_Grids(grid_input);
    grid=grid_input.Get_MAC_Grid();projection.Initialize_Grid(grid);
    if(strain) strain->Initialize_Grid(grid);
}
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Set_Body_Force(const bool use_force_input)
{
    use_force=use_force_input;
    if(use_force) force.Resize(grid,1);
    else force.Clean_Memory();
}
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Set_Variable_Surface_Tension(const bool use_variable_surface_tension_input)
{
    use_variable_surface_tension=use_variable_surface_tension_input;
    if(use_variable_surface_tension){
        nonzero_surface_tension=true;
        variable_surface_tension.Resize(grid.Cell_Indices(1));
        surface_tension=0;}
    else variable_surface_tension.Clean_Memory();
}
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Set_Variable_Viscosity(const bool use_variable_viscosity_input)
{
    use_variable_viscosity=use_variable_viscosity_input;
    if(use_variable_viscosity){
        nonzero_viscosity=true;
        variable_viscosity.Resize(grid.Cell_Indices(1));
        viscosity=0;}
    else variable_viscosity.Clean_Memory();
}
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Use_Variable_Vorticity_Confinement(GRID<TV>& grid,const bool use_variable_vorticity_confinement_input)
{
    use_variable_vorticity_confinement=use_variable_vorticity_confinement_input;
    if(use_variable_vorticity_confinement) variable_vorticity_confinement.Resize(grid.Cell_Indices(1));
    else variable_vorticity_confinement.Clean_Memory();
}
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Use_Variable_Vorticity_Confinement(const bool use_variable_vorticity_confinement_input)
{
    Use_Variable_Vorticity_Confinement(grid,use_variable_vorticity_confinement_input);
}
template<class TV> void INCOMPRESSIBLE_UNIFORM<TV>::
Use_Strain()
{
    delete strain;
    strain=new FLUID_STRAIN_UNIFORM<TV>(grid);
}
//#####################################################################
namespace PhysBAM{
template class INCOMPRESSIBLE_UNIFORM<VECTOR<float,1> >;
template class INCOMPRESSIBLE_UNIFORM<VECTOR<float,2> >;
template class INCOMPRESSIBLE_UNIFORM<VECTOR<float,3> >;
template class INCOMPRESSIBLE_UNIFORM<VECTOR<double,1> >;
template class INCOMPRESSIBLE_UNIFORM<VECTOR<double,2> >;
template class INCOMPRESSIBLE_UNIFORM<VECTOR<double,3> >;
}
