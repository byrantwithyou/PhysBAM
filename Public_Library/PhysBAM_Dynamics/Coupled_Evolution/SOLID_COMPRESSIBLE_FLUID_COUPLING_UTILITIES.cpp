//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry_Computations/BOX_BOX_INTERSECTION_AREA.h>
#include <PhysBAM_Geometry/Basic_Geometry_Computations/BOX_POLYGON_INTERSECTION_AREA.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/OBJECTS_IN_CELL.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/CUT_CELL_COMPUTATIONS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/RIGID_GEOMETRY_RASTERIZATION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
#include <PhysBAM_Dynamics/Forces_And_Torques/EULER_FLUID_FORCES.h>
#include <iomanip>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES(EULER_UNIFORM<T_GRID>& euler_input,MPI_UNIFORM_GRID<T_GRID>* mpi_grid_input):
    euler(euler_input),mpi_grid(mpi_grid_input),collision_bodies_affecting_fluid(0),thinshell(true),use_fast_marching(false),use_higher_order_solid_extrapolation(true),
    fluid_affects_solid(false),number_of_cells_to_extrapolate(7),solid_state(TV_DIMENSION()),euler_fluid_forces(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
~SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES()
{
    delete euler_fluid_forces;
}
//#####################################################################
// Initialize_Solid_Fluid_Coupling
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Initialize_Solid_Fluid_Coupling(GRID_BASED_COLLISION_GEOMETRY<T_GRID>* collision_bodies_affecting_fluid_input)
{
    phi_all_solids_negated.Resize(euler.grid.Domain_Indices(number_of_cells_to_extrapolate));
    outside_fluid.Resize(euler.grid.Domain_Indices(number_of_cells_to_extrapolate));
    uncovered_cells.Resize(euler.grid.Domain_Indices(1));
    collision_bodies_affecting_fluid=dynamic_cast<GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>*>(collision_bodies_affecting_fluid_input);
    collision_bodies_affecting_fluid->Initialize_Grids();
    collision_bodies_affecting_fluid->Rasterize_Objects();
    if(fluid_affects_solid){
        if(!euler.timesplit){
            euler.compute_pressure_fluxes=true;
            solid_fluid_face_time_n.Resize(euler.grid,0);
            cells_inside_fluid_time_n.Resize(euler.grid.Domain_Indices(0));
            pressure_at_faces.Resize(euler.grid.Domain_Indices(0));}
        if(euler.timesplit && thinshell){
            U_n.Resize(euler.grid.Domain_Indices(3));
            accumulated_flux.Resize(euler.grid.Domain_Indices());
            near_interface.Resize(euler.grid.Domain_Indices(1));near_interface.Fill(false);
            advection_velocities_n.Resize(euler.grid.Domain_Indices());
            cut_cells_n.Resize(euler.grid.Domain_Indices(1));cut_cells_n.Fill(0);
            cut_cells_np1.Resize(euler.grid.Domain_Indices(1));cut_cells_np1.Fill(0);
            cell_volumes_n.Resize(euler.grid.Domain_Indices(2));
            cell_volumes_np1.Resize(euler.grid.Domain_Indices(2));
            psi_n.Resize(euler.grid.Domain_Indices());
            psi_np1.Resize(euler.grid.Domain_Indices());

            uncovered_cells_n_p_half.Resize(euler.grid.Domain_Indices(1));
            cut_cells_n_p_half.Resize(euler.grid.Domain_Indices(1));
            cell_volumes_n_p_half.Resize(euler.grid.Domain_Indices(2));
            psi_n_p_half.Resize(euler.grid.Domain_Indices());}}
}
//#####################################################################
// Update_Cut_Out_Grid
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Update_Cut_Out_Grid()
{   
    Compute_Phi_Solids(0);
    if(euler.timesplit && !thinshell)
        for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            //euler.psi(cell_index)=!outside_fluid(cell_index);}
            euler.psi(cell_index)=phi_all_solids_negated(cell_index)<(T)2*euler.grid.dX.Max();}
    else
        for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            euler.psi(cell_index)=!outside_fluid(cell_index);}
}   
//#####################################################################
// Function Get_Neumann_Data 
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Get_Neumann_Data(const TV& location,const T max_distance,TV& normal_direction,T& object_velocity_normal_component,TV& reflected_point) const
{
    assert(!euler.timesplit || !thinshell);
    T distance;TV boundary_point;COLLISION_GEOMETRY_ID body_id;int simplex_id;
    boundary_point=collision_bodies_affecting_fluid->collision_geometry_collection.Closest_Boundary_Point(location,max_distance,distance,body_id,simplex_id);
    const COLLISION_GEOMETRY<TV>& collision_body=collision_bodies_affecting_fluid->collision_geometry_collection(body_id);
    reflected_point=location+((T)2*(boundary_point-location));
    normal_direction=boundary_point-location;normal_direction.Normalize();
    TV object_velocity=collision_body.Pointwise_Object_Velocity(simplex_id,boundary_point);
    object_velocity_normal_component=TV::Dot_Product(object_velocity,normal_direction);
}
//#####################################################################
// Function Get_Neumann_Data 
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Get_Neumann_Data(const TV& location,const T max_distance,TV& object_velocity,TV& reflected_point) const
{
    assert(!euler.timesplit || !thinshell);
    T distance;TV boundary_point;COLLISION_GEOMETRY_ID body_id;int simplex_id;
    boundary_point=collision_bodies_affecting_fluid->collision_geometry_collection.Closest_Boundary_Point(location,max_distance,distance,body_id,simplex_id);
    const COLLISION_GEOMETRY<TV>& collision_body=collision_bodies_affecting_fluid->collision_geometry_collection(body_id);
    reflected_point=location+((T)2*(boundary_point-location));
    object_velocity=collision_body.Pointwise_Object_Velocity(simplex_id,boundary_point);
}
//#####################################################################
// Function Fill_Uncovered_Cells 
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Fill_Uncovered_Cells()
{
    if(!euler.timesplit || !thinshell){
        T_LINEAR_INTERPOLATION_DIMENSION interpolation;
        T max_distance,object_velocity_normal_component;TV location,normal_direction,object_velocity,reflected_point;
        for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            if(uncovered_cells(cell_index)){location=euler.grid.Center(cell_index);
                max_distance=phi_all_solids_negated(cell_index)*(T)2;
                if(euler.use_solid_velocity_in_ghost_cells){Get_Neumann_Data(location,max_distance,object_velocity,reflected_point);normal_direction=object_velocity;}  // TODO(jontg): HACK
                else Get_Neumann_Data(location,max_distance,normal_direction,object_velocity_normal_component,reflected_point);
                euler.U(cell_index)=interpolation.Clamped_To_Array(euler.grid,euler.U,reflected_point);
                euler.conservation->object_boundary->Apply_Neumann_Boundary_Condition(euler.U(cell_index),normal_direction,object_velocity_normal_component);}}
        euler.Invalidate_Ghost_Cells();}
}
//#####################################################################
// Function Extrapolate_State_Into_Solids
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Extrapolate_State_Into_Solids(T_ARRAYS_SCALAR& phi_all_solids_negated,const int number_of_ghost_cells,const int number_of_cells_to_extrapolate)
{
    assert(!euler.timesplit || !thinshell);
    T_ARRAYS_DIMENSION_SCALAR U_extrapolated(euler.grid.Domain_Indices(number_of_ghost_cells));

    if(mpi_grid){
        BOUNDARY_UNIFORM<GRID<TV>,TV_DIMENSION> boundary; // Constant extrapolation at non-mpi faces.
        BOUNDARY_MPI<GRID<TV>,TV_DIMENSION> mpi_boundary(mpi_grid,boundary);
        mpi_boundary.Fill_Ghost_Cells(euler.grid,euler.U,U_extrapolated,(T)0,(T)0,number_of_ghost_cells);}
    else U_extrapolated=euler.U;

    T_EXTRAPOLATION_SCALAR_DIMENSION extrapolate(euler.grid,phi_all_solids_negated,U_extrapolated,number_of_ghost_cells);extrapolate.Set_Band_Width((T)number_of_cells_to_extrapolate);
    extrapolate.Extrapolate((T)0,false);
    T_ARRAYS_DIMENSION_SCALAR::Get(euler.U,U_extrapolated);

    T band_width=number_of_cells_to_extrapolate*euler.grid.dX.Max();
    T_LINEAR_INTERPOLATION_DIMENSION interpolation;
    T max_distance,object_velocity_normal_component;TV location,normal_direction,object_velocity,reflected_point;
    const RANGE<TV>& domain=RANGE<TV>::Intersect(euler.grid.Ghost_Domain(number_of_ghost_cells),mpi_grid?mpi_grid->global_grid.Domain():euler.grid.Domain());

    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();TV location=iterator.Location();
        if(outside_fluid(cell_index) && phi_all_solids_negated(cell_index)<band_width){
            max_distance=phi_all_solids_negated(cell_index)*(T)2;
            if(euler.use_solid_velocity_in_ghost_cells){Get_Neumann_Data(location,max_distance,object_velocity,reflected_point);normal_direction=object_velocity;}  // TODO(jontg): HACK
            else Get_Neumann_Data(location,max_distance,normal_direction,object_velocity_normal_component,reflected_point);
            COLLISION_GEOMETRY_ID body_id;
            if(use_higher_order_solid_extrapolation && domain.Inside(reflected_point,euler.grid.dX.Max()*(T).5) && 
                    !collision_bodies_affecting_fluid->Implicit_Geometry_Lazy_Inside_Any_Body(reflected_point,body_id)){
                euler.U(cell_index)=interpolation.Clamped_To_Array(euler.grid,U_extrapolated,reflected_point);}
            euler.conservation->object_boundary->Apply_Neumann_Boundary_Condition(euler.U(cell_index),normal_direction,object_velocity_normal_component);}}
    euler.Invalidate_Ghost_Cells();
}
//#####################################################################
// Function Compute_Phi_Solids
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Compute_Phi_Solids(const int number_of_ghost_cells)
{
    ARRAY<TV_INT> seed_indices;
    phi_all_solids_negated.Fill(-FLT_MAX);
    outside_fluid.Fill(false);
    T_FACE_ARRAYS_BOOL kinematic_faces(euler.grid.Domain_Indices(1));kinematic_faces.Fill(false);

    for(COLLISION_GEOMETRY_ID id(0);id<collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m;id++)
        if(collision_bodies_affecting_fluid->collision_geometry_collection.Is_Active(id)){
            COLLISION_GEOMETRY<TV>& collision_body=*(collision_bodies_affecting_fluid->collision_geometry_collection.bodies(id));
            T collision_thickness_over_two=(T).5*collision_bodies_affecting_fluid->collision_thickness;
            for(CELL_ITERATOR iterator(euler.grid,euler.grid.Clamp_To_Cell(collision_body.Axis_Aligned_Bounding_Box().Thickened(euler.grid.dX.Max()*(T)2),number_of_cells_to_extrapolate));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();T phi_value=collision_body.Implicit_Geometry_Extended_Value(iterator.Location());
                if(collision_body.Inside(iterator.Location(),collision_thickness_over_two)) outside_fluid(cell_index)=true;
                phi_all_solids_negated(cell_index)=max(-phi_value,phi_all_solids_negated(cell_index));}}

    if(use_fast_marching){
        collision_bodies_affecting_fluid->outside_fluid=&outside_fluid;
        UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV> iterator_info(*collision_bodies_affecting_fluid);
        iterator_info.Initialize_Collision_Aware_Face_Iterator(outside_fluid,kinematic_faces,7,false);

        for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV> iterator(iterator_info);iterator.Valid();iterator.Next()){
            seed_indices.Append(iterator.First_Cell_Index());seed_indices.Append(iterator.Second_Cell_Index());}

        typename LEVELSET_POLICY<GRID<TV> >::LEVELSET levelset(euler.grid,phi_all_solids_negated); // TODO(jontg): Make this a permanent member variable?
        FAST_MARCHING_METHOD_UNIFORM<GRID<TV> > fmm(levelset,number_of_cells_to_extrapolate);
        fmm.Fast_Marching_Method(phi_all_solids_negated,euler.grid.dX.Max()*(T)number_of_cells_to_extrapolate,&seed_indices,true);}
}
//#####################################################################
// Function Fill_Solid_Cells
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Fill_Solid_Cells(bool fill_pressure_only)
{
    if(euler.timesplit && thinshell) return;
    if(collision_bodies_affecting_fluid){
        if(thinshell) Fill_Uncovered_Cells();
        int number_of_ghost_cells=mpi_grid?3:0;
        Compute_Phi_Solids(number_of_ghost_cells);
        if(fill_pressure_only){LOG::cout<<"Defunct code"<<std::endl;exit(-1);}
        else Extrapolate_State_Into_Solids(phi_all_solids_negated,number_of_ghost_cells,number_of_cells_to_extrapolate);}
}
//#####################################################################
// Function Project_Fluid_Pressure_At_Neumann_Faces
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Project_Fluid_Pressure_At_Neumann_Faces(const T_ARRAYS_SCALAR& p_ghost,T_FACE_ARRAYS_SCALAR& p_face) const
{
    // Bp
    const RANGE<TV>& domain=euler.grid.domain;
    for(FACE_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        TV_INT face_index=iterator.Face_Index(),first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
        if(euler.euler_projection.elliptic_solver->psi_N.Component(axis)(face_index)){
            int direction;
            // TODO: This doesn't work with thin-shells
            bool first_cell_inside_solid=euler.euler_projection.elliptic_solver->psi_D(iterator.First_Cell_Index()) ||
                !domain.Lazy_Inside_Half_Open(iterator.First_Cell_Center());
            bool second_cell_inside_solid=euler.euler_projection.elliptic_solver->psi_D(iterator.Second_Cell_Index()) ||
                !domain.Lazy_Inside_Half_Open(iterator.Second_Cell_Center());

            if(!first_cell_inside_solid && second_cell_inside_solid) direction=1; // solid on right
            else if(first_cell_inside_solid && !second_cell_inside_solid) direction=-1; // solid on left
            else if(!first_cell_inside_solid && !second_cell_inside_solid){
                direction=0;
                LOG::cerr<<"Warning: Neumann face has fluid on both sides, face_index="<<face_index<<", axis="<<axis<<std::endl;}
            else direction=0; // solid on both sides
    
            if(direction!=0){
                if(direction==1) p_face.Component(axis)(face_index)=p_ghost(first_cell_index);
                else p_face.Component(axis)(face_index)=p_ghost(second_cell_index);}}}
}
//#####################################################################
// Function Apply_Isobaric_Fix
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Apply_Isobaric_Fix(const T dt,const T time)
{
    euler.Fill_Ghost_Cells(dt,time,3);
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(euler.psi(cell_index) && phi_all_solids_negated(cell_index)<0){
            bool encountered_neumann_face=false;TV_INT reference_point;
            for(int axis=0;axis<T_GRID::dimension;axis++){
                TV_INT first_face_index=iterator.First_Face_Index(axis),second_face_index=iterator.Second_Face_Index(axis);
                if(euler.euler_projection.elliptic_solver->psi_N.Component(axis)(first_face_index) &&
                        euler.euler_projection.elliptic_solver->psi_N.Component(axis)(second_face_index)) continue; // cant get a reference point
                else if(euler.euler_projection.elliptic_solver->psi_N.Component(axis)(first_face_index)){
                    encountered_neumann_face=true;reference_point=iterator.Cell_Neighbor(2*axis+1);}
                else if(euler.euler_projection.elliptic_solver->psi_N.Component(axis)(second_face_index)){
                    encountered_neumann_face=true;reference_point=iterator.Cell_Neighbor(2*axis);}}
            if(encountered_neumann_face){
                LOG::cout<<"ISOBARIC FIX: fixing cell "<<cell_index<<" with reference cell "<<reference_point<<std::endl;
                T rho=euler.U(cell_index)(1);
                TV velocity=EULER<T_GRID>::Get_Velocity(euler.U,cell_index);
                T e=EULER<T_GRID>::e(euler.U,cell_index);
                T p_cell=euler.eos->p(rho,e);

                T rho_reference=euler.U_ghost(reference_point)(1);
                T e_reference=EULER<T_GRID>::e(euler.U_ghost,reference_point);
                T p_reference=euler.eos->p(rho_reference,e_reference);

                if(p_cell>p_reference) rho=rho_reference*sqrt(p_cell/p_reference); //isobaric fix
                else rho=rho_reference;
                //rho=rho_reference*pow(p_cell/p_reference,(T)(1/1.4)); //isobaric fix (constant entropy)
                e=euler.eos->e_From_p_And_rho(p_cell,rho);
                EULER<T_GRID>::Set_Euler_State_From_rho_velocity_And_internal_energy(euler.U,cell_index,rho,velocity,e);}}}
    euler.Invalidate_Ghost_Cells();
}
//#####################################################################
// Function Extract_Time_N_Data_For_Explicit_Fluid_Forces
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Extract_Time_N_Data_For_Explicit_Fluid_Forces()
{
    if(!fluid_affects_solid || euler.timesplit) return;
    T_ARRAYS_SCALAR p_approx(euler.grid.Domain_Indices(1));
    for(CELL_ITERATOR iterator(euler.grid,1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        p_approx(cell_index)=euler.eos->p(euler.U_ghost(cell_index)(1),euler.e(euler.U_ghost,cell_index));}
    euler.euler_projection.Compute_Face_Pressure_From_Cell_Pressures(euler.grid,euler.U_ghost,euler.psi,pressure_at_faces,p_approx);

    solid_fluid_face_time_n.Fill(false);
    collision_bodies_affecting_fluid->Compute_Psi_N(solid_fluid_face_time_n,0);
    cells_inside_fluid_time_n=euler.psi;
}
//#####################################################################
// Function Snapshot_State
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Snapshot_State(const T_ARRAYS_DIMENSION_SCALAR& U_ghost)
{
    T_ARRAYS_DIMENSION_SCALAR::Put(U_ghost,U_n);
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next())
        if(!euler.psi(iterator.Cell_Index())) U_n(iterator.Cell_Index())=TV_DIMENSION();
        else advection_velocities_n(iterator.Cell_Index())=EULER<T_GRID>::Get_Velocity(U_n(iterator.Cell_Index()));

    accumulated_flux.Fill(TV_DIMENSION());
    TV_DIMENSION accumulated_material;
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next())
        accumulated_material+=U_n(iterator.Cell_Index())*cell_volumes_np1(iterator.Cell_Index());
    LOG::cout<<std::setprecision(16)<<"TOTAL ACCUMULATED MATERIAL = "<<accumulated_material<<std::endl;

    near_interface.Fill(false);
}
//#####################################################################
// Function Initialize_Collision_Data
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Initialize_Collision_Data()
{
    CUT_CELL_COMPUTATIONS::Compute_Cut_Geometries(euler.grid,1,*collision_bodies_affecting_fluid,cut_cells_np1);

    psi_np1.Fill(true);
    for(COLLISION_GEOMETRY_ID id(0);id<collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m;id++)
        if(collision_bodies_affecting_fluid->collision_geometry_collection.Is_Active(id)){
            COLLISION_GEOMETRY<TV>& collision_body=*(collision_bodies_affecting_fluid->collision_geometry_collection.bodies(id));
            T collision_thickness_over_two=(T).5*collision_bodies_affecting_fluid->collision_thickness;
            for(CELL_ITERATOR iterator(euler.grid,euler.grid.Clamp_To_Cell(collision_body.Axis_Aligned_Bounding_Box().Thickened(euler.grid.dX.Max()*(T)2),0));iterator.Valid();iterator.Next()){
                if(collision_body.Inside(iterator.Location(),collision_thickness_over_two)) psi_np1(iterator.Cell_Index())=false;}}

    cell_volumes_np1.Fill(euler.grid.Cell_Size());
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(!psi_np1(cell_index) || (cut_cells_np1(cell_index) && !cut_cells_np1(cell_index)->dominant_element)) cell_volumes_np1(cell_index)=0;
        if(cut_cells_np1(cell_index))
            for(int poly=0;poly<cut_cells_np1(cell_index)->geometry.Size();poly++){
                POLYGON<TV>& polygon=cut_cells_np1(cell_index)->geometry(poly);
                T volume=polygon.Area();
                if(poly==cut_cells_np1(cell_index)->dominant_element){
                    cell_volumes_np1(cell_index)=volume;continue;}
                if(!cut_cells_np1(cell_index)->visibility(poly).Size()){
                    LOG::cout<<"No visible neighbors for cut cell ";for(int i=0;i<polygon.X.Size();i++) LOG::cout<<polygon.X(i)<<", ";LOG::cout<<"discarding "<<volume<<" volume"<<std::endl;
                    continue;}
                volume /= cut_cells_np1(cell_index)->visibility(poly).Size();
                for(int n=0;n<cut_cells_np1(cell_index)->visibility(poly).Size();n++)
                    cell_volumes_np1(cut_cells_np1(cell_index)->visibility(poly)(n)) += volume;}}

    uncovered_cells.Fill(false);
}
//#####################################################################
// Function Update_Np1_Collision_Data
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Update_Np1_Collision_Data(const T dt)
{
    collision_bodies_affecting_fluid->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
    collision_bodies_affecting_fluid->Rasterize_Objects(); // non-swept
    collision_bodies_affecting_fluid->Compute_Occupied_Blocks(false,(T)2*euler.grid.Minimum_Edge_Length(),5);  // static occupied blocks

    cut_cells_n.Delete_Pointers_And_Clean_Memory();cut_cells_n.Resize(euler.grid.Domain_Indices(1)); // TODO(jontg): Swap instead of Deleting, Resizing and Putting
    T_ARRAYS_CUT_CELLS::Put(cut_cells_np1,cut_cells_n);
    CUT_CELL_COMPUTATIONS::Compute_Cut_Geometries(euler.grid,1,*collision_bodies_affecting_fluid,cut_cells_np1);

    T_ARRAYS_BOOL::Put(psi_np1,psi_n); psi_np1.Fill(true);
    for(COLLISION_GEOMETRY_ID id(0);id<collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m;id++)
        if(collision_bodies_affecting_fluid->collision_geometry_collection.Is_Active(id)){
            COLLISION_GEOMETRY<TV>& collision_body=*(collision_bodies_affecting_fluid->collision_geometry_collection.bodies(id));
            T collision_thickness_over_two=(T).5*collision_bodies_affecting_fluid->collision_thickness;
            for(CELL_ITERATOR iterator(euler.grid,euler.grid.Clamp_To_Cell(collision_body.Axis_Aligned_Bounding_Box().Thickened(euler.grid.dX.Max()*(T)2),0));iterator.Valid();iterator.Next()){
                if(collision_body.Inside(iterator.Location(),collision_thickness_over_two)) psi_np1(iterator.Cell_Index())=false;}}

    T_ARRAYS_SCALAR::Put(cell_volumes_np1,cell_volumes_n); cell_volumes_np1.Fill(euler.grid.Cell_Size());
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(!psi_np1(cell_index) || (cut_cells_np1(cell_index) && !cut_cells_np1(cell_index)->dominant_element)) cell_volumes_np1(cell_index)=0;
        if(cut_cells_np1(cell_index))
            for(int poly=0;poly<cut_cells_np1(cell_index)->geometry.Size();poly++){
                POLYGON<TV>& polygon=cut_cells_np1(cell_index)->geometry(poly);
                T volume=polygon.Area();
                if(poly==cut_cells_np1(cell_index)->dominant_element){cell_volumes_np1(cell_index)=volume;continue;}
                if(!cut_cells_np1(cell_index)->visibility(poly).Size()){
                    LOG::cout<<"No visible neighbors for cut cell ";for(int i=0;i<polygon.X.Size();i++) LOG::cout<<polygon.X(i)<<", ";LOG::cout<<"discarding "<<volume<<" volume"<<std::endl;
                    continue;}
                volume /= cut_cells_np1(cell_index)->visibility(poly).Size();
                for(int n=0;n<cut_cells_np1(cell_index)->visibility(poly).Size();n++) cell_volumes_np1(cut_cells_np1(cell_index)->visibility(poly)(n)) += volume;}}

    uncovered_cells.Fill(false);
    collision_bodies_affecting_fluid->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
    for(COLLISION_GEOMETRY_ID id(0);id<collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m;id++)
        if(collision_bodies_affecting_fluid->collision_geometry_collection.Is_Active(id)){
            COLLISION_GEOMETRY<TV>& collision_body=*(collision_bodies_affecting_fluid->collision_geometry_collection.bodies(id));
            for(CELL_ITERATOR iterator(euler.grid,euler.grid.Clamp_To_Cell(collision_body.Axis_Aligned_Bounding_Box().Thickened(euler.grid.dX.Max()*(T)2),1));iterator.Valid();iterator.Next()){
                if(collision_body.Any_Simplex_Crossover(iterator.Location(),iterator.Location(),dt)) uncovered_cells(iterator.Cell_Index())=true;}}
}
//#####################################################################
// Function Compute_Intermediate_Solid_Position_Data
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Compute_Intermediate_Solid_Position_Data(const T dt)
{
    collision_bodies_affecting_fluid->Average_States(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_SAVED_NEW_STATE,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE, (T).5);
    collision_bodies_affecting_fluid->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
    collision_bodies_affecting_fluid->Rasterize_Objects(); // non-swept
    collision_bodies_affecting_fluid->Compute_Occupied_Blocks(false,(T)2*euler.grid.Minimum_Edge_Length(),5);  // static occupied blocks

    cut_cells_n_p_half.Delete_Pointers_And_Clean_Memory(); cut_cells_n_p_half.Resize(euler.grid.Domain_Indices(1));
    CUT_CELL_COMPUTATIONS::Compute_Cut_Geometries(euler.grid,1,*collision_bodies_affecting_fluid,cut_cells_n_p_half);

    psi_n_p_half.Fill(true);
    for(COLLISION_GEOMETRY_ID id(0);id<collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m;id++)
        if(collision_bodies_affecting_fluid->collision_geometry_collection.Is_Active(id)){
            COLLISION_GEOMETRY<TV>& collision_body=*(collision_bodies_affecting_fluid->collision_geometry_collection.bodies(id));
            T collision_thickness_over_two=(T).5*collision_bodies_affecting_fluid->collision_thickness;
            for(CELL_ITERATOR iterator(euler.grid,euler.grid.Clamp_To_Cell(collision_body.Axis_Aligned_Bounding_Box().Thickened(euler.grid.dX.Max()*(T)2),0));iterator.Valid();iterator.Next()){
                if(collision_body.Inside(iterator.Location(),collision_thickness_over_two)) psi_n_p_half(iterator.Cell_Index())=false;}}

    cell_volumes_n_p_half.Fill(euler.grid.Cell_Size());
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(!psi_n_p_half(cell_index) || (cut_cells_n_p_half(cell_index) && !cut_cells_n_p_half(cell_index)->dominant_element)) cell_volumes_n_p_half(cell_index)=0;
        if(cut_cells_n_p_half(cell_index))
            for(int poly=0;poly<cut_cells_n_p_half(cell_index)->geometry.Size();poly++){
                POLYGON<TV>& polygon=cut_cells_n_p_half(cell_index)->geometry(poly);
                T volume=polygon.Area();
                if(poly==cut_cells_n_p_half(cell_index)->dominant_element){cell_volumes_n_p_half(cell_index)=volume;continue;}
                if(!cut_cells_n_p_half(cell_index)->visibility(poly).Size()){
                    LOG::cout<<"No visible neighbors for cut cell ";for(int i=0;i<polygon.X.Size();i++) LOG::cout<<polygon.X(i)<<", ";LOG::cout<<"discarding "<<volume<<" volume"<<std::endl;
                    continue;}
                volume /= cut_cells_n_p_half(cell_index)->visibility(poly).Size();
                for(int n=0;n<cut_cells_n_p_half(cell_index)->visibility(poly).Size();n++) cell_volumes_n_p_half(cut_cells_n_p_half(cell_index)->visibility(poly)(n)) += volume;}}

    uncovered_cells_n_p_half.Fill(false);
    collision_bodies_affecting_fluid->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
    for(COLLISION_GEOMETRY_ID id(0);id<collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m;id++)
        if(collision_bodies_affecting_fluid->collision_geometry_collection.Is_Active(id)){
            COLLISION_GEOMETRY<TV>& collision_body=*(collision_bodies_affecting_fluid->collision_geometry_collection.bodies(id));
            for(CELL_ITERATOR iterator(euler.grid,euler.grid.Clamp_To_Cell(collision_body.Axis_Aligned_Bounding_Box().Thickened(euler.grid.dX.Max()*(T)2),1));iterator.Valid();iterator.Next()){
                if(collision_body.Any_Simplex_Crossover(iterator.Location(),iterator.Location(),(T).5*dt)) uncovered_cells_n_p_half(iterator.Cell_Index())=true;}}

    collision_bodies_affecting_fluid->Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_SAVED_NEW_STATE);
    collision_bodies_affecting_fluid->Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
}
//#####################################################################
// Function Revert_Cells_Near_Interface
//#####################################################################
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Revert_Cells_Near_Interface(const int iteration_number)
{
    if(iteration_number==1){
        for(CELL_ITERATOR iterator(euler.grid,1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            near_interface(cell_index) = (cut_cells_n(cell_index)!=0 || cut_cells_n_p_half(cell_index)!=0 || cut_cells_np1(cell_index)!=0);}
        for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            for(int dim=0;dim<TV::dimension;dim++) for(int n=-1;n<=1;n+=2){
                TV_INT neighbor_index=cell_index+n*TV_INT::Axis_Vector(dim);
                if((cut_cells_n(cell_index) && (!cut_cells_n(cell_index)->dominant_element ||
                                                !cut_cells_n(cell_index)->visibility(cut_cells_n(cell_index)->dominant_element).Contains(neighbor_index)))
                   || (cut_cells_n_p_half(cell_index) && (!cut_cells_n_p_half(cell_index)->dominant_element ||
                                                          !cut_cells_n_p_half(cell_index)->visibility(cut_cells_n_p_half(cell_index)->dominant_element).Contains(neighbor_index)))
                   || (cut_cells_np1(cell_index) && (!cut_cells_np1(cell_index)->dominant_element ||
                                                          !cut_cells_np1(cell_index)->visibility(cut_cells_np1(cell_index)->dominant_element).Contains(neighbor_index))))
                    near_interface(neighbor_index)=true;}}}
}
//#####################################################################
// Function Update_Cells_Near_Interface
//#####################################################################
namespace {
template<class T,int d>
void Add_Weight_To_Advection(const T weight, const VECTOR<int,d>& donor_cell, const VECTOR<int,d>& receiver_cell,ARRAY<PAIR<T,VECTOR<int,d> > >& weight_array,ARRAY<ARRAY<int>,VECTOR<int,d> >& donor_array,ARRAY<ARRAY<int>,VECTOR<int,d> >& receiver_array,ARRAY<T,VECTOR<int,d> >& sigma){
    int index = weight_array.Append(PAIR<T, VECTOR<int,d> >(weight,receiver_cell));
    sigma(donor_cell)+=weight;
    donor_array(donor_cell).Append(index);
    receiver_array(receiver_cell).Append(index);
    // LOG::cout<<"Adding weight "<<weight<<" from "<<donor_cell<<" to "<<receiver_cell<<"; sigma("<<donor_cell<<") = "<<sigma(donor_cell)<<std::endl;
}

template<class TV> void Fill_Fluid_Velocity(const GRID<TV>& grid, const ARRAY<bool,VECTOR<int,TV::dimension> >& psi, const ARRAY<bool,VECTOR<int,TV::dimension> >& psi_new, const ARRAY<bool,VECTOR<int,TV::dimension> >& swept_cells,
                                            const ARRAY<TV,VECTOR<int,TV::dimension> >& V_n, const ARRAY<VECTOR<typename TV::SCALAR,TV::dimension+2>,VECTOR<int,TV::dimension> >& U_new,ARRAY<TV,VECTOR<int,TV::dimension> >& velocity_field)
{
    typedef GRID<TV> T_GRID; typedef typename T_GRID::VECTOR_INT TV_INT; typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(psi(cell_index)) velocity_field(cell_index)=EULER<T_GRID>::Get_Velocity(U_new(cell_index));
        else if(swept_cells(cell_index)) velocity_field(cell_index)=V_n(cell_index);}
}

template<class TV> void Advect_Near_Interface_Data(const GRID<TV>& grid,const typename TV::SCALAR collision_thickness,const typename TV::SCALAR dt,const ARRAY<TV,VECTOR<int,TV::dimension> >& advection_velocity,const ARRAY<bool,VECTOR<int,TV::dimension> >& near_interface_mask,
                                                   const ARRAY<bool,VECTOR<int,TV::dimension> >& swept_cells,const ARRAY<bool,VECTOR<int,TV::dimension> >& psi_np1,ARRAY<VECTOR<typename TV::SCALAR,TV::dimension+2>,FACE_INDEX<TV::dimension> > flux_boundary_conditions,
                                                   const ARRAY<CUT_CELLS<typename TV::SCALAR,TV::dimension>*,VECTOR<int,TV::dimension> >& cut_cells_n,  const ARRAY<typename TV::SCALAR,VECTOR<int,TV::dimension> >& cell_volumes_n,  const ARRAY<VECTOR<typename TV::SCALAR,TV::dimension+2>,VECTOR<int,TV::dimension> >& U_n,
                                                   const ARRAY<CUT_CELLS<typename TV::SCALAR,TV::dimension>*,VECTOR<int,TV::dimension> >& cut_cells_np1,const ARRAY<typename TV::SCALAR,VECTOR<int,TV::dimension> >& cell_volumes_np1,      ARRAY<VECTOR<typename TV::SCALAR,TV::dimension+2>,VECTOR<int,TV::dimension> >& U_np1)
{   // TODO(jontg): Sparse data representation.
    typedef typename TV::SCALAR T; typedef VECTOR<T,TV::dimension+2> TV_DIMENSION; typedef GRID<TV> T_GRID;
    typedef typename T_GRID::VECTOR_INT TV_INT; typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;

    ARRAY<bool,VECTOR<int,TV::dimension> > near_interface(near_interface_mask);
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(!psi_np1(iterator.Cell_Index())) near_interface(iterator.Cell_Index())=false;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(near_interface(iterator.Cell_Index())) U_np1(iterator.Cell_Index())=TV_DIMENSION();

    ARRAY<TRIPLE<FACE_INDEX<TV::dimension>,TV_INT,TV_INT> > hybrid_boundary_flux;
    for(FACE_ITERATOR iterator(grid,0,T_GRID::INTERIOR_REGION);iterator.Valid();iterator.Next()){
        TV_INT first_cell_index=iterator.First_Cell_Index(), second_cell_index=iterator.Second_Cell_Index();
        if((near_interface(first_cell_index) ^ near_interface(second_cell_index)) && psi_np1(first_cell_index) && psi_np1(second_cell_index))
            hybrid_boundary_flux.Append(TRIPLE<FACE_INDEX<TV::dimension>,TV_INT,TV_INT>(iterator.Full_Index(),first_cell_index,second_cell_index));}

    for(int variable_index=0;variable_index<TV_DIMENSION::dimension;variable_index++){
        ARRAY<PAIR<T,TV_INT> > weights;
        ARRAY<ARRAY<int>,TV_INT> donors(grid.Domain_Indices(1));
        ARRAY<ARRAY<int>,TV_INT> receivers(grid.Domain_Indices(1));
        ARRAY<T,TV_INT> sigma(grid.Domain_Indices(1));

        {
            LOG::SCOPE scope("Converting Flux into weights.");
            for(int i=0;i<hybrid_boundary_flux.Size();i++){
                TV_INT first_cell_index=hybrid_boundary_flux(i).y,second_cell_index=hybrid_boundary_flux(i).z;
                T weight=dt*flux_boundary_conditions(hybrid_boundary_flux(i).x)(variable_index);
                Add_Weight_To_Advection(abs(weight), (weight >= 0 ? first_cell_index : second_cell_index), (weight >= 0 ? second_cell_index : first_cell_index), weights, donors, receivers, sigma);}
        }
 
        {
            LOG::SCOPE scope("Swept / uncovered cells.");
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();
            if(swept_cells(cell_index)){
                LOG::cout<<"Swept cell "<<cell_index<<" encountered"<<std::endl;
                if(near_interface(cell_index)){ // Backward cast to populate a newly uncovered cell
                    CUT_CELLS<T,TV::dimension>* cut_cells=cut_cells_np1(cell_index);
                    if(cut_cells) for(int poly=0;poly<cut_cells->geometry.Size();poly++){
                        if(cut_cells->visibility(poly).Contains(cell_index)){
                            if(cut_cells->visibility(poly).Size()<=1) LOG::cout<<"Lost material in uncovered cell "<<cell_index<<" which has no visible neighbors at t_np1"<<std::endl;
                            else{
                                T volume_fraction=cut_cells->geometry(poly).Area()/((T)(cut_cells->visibility(poly).Size()-1));
                                for(int visible_neighbor=0;visible_neighbor<cut_cells->visibility(poly).Size();visible_neighbor++){
                                    TV_INT neighbor=cut_cells->visibility(poly)(visible_neighbor);
                                    if(neighbor==cell_index) continue;
                                    else{
                                        T weight=U_n(neighbor)(variable_index)*volume_fraction;
                                        Add_Weight_To_Advection(weight, neighbor, cell_index, weights, donors, receivers, sigma);}}}}}
                    else LOG::cout<<"Cell "<<cell_index<<" is swept but is not a cut cell at t_np1!"<<std::endl;}

                CUT_CELLS<T,TV::dimension>* cut_cells=cut_cells_n(cell_index); // Push out any data from a newly covered cell
                if(cut_cells) for(int poly=0;poly<cut_cells->geometry.Size();poly++){
                        if(cut_cells->visibility(poly).Contains(cell_index)){
                            if(cut_cells->visibility(poly).Size()==1) LOG::cout<<"Lost material in uncovered cell "<<cell_index<<" which has no visible neighbors at t_n"<<std::endl;
                            else{
                                T volume_fraction=cut_cells->geometry(poly).Area()/((T)(cut_cells->visibility(poly).Size()-1));
                                for(int visible_neighbor=0;visible_neighbor<cut_cells->visibility(poly).Size();visible_neighbor++){
                                    TV_INT neighbor=cut_cells->visibility(poly)(visible_neighbor);
                                    if(neighbor==cell_index) continue;
                                    else{
                                        T weight=U_n(cell_index)(variable_index)*volume_fraction;
                                        Add_Weight_To_Advection(weight, cell_index, neighbor, weights, donors, receivers, sigma);}}}}}
                else LOG::cout<<"Cell "<<cell_index<<" is swept but has no place to put its t_n data!"<<std::endl;}}
        }

        // Backward-cast weights
        {
            LOG::SCOPE scope("Backward-cast weights");
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT receiver_cell=iterator.Cell_Index();
            if(near_interface(receiver_cell) && !swept_cells(receiver_cell)){
                TV cell_preimage_center=iterator.Location()-dt*advection_velocity(receiver_cell);
                // LOG::cout<<"Looking back from "<<receiver_cell<<" by "<<-dt*advection_velocity(receiver_cell)<<std::endl;
                RANGE<TV> cell_preimage(iterator.Bounding_Box()-iterator.Location()+cell_preimage_center);
                RANGE<TV_INT> affected_cells(grid.Index(cell_preimage.min_corner-collision_thickness),grid.Index(cell_preimage.max_corner+collision_thickness)+TV_INT::All_Ones_Vector());
                ARRAY<TV_INT>* visibility=0;
                if(cut_cells_np1(receiver_cell) && cut_cells_np1(receiver_cell)->dominant_element) visibility=&cut_cells_np1(receiver_cell)->visibility(cut_cells_np1(receiver_cell)->dominant_element);

                for(CELL_ITERATOR intersecting_iter(grid,affected_cells);intersecting_iter.Valid();intersecting_iter.Next()){
                    TV_INT donor_cell=intersecting_iter.Cell_Index();
                    CUT_CELLS<T,TV::dimension>* cut_cells=cut_cells_n(donor_cell);
                    if(!cut_cells){
                        if(swept_cells(donor_cell)) continue;// already handled in *forward* step above
                        if(near_interface(donor_cell) && (!visibility || visibility->Contains(donor_cell))){
                            T weight=U_n(donor_cell)(variable_index) * INTERSECTION::Intersection_Area(intersecting_iter.Bounding_Box(),cell_preimage);
                            Add_Weight_To_Advection(weight, donor_cell, receiver_cell, weights, donors, receivers, sigma);}}
                    else{
                        if(cut_cells->dominant_element && cut_cells->visibility(cut_cells->dominant_element).Contains(receiver_cell)){
                            if(swept_cells(donor_cell)) continue;// already handled in *forward* step above
                            T weight=U_n(donor_cell)(variable_index) * INTERSECTION::Intersection_Area(cell_preimage,cut_cells->geometry(cut_cells->dominant_element));
                            Add_Weight_To_Advection(weight, donor_cell, receiver_cell, weights, donors, receivers, sigma);}
                        else for(int poly=0;poly<cut_cells->geometry.Size();poly++){
                            if(cut_cells->visibility(poly).Contains(receiver_cell) && INTERSECTION::Intersection_Area(cell_preimage,cut_cells->geometry(poly)) >= (T)1e-6*cut_cells->geometry(poly).Area()){
                                T volume=cut_cells->geometry(poly).Area() / (T)cut_cells->visibility(poly).Size();
                                for(int visible_neighbor=0;visible_neighbor<cut_cells->visibility(poly).Size();visible_neighbor++){
                                    TV_INT neighbor=cut_cells->visibility(poly)(visible_neighbor);
                                    Add_Weight_To_Advection(U_n(neighbor)(variable_index) * volume, neighbor, receiver_cell, weights, donors, receivers, sigma);}}}}}}}
        }

        // Clamp and forward-cast weights
        {
            LOG::SCOPE scope("Forward-cast weights");
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT donor_cell=iterator.Cell_Index();
            if(near_interface(donor_cell) && !swept_cells(donor_cell)){
                T cell_stuff = U_n(donor_cell)(variable_index) * cell_volumes_n(donor_cell);
                T sigma_cell=sigma(donor_cell);
                if(abs(sigma_cell) > abs(cell_stuff)){
                    for(int i=0;i<donors(donor_cell).Size();i++) if(!near_interface(weights(donors(donor_cell)(i)).y)){cell_stuff -= weights(donors(donor_cell)(i)).x; sigma_cell -= weights(donors(donor_cell)(i)).x;}
                    if(sigma_cell <= (T)1e-10) LOG::cout<<"Not clamping for tiny sigma "<<sigma_cell<<" from "<<sigma(donor_cell)<<" total cell stuff (minus boundary conditions) = "<<cell_stuff<<std::endl;
                    else{
                        T one_over_sigma = cell_stuff / sigma_cell;
                        // LOG::cout<<"Clamping weights leaving "<<donor_cell<<" by factor "<<one_over_sigma<<" from "<<sigma(donor_cell)<<" total cell stuff (minus boundary conditions) = "<<cell_stuff<<std::endl;
                        for(int i=0;i<donors(donor_cell).Size();i++) if(near_interface(weights(donors(donor_cell)(i)).y)) weights(donors(donor_cell)(i)).x *= one_over_sigma;}}
                else{
                    T remainder = (cell_stuff - sigma_cell);
                    TV cell_postimage_center=iterator.Location()+dt*advection_velocity(donor_cell);
                    // LOG::cout<<"Looking forward from "<<donor_cell<<" by "<<dt*advection_velocity(donor_cell)<<std::endl;
                    RANGE<TV> cell_postimage(iterator.Bounding_Box()-iterator.Location()+cell_postimage_center);
                    RANGE<TV_INT> affected_cells(grid.Index(cell_postimage.min_corner-collision_thickness),grid.Index(cell_postimage.max_corner+collision_thickness)+TV_INT::All_Ones_Vector());
                    ARRAY<TV_INT>* visibility=0;
                    if(cut_cells_n(donor_cell) && cut_cells_n(donor_cell)->dominant_element)
                         visibility=&cut_cells_n(donor_cell)->visibility(cut_cells_n(donor_cell)->dominant_element);

                    ARRAY<PAIR<TV_INT,T> > forward_weights;
                    T distributed_volume=0;
                    for(CELL_ITERATOR intersecting_iter(grid,affected_cells);intersecting_iter.Valid();intersecting_iter.Next()){
                        TV_INT receiver_cell=intersecting_iter.Cell_Index();
                        if(swept_cells(receiver_cell)) continue; // already handled in *backward* step above
                        CUT_CELLS<T,TV::dimension>* cut_cells=cut_cells_np1(receiver_cell);
                        if(!cut_cells){
                            if(near_interface(receiver_cell) && (!visibility || visibility->Contains(donor_cell))){
                                T weight=INTERSECTION::Intersection_Area(intersecting_iter.Bounding_Box(),cell_postimage);
                                distributed_volume+=weight; forward_weights.Append(PAIR<TV_INT,T>(receiver_cell,weight));}}
                        else{
                            if(cut_cells->dominant_element && cut_cells->visibility(cut_cells->dominant_element).Contains(donor_cell)){
                               T weight=INTERSECTION::Intersection_Area(cell_postimage,cut_cells->geometry(cut_cells->dominant_element));
                               distributed_volume+=weight; forward_weights.Append(PAIR<TV_INT,T>(receiver_cell,weight));}
                            else for(int poly=0;poly<cut_cells->geometry.Size();poly++){
                                if(cut_cells->visibility(poly).Contains(donor_cell)){
                                    T volume=INTERSECTION::Intersection_Area(cell_postimage,cut_cells->geometry(poly)) / (T)cut_cells->visibility(poly).Size();
                                    for(int visible_neighbor=0;visible_neighbor<cut_cells->visibility(poly).Size();visible_neighbor++){
                                        TV_INT neighbor=cut_cells->visibility(poly)(visible_neighbor);
                                        distributed_volume+=volume; forward_weights.Append(PAIR<TV_INT,T>(neighbor,volume));}}}}}

                    if(distributed_volume <= (T)1e-4*cell_volumes_n(donor_cell)){LOG::cout<<"Forward-cast failed to distribute a significant fraction of its volume..."<<std::endl;assert(false);} // TODO(jontg): Handle this case somehow...
                    T dilation_factor=remainder / distributed_volume;

                    for(int i=0;i<forward_weights.Size();i++) Add_Weight_To_Advection(dilation_factor*forward_weights(i).y, donor_cell, forward_weights(i).x, weights, donors, receivers, sigma);}}}
        }

        for(int i=0;i<weights.Size();i++){
            PAIR<T,TV_INT>& weight_pair=weights(i);
            if(near_interface(weight_pair.y)) U_np1(weight_pair.y)(variable_index) += weight_pair.x;}}

    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(near_interface(cell_index)){U_np1(cell_index) /= cell_volumes_np1(cell_index);
            LOG::cout<<"Cell "<<cell_index<<" of volume "<<cell_volumes_np1(cell_index)<<" set to "<<U_np1(cell_index)<<std::endl;}}
}
};
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Update_Cells_Near_Interface(const T dt,const int rk_order,const int rk_substep)
{
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(near_interface(cell_index) && !psi_np1(cell_index)) near_interface(cell_index)=false;}

    const int saved_order=euler.conservation->order;

    T_FACE_ARRAYS_INT highest_eno_order(euler.grid,saved_order);highest_eno_order.Fill(saved_order);
    for(FACE_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
        if(near_interface(iterator.First_Cell_Index()) ^ near_interface(iterator.Second_Cell_Index())){
            highest_eno_order(iterator.Full_Index())=1;
            const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
            for(int offset=1;offset<saved_order;++offset){
                highest_eno_order.Component(axis)(face_index-offset*TV_INT::Axis_Vector(axis))=min(offset+1,highest_eno_order.Component(axis)(face_index-offset*TV_INT::Axis_Vector(axis)));
                highest_eno_order.Component(axis)(face_index+offset*TV_INT::Axis_Vector(axis))=min(offset+1,highest_eno_order.Component(axis)(face_index+offset*TV_INT::Axis_Vector(axis)));}}}

    for(FACE_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face_index=iterator.Full_Index();
        if(highest_eno_order(face_index)==saved_order) accumulated_flux(face_index) += euler.conservation->fluxes(face_index);}

    for(int ord=saved_order-2;ord>=0;--ord){
        euler.conservation->Set_Order(ord);
        euler.conservation->Update_Conservation_Law(euler.grid,euler.U,euler.U_ghost,euler.psi,dt,euler.eigensystems,euler.eigensystems_default,euler.euler_projection.elliptic_solver->psi_N,
                                                    euler.euler_projection.face_velocities,false,euler.open_boundaries);
        for(FACE_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face_index=iterator.Full_Index();
            if(highest_eno_order(face_index)==ord && !(near_interface(iterator.First_Cell_Index()) && near_interface(iterator.Second_Cell_Index())))
                accumulated_flux(face_index) += euler.conservation->fluxes(face_index);}}
    euler.conservation->Set_Order(saved_order);

    if(rk_order==2 && rk_substep==2) accumulated_flux *= (T).5;
    else if(rk_order==3){
        if(rk_substep==2) accumulated_flux *= (T).5;
        if(rk_substep==3) accumulated_flux *= (T)2./3;}

    const TV dt_over_dx=(rk_order==3 && rk_substep==2 ? dt/(T)2 : dt)*euler.grid.One_Over_DX();
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(euler.psi(cell_index)) euler.U(cell_index)=U_n(cell_index);}
    for(FACE_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
        TV_INT first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
        if(euler.U.Valid_Index(first_cell)  && !near_interface(first_cell))  euler.U(first_cell)  -= accumulated_flux(iterator.Full_Index())*dt_over_dx(iterator.Axis());
        if(euler.U.Valid_Index(second_cell) && !near_interface(second_cell)) euler.U(second_cell) += accumulated_flux(iterator.Full_Index())*dt_over_dx(iterator.Axis());}

    if(rk_order==3 && rk_substep==2)
        Advect_Near_Interface_Data(euler.grid,collision_bodies_affecting_fluid->collision_thickness,dt/(T)2,advection_velocities_n,near_interface,uncovered_cells_n_p_half,psi_n_p_half,accumulated_flux,
                                   cut_cells_n,cell_volumes_n,U_n,cut_cells_n_p_half,cell_volumes_n_p_half, euler.U);
    else
        Advect_Near_Interface_Data(euler.grid,collision_bodies_affecting_fluid->collision_thickness,dt,advection_velocities_n,near_interface,uncovered_cells,psi_np1,accumulated_flux,
                                   cut_cells_n,cell_volumes_n,U_n,cut_cells_np1,cell_volumes_np1, euler.U);

    euler.Invalidate_Ghost_Cells();
    if(rk_order==3 && rk_substep==2) accumulated_flux *= (T).5;
}
template<class TV> void SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>::
Compute_Post_Advected_Variables()
{
    TV_DIMENSION accumulated_material=TV_DIMENSION();
    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next())
        if(psi_np1(iterator.Cell_Index())) accumulated_material+=euler.U(iterator.Cell_Index())*cell_volumes_np1(iterator.Cell_Index());
    LOG::cout<<std::setprecision(16)<<"TOTAL ACCUMULATED MATERIAL POST FSI ADVECTION = "<<accumulated_material<<std::endl;

    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(near_interface(cell_index)){
            if(uncovered_cells(cell_index)){T one_over_c=euler.eos->one_over_c(euler.U(cell_index)(1),euler.e(euler.U(cell_index)));
                LOG::cout<<"Computing 1/(rho c^2) for uncovered cell "<<cell_index<<" = "<<one_over_c*one_over_c/(euler.U(cell_index)(1))<<std::endl;
                euler.euler_projection.one_over_rho_c_squared(cell_index) = one_over_c*one_over_c/(euler.U(cell_index)(1));}
            // LOG::cout<<"Updated Cell "<<cell_index<<" from "<<U_n(cell_index)<<" to "<<euler.U(cell_index)<<"\t\t 1/(rho c^2) = "<<euler.euler_projection.one_over_rho_c_squared(cell_index)<<std::endl;
        }}
}
//#####################################################################
template class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<VECTOR<float,1> >;
template class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<VECTOR<float,2> >;
template class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<VECTOR<double,1> >;
template class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<VECTOR<double,2> >;
template class SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<VECTOR<double,3> >;
#endif
