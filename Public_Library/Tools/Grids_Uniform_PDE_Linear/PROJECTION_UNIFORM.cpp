//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PROJECTION_UNIFORM<TV>::
PROJECTION_UNIFORM(const GRID<TV>& mac_grid,const bool use_variable_beta,const bool use_poisson,THREAD_QUEUE* thread_queue_input)
    :use_divergence_multiplier(false),thread_queue(thread_queue_input)
{
    if(use_variable_beta || use_poisson){
        poisson=new POISSON_UNIFORM<TV>(p_grid,p,true,false,true);
        if(use_variable_beta) poisson->Set_Variable_beta();elliptic_solver=poisson;laplace=0;}
    else{laplace=new LAPLACE_UNIFORM<TV>(p_grid,p,true,true,thread_queue);elliptic_solver=laplace;poisson=0;}
    Initialize_Grid(mac_grid);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PROJECTION_UNIFORM<TV>::
PROJECTION_UNIFORM(THREAD_QUEUE* thread_queue_input)
    :elliptic_solver(0),laplace(0),poisson(0),use_divergence_multiplier(false),thread_queue(thread_queue_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PROJECTION_UNIFORM<TV>::
~PROJECTION_UNIFORM()
{
    delete elliptic_solver;
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class TV> void PROJECTION_UNIFORM<TV>::
Initialize_Grid(const GRID<TV>& mac_grid)
{
    assert(mac_grid.Is_MAC_Grid());p_grid=mac_grid;
    if(poisson) poisson->Initialize_Grid(p_grid);else laplace->Initialize_Grid(p_grid);
    Use_Non_Zero_Divergence(use_non_zero_divergence); // call this since the grid changed
    p.Resize(p_grid.Domain_Indices(1));p_save_for_projection.Resize(p_grid.Domain_Indices(1));face_velocities_save_for_projection.Resize(p_grid);
}
//#####################################################################
// Function Make_Divergence_Free
//#####################################################################
template<class TV> void PROJECTION_UNIFORM<TV>::
Make_Divergence_Free(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    // find f - divergence of the velocity
    Compute_Divergence(T_FACE_LOOKUP(face_velocities),elliptic_solver);

    // find the pressure
    elliptic_solver->Find_Solution_Regions(); // flood fill
    elliptic_solver->Compute_beta_And_Add_Jumps_To_b(dt,time); // only does something for poisson solver
    if(elliptic_solver->solve_neumann_regions) Enforce_Velocity_Compatibility(face_velocities); // make all the right hand sides compatible
    elliptic_solver->Solve(time,true); // solve all regions

    Apply_Pressure(face_velocities,dt,time);
}
//#####################################################################
// Function Zero_Out_Neumann_Pocket_Velocities
//#####################################################################
template<class TV> void PROJECTION_UNIFORM<TV>::
Zero_Out_Neumann_Pocket_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities)
{
    // zero out the velocities in neumann pockets to prevent gravity from accumulating
    T_FACE_ARRAYS_BOOL &psi_N=elliptic_solver->psi_N;
    if(!elliptic_solver->solve_neumann_regions){ 
        for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
            int color1=elliptic_solver->filled_region_colors(first_cell),color2=elliptic_solver->filled_region_colors(second_cell);
            if(color1==-1 || color2==-1 || psi_N.Component(axis)(face_index)) continue;
            if(!elliptic_solver->filled_region_touches_dirichlet(color1) && !elliptic_solver->filled_region_touches_dirichlet(color2))
                face_velocities.Component(axis)(face_index)=0;}}
}
//#####################################################################
// Function Apply_Pressure
//#####################################################################
template<class TV> void PROJECTION_UNIFORM<TV>::
Apply_Pressure(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,bool scale_by_dt)
{
    Zero_Out_Neumann_Pocket_Velocities(face_velocities);
    // find divergence free u, v and w 
    TV dx=p_grid.dX,one_over_dx=Inverse(dx);
    ARRAY<bool,TV_INT>& psi_D=elliptic_solver->psi_D;
    T_FACE_ARRAYS_BOOL& psi_N=elliptic_solver->psi_N;
    if(scale_by_dt) p*=dt;
    if(laplace) 
        for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
            if(!psi_N.Component(axis)(face_index) && !(psi_D(first_cell) && psi_D(second_cell)))
                face_velocities.Component(axis)(face_index)-=(p(second_cell)-p(first_cell))*one_over_dx[axis];}
    else if(poisson)
        for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
            if(!psi_N.Component(axis)(face_index) && !(psi_D(first_cell) && psi_D(second_cell)))
                face_velocities.Component(axis)(face_index)-=poisson->beta_face.Component(axis)(face_index)*(p(second_cell)-p(first_cell))*one_over_dx[axis];}
    if(scale_by_dt) p*=1/dt;
}
//#####################################################################
// Function Compute_Divergence
//#####################################################################
template<class TV> void PROJECTION_UNIFORM<TV>::
Compute_Divergence(const T_FACE_LOOKUP& face_lookup,LAPLACE_UNIFORM<TV>* solver)
{
    DOMAIN_ITERATOR_THREADED_ALPHA<PROJECTION_UNIFORM<TV>,TV>(p_grid.Domain_Indices(),thread_queue).template Run<const T_FACE_LOOKUP&,LAPLACE_UNIFORM<TV>*>(*this,&PROJECTION_UNIFORM<TV>::Compute_Divergence_Threaded,face_lookup,solver);
}
template<class TV> void PROJECTION_UNIFORM<TV>::
Compute_Divergence_Threaded(RANGE<TV_INT>& domain,const T_FACE_LOOKUP& face_lookup,LAPLACE_UNIFORM<TV>* solver)
{
    TV one_over_dx=p_grid.one_over_dX;
    for(CELL_ITERATOR<TV> iterator(p_grid,domain);iterator.Valid();iterator.Next()){
        const typename T_FACE_LOOKUP::LOOKUP& lookup=face_lookup.Starting_Point_Cell(iterator.Cell_Index());T divergence=0;
        for(int axis=0;axis<TV::m;axis++)divergence+=(lookup(axis,iterator.Second_Face_Index(axis))-lookup(axis,iterator.First_Face_Index(axis)))*one_over_dx[axis];
        solver->f(iterator.Cell_Index())=divergence;}

    if(use_non_zero_divergence) for(CELL_ITERATOR<TV> iterator(p_grid,domain);iterator.Valid();iterator.Next())
        solver->f(iterator.Cell_Index())-=divergence(iterator.Cell_Index());
    if(use_divergence_multiplier) for(CELL_ITERATOR<TV> iterator(p_grid,domain);iterator.Valid();iterator.Next())
        solver->f(iterator.Cell_Index())*=divergence_multiplier(iterator.Cell_Index());
}
//#####################################################################
// Function Enforce_Velocity_Compatibility
//#####################################################################
// For each Neumann region, modifies divergence in cells touching its boundary so that the sum of f in the region is zero. Also modifies each bounding face velocity to the desired velocity
// (or the average of the desired velocity if the face joins two different Neumann regions).
template<class TV> void PROJECTION_UNIFORM<TV>::
Enforce_Velocity_Compatibility(T_FACE_ARRAYS_SCALAR& face_velocities)
{
    bool found_neumann_region=false;
    for(int color=0;color<elliptic_solver->number_of_regions;color++)if(!elliptic_solver->filled_region_touches_dirichlet(color)){found_neumann_region=true;break;}
    if(!found_neumann_region) return; // don't need to do any of the following

    TV one_over_dx=Inverse(p_grid.dX);T cell_size=p_grid.Cell_Size();

    ARRAY<double> compatibility_fraction(elliptic_solver->number_of_regions),boundary_size(elliptic_solver->number_of_regions);
    // calculate the compatibility errors
    for(CELL_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
        int color=elliptic_solver->filled_region_colors(iterator.Cell_Index());
        if(color>0 && !elliptic_solver->filled_region_touches_dirichlet(color)) compatibility_fraction(color)+=elliptic_solver->f(iterator.Cell_Index());}

    for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
        int color1=elliptic_solver->filled_region_colors(iterator.First_Cell_Index()),color2=elliptic_solver->filled_region_colors(iterator.Second_Cell_Index());
        if(color1==color2) continue;T face_size=iterator.Face_Size();
        if(color1>0 && !elliptic_solver->filled_region_touches_dirichlet(color1)) boundary_size(color1)+=face_size;
        if(color2>0 && !elliptic_solver->filled_region_touches_dirichlet(color2)) boundary_size(color2)+=face_size;}

    // sum up the boundary sizes and errors with the other nodes
    if(elliptic_solver->mpi_grid){
        bool found_global_neumann_region=false;
        for(int color=0;color<elliptic_solver->laplace_mpi->number_of_global_regions;color++)
            if(!elliptic_solver->filled_region_touches_dirichlet(color)){found_global_neumann_region=true;break;}
        if(found_global_neumann_region){
            ARRAY<T> global_compatibility_fraction(elliptic_solver->laplace_mpi->number_of_global_regions),global_boundary_size(elliptic_solver->laplace_mpi->number_of_global_regions);
            for(int i=0;i<elliptic_solver->laplace_mpi->number_of_global_regions;i++){global_compatibility_fraction(i)=(T)compatibility_fraction(i);global_boundary_size(i)=(T)boundary_size(i);}
            ARRAY<T> global_compatibility_fraction_result(elliptic_solver->laplace_mpi->number_of_global_regions),global_boundary_size_result(elliptic_solver->laplace_mpi->number_of_global_regions);
            elliptic_solver->mpi_grid->Reduce_Add(global_compatibility_fraction,global_compatibility_fraction_result);
            elliptic_solver->mpi_grid->Reduce_Add(global_boundary_size,global_boundary_size_result);
            for(int i=0;i<elliptic_solver->laplace_mpi->number_of_global_regions;i++){compatibility_fraction(i)=global_compatibility_fraction_result(i);boundary_size(i)=global_boundary_size_result(i);}}}

    for(int i=0;i<elliptic_solver->number_of_regions;i++) if(boundary_size(i)) compatibility_fraction(i)*=cell_size/boundary_size(i);

    // adjust the compatibility error to zero
    for(FACE_ITERATOR<TV> iterator(p_grid);iterator.Valid();iterator.Next()){
        int color1=elliptic_solver->filled_region_colors(iterator.First_Cell_Index()),color2=elliptic_solver->filled_region_colors(iterator.Second_Cell_Index());
        if(color1<=0 || elliptic_solver->filled_region_touches_dirichlet(color1)) color1=0;
        if(color2<=0 || elliptic_solver->filled_region_touches_dirichlet(color2)) color2=0;
        if(color1!=color2){T u_delta1=0,u_delta2=0;
            if(color1){u_delta1=-(T)compatibility_fraction(color1);elliptic_solver->f(iterator.First_Cell_Index())+=u_delta1*one_over_dx[iterator.Axis()];}
            if(color2){u_delta2=(T)compatibility_fraction(color2);elliptic_solver->f(iterator.Second_Cell_Index())-=u_delta2*one_over_dx[iterator.Axis()];}
            T u_delta=u_delta1+u_delta2;if(color1&&color2)u_delta*=(T).5;face_velocities.Component(iterator.Axis())(iterator.Face_Index())+=u_delta;}}
}
//#####################################################################
// Function Set_Up_For_Projection
//#####################################################################
template<class TV> void PROJECTION_UNIFORM<TV>::
Set_Up_For_Projection(T_FACE_ARRAYS_SCALAR& face_velocities)
{
    face_velocities_save_for_projection.Copy(face_velocities);
}
//#####################################################################
// Function Restore_After_Projection
//#####################################################################
template<class TV> void PROJECTION_UNIFORM<TV>::
Restore_After_Projection(T_FACE_ARRAYS_SCALAR& face_velocities)
{
    T_FACE_ARRAYS_SCALAR::Exchange(face_velocities,face_velocities_save_for_projection);
}
//#####################################################################
// Function Exchange_Pressures_For_Projection
//#####################################################################
template<class TV> void PROJECTION_UNIFORM<TV>::
Exchange_Pressures_For_Projection()
{
    ARRAY<T,TV_INT>::Exchange(p,p_save_for_projection);
}
//#####################################################################
// Function Calculate_Kinetic_Energy_Error
//#####################################################################
template<class TV> void PROJECTION_UNIFORM<TV>::
Calculate_Kinetic_Energy_Error(T_FACE_ARRAYS_SCALAR& face_velocities,ARRAY<TV,TV_INT>& kinetic_energy_error)
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
namespace PhysBAM{
template class PROJECTION_UNIFORM<VECTOR<float,1> >;
template class PROJECTION_UNIFORM<VECTOR<float,2> >;
template class PROJECTION_UNIFORM<VECTOR<float,3> >;
template class PROJECTION_UNIFORM<VECTOR<double,1> >;
template class PROJECTION_UNIFORM<VECTOR<double,2> >;
template class PROJECTION_UNIFORM<VECTOR<double,3> >;
}
