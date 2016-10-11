//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Frank Losasso, Nick Rasmussen, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Grid_PDE/Interpolation/AVERAGING_UNIFORM.h>
#include <Grid_PDE/Poisson/LAPLACE_UNIFORM.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/LAPLACE_COLLIDABLE_UNIFORM.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/IMPLICIT_VISCOSITY_UNIFORM.h>
#include <Dynamics/Heat_Flows/HEAT_LAPLACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> IMPLICIT_VISCOSITY_UNIFORM<TV>::
IMPLICIT_VISCOSITY_UNIFORM(LAPLACE_UNIFORM<TV>& elliptic_solver_input,const ARRAY<T,TV_INT>& variable_viscosity_input,const T density_input,const T viscosity_input,T_MPI_GRID* mpi_grid_input,
    const int axis_input,bool use_variable_viscosity_input,bool use_psi_R_input)
    :elliptic_solver(elliptic_solver_input),variable_viscosity(variable_viscosity_input),density(density_input),viscosity(viscosity_input),
    mpi_grid(mpi_grid_input),axis(axis_input),heat_solver(0),use_variable_viscosity(use_variable_viscosity_input),use_psi_R(use_psi_R_input)
{
    if(mpi_grid) face_grid=mpi_grid->Get_Non_Overlapping_Face_Grid(axis);
    else face_grid=elliptic_solver.grid.Get_Face_MAC_Grid(axis);
    u.Resize(face_grid.Domain_Indices(mpi_grid?1:0));
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_VISCOSITY_UNIFORM<TV>::
~IMPLICIT_VISCOSITY_UNIFORM()
{
    delete heat_solver;
}
//#####################################################################
// Function Viscous_Update
//#####################################################################
template<class TV> void IMPLICIT_VISCOSITY_UNIFORM<TV>::
Viscous_Update(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,const T dt,const T time,const int maximum_implicit_viscosity_iterations)
{
    const T_ARRAYS_BASE& velocity_component_ghost=face_velocities_ghost.Component(axis);

    if(!heat_solver){Allocate_Heat_Solver();heat_solver->mpi_grid=mpi_grid;}
    heat_solver->pcg.Set_Maximum_Iterations(maximum_implicit_viscosity_iterations);

    Setup_Viscosity(dt);
    Setup_Boundary_Conditions(face_velocities_ghost);

    // copy velocities into u so unchanged velocities can be copied out unconditionally...
    // TODO: make shared face copied in a more clever way to avoid this copy 
    for(FACE_ITERATOR<TV> iterator(grid,0,GRID<TV>::WHOLE_REGION,-1,axis);iterator.Valid();iterator.Next()){
        TV_INT p_face_index=iterator.Face_Index(),cell_index=p_face_index;
        u(cell_index)=velocity_component_ghost(p_face_index);}

    for(CELL_ITERATOR<TV> iterator(face_grid);iterator.Valid();iterator.Next()) heat_solver->f(iterator.Cell_Index())=velocity_component_ghost(iterator.Cell_Index());
    heat_solver->f*=-1;
    heat_solver->Compute_beta_And_Add_Jumps_To_b(dt,time);
    heat_solver->Solve(time);

    // update face_velocities
    T_ARRAYS_BASE& velocity_component=face_velocities.Component(axis);
    for(FACE_ITERATOR<TV> iterator(grid,0,GRID<TV>::WHOLE_REGION,-1,axis);iterator.Valid();iterator.Next()){
        TV_INT index=iterator.Face_Index();
        velocity_component(index)=u(index);}
}
//#####################################################################
// Function Variable_Viscosity_Explicit_Part
//#####################################################################
template<class TV> void IMPLICIT_VISCOSITY_UNIFORM<TV>::
Variable_Viscosity_Explicit_Part(const T density,const ARRAY<T,TV_INT>& variable_viscosity,const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,const T dt,const T time)
{
    TV one_over_DX=grid.one_over_dX;
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index(),first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
        for(int term_axis=0;term_axis<TV::m;term_axis++){
            TV_INT positive_base_index=face_index+TV_INT::Axis_Vector(term_axis);
            T velocity_deriv_plus=one_over_DX[axis]*(face_velocities_ghost(term_axis,positive_base_index)-face_velocities_ghost(term_axis,positive_base_index-TV_INT::Axis_Vector(axis))),
                velocity_deriv_minus=one_over_DX[axis]*(face_velocities_ghost(term_axis,face_index)-face_velocities_ghost(term_axis,face_index-TV_INT::Axis_Vector(axis)));
            T viscosity_plus,viscosity_minus;
            if(term_axis==axis){viscosity_plus=variable_viscosity(second_cell);viscosity_minus=variable_viscosity(first_cell);}
            else{ // stencil points for derivative lie on node in 2D or edge in 3D.
                viscosity_plus=(T).25*(variable_viscosity(first_cell)+variable_viscosity(second_cell)+
                    variable_viscosity(first_cell+TV_INT::Axis_Vector(term_axis))+variable_viscosity(second_cell+TV_INT::Axis_Vector(term_axis)));
                viscosity_minus=(T).25*(variable_viscosity(first_cell)+variable_viscosity(second_cell)+
                    variable_viscosity(first_cell-TV_INT::Axis_Vector(term_axis))+variable_viscosity(second_cell-TV_INT::Axis_Vector(term_axis)));}
            T term=dt/density*one_over_DX[term_axis]*(viscosity_plus*velocity_deriv_plus-viscosity_minus*velocity_deriv_minus);
            face_velocities(axis,face_index)+=term;}}
}
//#####################################################################
// Function Allocate_Heat_Solver
//#####################################################################
template<class TV> void IMPLICIT_VISCOSITY_UNIFORM<TV>::
Allocate_Heat_Solver()
{
    if(use_variable_viscosity) heat_solver=new HEAT_LAPLACE<POISSON_COLLIDABLE_UNIFORM<TV> >(face_grid,u);
    else heat_solver=new HEAT_LAPLACE<LAPLACE_COLLIDABLE_UNIFORM<TV> >(face_grid,u);
    heat_solver->Use_Psi_R();
    heat_solver->mpi_grid=mpi_grid;
    heat_solver->pcg.show_results=elliptic_solver.pcg.show_results;
}
//#####################################################################
// Function Setup_Viscosity
//#####################################################################
template<class TV> void IMPLICIT_VISCOSITY_UNIFORM<TV>::
Setup_Viscosity(const T dt)
{
    if(HEAT_LAPLACE<POISSON_COLLIDABLE_UNIFORM<TV> >* heat_poisson=dynamic_cast<HEAT_LAPLACE<POISSON_COLLIDABLE_UNIFORM<TV> >*>(heat_solver)){
        heat_poisson->Set_Variable_beta();
        T half_dt_over_density=(T).5*dt/density;
        for(CELL_ITERATOR<TV> iterator(face_grid);iterator.Valid();iterator.Next())
            heat_poisson->variable_beta(iterator.Cell_Index())=half_dt_over_density*(variable_viscosity(iterator.Cell_Index())+variable_viscosity(iterator.Cell_Index()-TV_INT::Axis_Vector(axis)));}
    else if(HEAT_LAPLACE<LAPLACE_COLLIDABLE_UNIFORM<TV> >* heat_laplace=dynamic_cast<HEAT_LAPLACE<LAPLACE_COLLIDABLE_UNIFORM<TV> >*>(heat_solver)) {
        T dt_times_kinematic_viscosity=dt*viscosity/density;
        heat_laplace->coefficient=dt_times_kinematic_viscosity;}
}
//#####################################################################
// Function Setup_Boundary_Conditions
//#####################################################################
template<class TV> void IMPLICIT_VISCOSITY_UNIFORM<TV>::
Setup_Boundary_Conditions(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities)
{
    const ARRAY<bool,FACE_INDEX<TV::m> >& p_psi_N=elliptic_solver.psi_N;const ARRAY<bool,TV_INT>& p_psi_D=elliptic_solver.psi_D;
    ARRAY<bool,FACE_INDEX<TV::m> >& psi_N=heat_solver->psi_N;ARRAY<bool,TV_INT>& psi_D=heat_solver->psi_D;
    const ARRAY<T,FACE_INDEX<TV::m> >& p_psi_R=elliptic_solver.psi_R;
    ARRAY<T,FACE_INDEX<TV::m> >& psi_R=heat_solver->psi_R;

    // dirichlet velocity boundary condition where we have neumann condition for the pressure
    for(CELL_ITERATOR<TV> iterator(face_grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index(),p_face_index=cell_index;
        u(cell_index)=face_velocities(axis,p_face_index); // always copy over current velocity to get a good initial guess
        if(p_psi_N(axis,p_face_index)) psi_D(cell_index)=true;} // TODO -- make sure face velocity bc are correct after advection
    if(heat_solver->mpi_grid) heat_solver->mpi_grid->Exchange_Boundary_Cell_Data(*heat_solver->mpi_grid,face_grid,psi_D,1,false);

    // set neumann b.c. around a cell if its corresponding face in the pressure grid has psi_D set on both sides
    for(CELL_ITERATOR<TV> iterator(face_grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index(),p_face_index=cell_index;
        TV_INT p_cell_index_1,p_cell_index_2;GRID<TV>::Cells_Touching_Face(axis,p_face_index,p_cell_index_1,p_cell_index_2);
        if(p_psi_D(p_cell_index_1)&&p_psi_D(p_cell_index_2)) for(int face_axis=0;face_axis<TV::m;face_axis++){
            psi_N(face_axis,face_grid.First_Face_Index_In_Cell(face_axis,cell_index))=psi_N(face_axis,face_grid.Second_Face_Index_In_Cell(face_axis,cell_index))=true;}}

    // do the same for the sides tangential to the axis being iterated over
    for(int grid_axis=0;grid_axis<TV::m;grid_axis++)if(grid_axis!=axis)for(int side=0;side<2;side++){
        for(CELL_ITERATOR<TV> iterator(face_grid,1,GRID<TV>::GHOST_REGION,2*grid_axis+side);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index(),p_face_index=cell_index;
            TV_INT p_cell_index_1,p_cell_index_2;GRID<TV>::Cells_Touching_Face(axis,p_face_index,p_cell_index_1,p_cell_index_2);
            if(p_psi_D.Valid_Index(p_cell_index_1)&&p_psi_D.Valid_Index(p_cell_index_2)&&p_psi_D(p_cell_index_1)&&p_psi_D(p_cell_index_2))
                for(int face_axis=0;face_axis<TV::m;face_axis++){
                    psi_N(face_axis,face_grid.First_Face_Index_In_Cell(face_axis,cell_index))=psi_N(face_axis,face_grid.Second_Face_Index_In_Cell(face_axis,cell_index))=true;}}}

    // set neumann for the faces in the same axis that are on dirichlet primal cells
    const TV_INT axis_offset=TV_INT::Axis_Vector(axis);
    for(FACE_ITERATOR<TV> iterator(face_grid,0,GRID<TV>::WHOLE_REGION,-1,axis);iterator.Valid();iterator.Next())if(p_psi_D(iterator.Face_Index()-axis_offset)) 
        psi_N(axis,iterator.Face_Index())=true;

    // set slip boundary conditions for tangential walls
    for(int other_axis=0;other_axis<TV::m;other_axis++) if(other_axis!=axis){
        for(FACE_ITERATOR<TV> iterator(face_grid,0,GRID<TV>::BOUNDARY_REGION,-1,other_axis);iterator.Valid();iterator.Next()){TV_INT face_index=iterator.Face_Index(),p_node_index=face_index;
            if(p_psi_N(other_axis,p_node_index-axis_offset)||p_psi_N(other_axis,p_node_index)) psi_N(other_axis,face_index)=true;}}

    // set slip boundary conditions for tangential walls
    if(p_psi_R.array.base_pointer)
        for(int other_axis=0;other_axis<TV::m;other_axis++) if(other_axis!=axis){
            for(FACE_ITERATOR<TV> iterator(face_grid,0,GRID<TV>::BOUNDARY_REGION,-1,other_axis);iterator.Valid();iterator.Next()){TV_INT face_index=iterator.Face_Index(),p_node_index=face_index;
                T a=p_psi_R(other_axis,p_node_index-axis_offset),b=p_psi_R(other_axis,p_node_index);
                if(a || b) psi_R(other_axis,face_index)=(T).5*(a+b);}}
}
//#####################################################################
namespace PhysBAM{
template class IMPLICIT_VISCOSITY_UNIFORM<VECTOR<float,1> >;
template class IMPLICIT_VISCOSITY_UNIFORM<VECTOR<float,2> >;
template class IMPLICIT_VISCOSITY_UNIFORM<VECTOR<float,3> >;
template class IMPLICIT_VISCOSITY_UNIFORM<VECTOR<double,1> >;
template class IMPLICIT_VISCOSITY_UNIFORM<VECTOR<double,2> >;
template class IMPLICIT_VISCOSITY_UNIFORM<VECTOR<double,3> >;
}
