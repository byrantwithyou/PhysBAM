//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP
//#####################################################################
#include <Core/Matrices/MATRIX.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
#include <Compressible/Euler_Equations/EULER_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>::
BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP(EULER_UNIFORM<TV>* euler_input,const T_FACE_VECTOR rho_far_field,
    const T_FACE_VECTOR p_far_field,const TV_FACE_VECTOR velocity_far_field,const T inflow_attenuation_input,
    const TV_SIDES& constant_extrapolation,const bool always_attenuate_input,const T_FACE_VECTOR linear_attenuations_input,
    const T_FACE_VECTOR_BOOL linear_attenuation_faces_input)
{
    euler=euler_input;
    Set_Constant_Extrapolation(constant_extrapolation);
    inflow_attenuation=inflow_attenuation_input;
    always_attenuate=always_attenuate_input;
    linear_attenuations=linear_attenuations_input;
    linear_attenuation_faces=linear_attenuation_faces_input;

    T gamma=dynamic_cast<EOS_GAMMA<T>*>(euler->eos)->gamma;
    for(int side=0;side<2*TV::m;side++){
        T e_far_field=euler->eos->e_From_p_And_rho(p_far_field(side),rho_far_field(side));
        T c_far_field=euler->eos->c(rho_far_field(side),e_far_field);
        S_far_field(side)=euler->eos->S(rho_far_field(side),e_far_field);
        iL_far_field(side)=-velocity_far_field(side)[0]+2*c_far_field/(gamma-1);
        iR_far_field(side)=velocity_far_field(side)[0]+2*c_far_field/(gamma-1);    
        U_far_field(side)=EULER<TV>::Get_Euler_State_From_rho_velocity_And_internal_energy(rho_far_field(side),velocity_far_field(side),e_far_field);}
}
//#####################################################################
// Function Attenuate_To_Far_Field_Values_Using_Riemann_Invariants
//#####################################################################
template<class TV> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>::
Attenuate_To_Far_Field_Values_Using_Riemann_Invariants(const T_ARRAYS_DIMENSION_BASE& u_ghost,const TV_INT& node_index,const int side,TV_DIMENSION &U,const T dt) const
{
    // TODO: Implement for multi-dimension.
    T gamma=dynamic_cast<EOS_GAMMA<T>*>(euler->eos)->gamma;
    T net_inflow_attenuation=exp((inflow_attenuation-1)*dt);

    T rho=u_ghost(node_index)(0);
    T u_velocity=EULER<TV>::Get_Velocity_Component(u_ghost(node_index),1);
    T e=EULER<TV>::e(u_ghost(node_index));

    T c=euler->eos->c(rho,e);
    T S=euler->eos->S(rho,e);
    T iL=-u_velocity+2*c/(gamma-1);
    T iR=u_velocity+2*c/(gamma-1);

    int flip=side&1?-1:1;
    if(flip*u_velocity > 0) S=S_far_field(side)+net_inflow_attenuation*(S-S_far_field(side));
    if(flip*(u_velocity-c) > 0) iL=iL_far_field(side)+net_inflow_attenuation*(iL-iL_far_field(side));
    if(flip*(u_velocity+c) > 0) iR=iR_far_field(side)+net_inflow_attenuation*(iR-iR_far_field(side));

    u_velocity=(iR-iL)/2;
    c=(T)((iR-u_velocity)*(gamma-1)*.5);
    rho=pow(sqr(c)/(gamma*S),1/(gamma-1));
    e=euler->eos->e_From_S_And_rho(S,rho);

    U[0]=rho;
    U[1]=rho*u_velocity;
    U[2]=rho*(e+sqr(u_velocity)/2);
}
//#####################################################################
// Function Attenuate_To_Far_Field_Values_Using_Characteristics
//#####################################################################
template<class TV> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>::
Attenuate_To_Far_Field_Values_Using_Characteristics(const T_ARRAYS_DIMENSION_BASE& u_ghost,const TV_INT& node_index,const int side,TV_DIMENSION &U,const T dt) const
{
    T net_inflow_attenuation=exp((inflow_attenuation-1)*dt);
    if(linear_attenuation_faces[side]){net_inflow_attenuation=(1-linear_attenuations(side));}

    if(always_attenuate){
        U=u_ghost(node_index);
        U=U_far_field(side)+net_inflow_attenuation*(U-U_far_field(side));
        return;}
   
    int axis=side/2;
    MATRIX<T,d,d> L,R;
    TV_DIMENSION LU,LU_far_field;
    VECTOR<T,d> lambda,lambda_left,lambda_right;

    // extract a 1d array
    int min_index=u_ghost.Domain_Indices().min_corner[axis],max_index=u_ghost.Domain_Indices().max_corner[axis];
    TV_INT current_index=node_index;
    ARRAY<TV_DIMENSION,VECTOR<int,1> > U_1d(min_index,max_index);
    for(int i=min_index;i<max_index;i++){current_index[axis]=i;U_1d(i)=u_ghost(current_index);}

    // apply L
    euler->eigensystems_default[axis]->Eigenvectors(U_1d,node_index[axis],L,R);
    euler->eigensystems_default[axis]->Eigenvalues(U_1d,node_index[axis],lambda,lambda_left,lambda_right);
    LU=L*u_ghost(node_index);
    LU_far_field=L*U_far_field(side);

    // attenuate
    int flip=side&1?-1:1;
    for(int k=0;k<d;k++) if(flip*lambda_left(k)>0) LU(k)=LU_far_field(k)+net_inflow_attenuation*(LU(k)-LU_far_field(k));

    // apply R
    U=R.Transpose_Times(LU);
}
//#####################################################################
// Function Attenuate_To_Far_Field_Values
//#####################################################################
template<class TV> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>::
Attenuate_To_Far_Field_Values(const T_ARRAYS_DIMENSION_BASE& u_ghost,const TV_INT& node_index,const int side,TV_DIMENSION &U,const T dt) const
{
    U=u_ghost(node_index);
    T net_inflow_attenuation=exp((inflow_attenuation-1)*dt);
    U=U_far_field(side)+net_inflow_attenuation*(U-U_far_field(side));
}
//#####################################################################
// Function Fill_Single_Ghost_Region
//#####################################################################
template<class TV> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>::
Fill_Single_Ghost_Region(const GRID<TV>& grid,T_ARRAYS_DIMENSION_BASE& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells) const
{
    int total_number_of_ghost_cells=u_ghost.Domain_Indices().max_corner.x-grid.Domain_Indices().max_corner.x;
    int ghost_cells_to_fill=number_of_ghost_cells;
    bool attenuate_dirichlet_cells=true;
    int axis=side/2;
    if(Constant_Extrapolation(side)){
        RANGE<TV_INT> region_new(region);
        if(!attenuate_dirichlet_cells){
            int already_filled_cells=total_number_of_ghost_cells-ghost_cells_to_fill;
            side&1?region_new.min_corner(axis)+=already_filled_cells:region_new.max_corner(axis)-=already_filled_cells;}
        int boundary=Boundary(side,region_new);
        if(attenuate_dirichlet_cells){
            Fill_Single_Ghost_Region(grid,u_ghost,side,region_new); // extrapolate values first, so that eigenvectors calculated at boundary flux values are calculated correctly.
            for(CELL_ITERATOR<TV> iterator(grid,region_new);iterator.Valid();iterator.Next()){
                TV_INT cell=iterator.Cell_Index(),boundary_node=cell;boundary_node[axis]=boundary;TV_DIMENSION U_boundary;
                Attenuate_To_Far_Field_Values_Using_Characteristics(u_ghost,boundary_node,side,U_boundary,dt);
                u_ghost(cell)=U_boundary;}}
        else{
            for(CELL_ITERATOR<TV> iterator(grid,region_new);iterator.Valid();iterator.Next()){
                TV_INT cell=iterator.Cell_Index(),boundary_node=cell;boundary_node[axis]=boundary;
                u_ghost(cell)=u_ghost(boundary_node);}}}
    else{
        int boundary=Boundary(side,region);
        int reflection_times_two;
        if(grid.Is_MAC_Grid())
            reflection_times_two=2*boundary+(side&1?1:-1);
        else reflection_times_two=2*boundary+2;
        for(CELL_ITERATOR<TV> iterator(grid,region);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            TV_INT reflected_node=cell;reflected_node[axis]=reflection_times_two-cell[axis];
            T rho=u_ghost(reflected_node)(0);
            TV velocity=EULER<TV>::Get_Velocity(u_ghost(reflected_node));
            velocity(axis)*=-1;
            T e=EULER<TV>::e(u_ghost(reflected_node));
            EULER<TV>::Set_Euler_State_From_rho_velocity_And_internal_energy(u_ghost,cell,rho,velocity,e);}}
}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>::
Fill_Ghost_Cells(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_BASE& u,T_ARRAYS_DIMENSION_BASE& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    T_ARRAYS_DIMENSION_BASE::Put(u,u_ghost); // interior
    VECTOR<RANGE<TV_INT>,2*TV::m> regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int side=0;side<2*TV::m;side++){
        Fill_Single_Ghost_Region(grid,u_ghost,regions(side),side,dt,time,number_of_ghost_cells);}
}
//#####################################################################
// Function Apply_Boundary_Condition_Single_Side
//#####################################################################
template<class TV> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>::
Apply_Boundary_Condition_Single_Side(const GRID<TV>& grid,T_ARRAYS_DIMENSION_BASE& u,const int side,const T time) const 
{
    if(grid.Is_MAC_Grid()) return;
    int axis=side/2;
    int axis_side=side%2;
    if(!euler->mpi_grid||!euler->mpi_grid->Neighbor(axis,axis_side)){
        if(!Constant_Extrapolation(side)) for(CELL_ITERATOR<TV> iterator(grid,0,GRID<TV>::BOUNDARY_INTERIOR_REGION,side);iterator.Valid();iterator.Next()){
            TV_INT boundary_node=iterator.Cell_Index();
            T rho=u(boundary_node)(0);
            TV velocity=EULER<TV>::Get_Velocity(u(boundary_node));
            T e=EULER<TV>::e(u(boundary_node));
            velocity[axis]=0; // Boundary condition
            EULER<TV>::Set_Euler_State_From_rho_velocity_And_internal_energy(u,boundary_node,rho,velocity,e);}}
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<TV>::
Apply_Boundary_Condition(const GRID<TV>& grid,T_ARRAYS_DIMENSION_BASE& u,const T time) const
{
    if(grid.Is_MAC_Grid()) return;
    for(int side=0;side<2*TV::m;side++) Apply_Boundary_Condition_Single_Side(grid,u,side,time);
}
//#####################################################################
namespace PhysBAM{
template class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<VECTOR<float,1> >;
template class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<VECTOR<float,2> >;
template class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<VECTOR<float,3> >;
template class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<VECTOR<double,1> >;
template class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<VECTOR<double,2> >;
template class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<VECTOR<double,3> >;
}
