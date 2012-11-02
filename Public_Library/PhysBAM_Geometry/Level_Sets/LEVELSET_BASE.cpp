//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Sergey Koltakov, Neil Molino, Igor Neverov, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_BASE.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
template<class TV> typename LEVELSET_BASE<TV>::T_LINEAR_INTERPOLATION_SCALAR LEVELSET_BASE<TV>::interpolation_default;
template<class TV> typename LEVELSET_BASE<TV>::T_LINEAR_INTERPOLATION_VECTOR LEVELSET_BASE<TV>::normal_interpolation_default;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVELSET_BASE<TV>::
LEVELSET_BASE(GRID<TV>& grid_input,T_ARRAYS_SCALAR& phi_input,const int number_of_ghost_cells_input)
    :levelset_callbacks(0),collision_body_list(0),face_velocities_valid_mask_current(0),clamp_phi_with_collision_bodies(true),boundary_default(*new BOUNDARY_UNIFORM<GRID<TV>,T>),
    collision_aware_interpolation_plus(0),collision_aware_interpolation_minus(0),collision_unaware_interpolation(0),collidable_phi_replacement_value((T)1e-5),grid(grid_input),
    phi(phi_input),normals(0),curvature(0),cell_range(0),thread_queue(0),number_of_ghost_cells(number_of_ghost_cells_input)
{
    Set_Small_Number();
    Set_Max_Time_Step();
    curvature_motion=false; // default is no curvature motion
    Initialize_FMM_Initialization_Iterative_Solver();
    boundary=&boundary_default;
    interpolation=&interpolation_default;
    curvature_interpolation=&interpolation_default;
    normal_interpolation=&normal_interpolation_default;
    secondary_interpolation=0;
    Set_Band_Width();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_BASE<TV>::
~LEVELSET_BASE()
{
    assert(!collision_unaware_interpolation);
    delete collision_aware_interpolation_plus;
    delete collision_aware_interpolation_minus;
    delete &boundary_default;
    delete normals;
    delete curvature;
    delete cell_range;
}
//#####################################################################
// Function Set_Collision_Body_List
//#####################################################################
template<class TV> void LEVELSET_BASE<TV>::
Set_Collision_Body_List(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collision_body_list_input,const bool set_secondary_interpolation)
{
    collision_body_list=&collision_body_list_input;
    delete collision_aware_interpolation_plus;delete collision_aware_interpolation_minus;
    collision_aware_interpolation_plus=new T_LINEAR_INTERPOLATION_SCALAR;
    collision_aware_interpolation_minus=new LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<TV>,T>(*collision_body_list,&valid_mask_current,collidable_phi_replacement_value);
    if(set_secondary_interpolation) secondary_interpolation=collision_aware_interpolation_minus;
}
//#####################################################################
// Function Collision_Aware_Phi
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_BASE<TV>::
Collision_Aware_Phi(const TV& location) const
{
    assert(collision_body_list);
    return LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<TV>,T>(*collision_body_list,&valid_mask_current,collidable_phi_replacement_value).Clamped_To_Array(grid,phi,location);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_BASE<TV>::
CFL(const T_FACE_ARRAYS_SCALAR& face_velocities) const
{
    T dt_convection=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        T local_V_norm=0;
        for(int axis=0;axis<TV::dimension;axis++)
            local_V_norm+=grid.one_over_dX[axis]*maxabs(face_velocities(axis,grid.First_Face_Index_In_Cell(axis,cell)),face_velocities(axis,grid.Second_Face_Index_In_Cell(axis,cell)));
        dt_convection=max(dt_convection,local_V_norm);}
    T dt_curvature=(curvature_motion && TV::dimension>1)?sigma*2*grid.one_over_dX.Magnitude_Squared():0;
    T dt_overall=dt_convection+dt_curvature;
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_BASE<TV>::
CFL(const ARRAY<TV,TV_INT>& velocity) const
{
    T dt_convection=0;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        dt_convection=max(dt_convection,TV::Dot_Product(velocity(cell),grid.one_over_dX));}
    T dt_curvature=(curvature_motion && TV::dimension>1)?sigma*2*grid.one_over_dX.Magnitude_Squared():0;
    T dt_overall=dt_convection+dt_curvature;
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function Iterative_Find_Interface
//#####################################################################
template<class TV> TV LEVELSET_BASE<TV>::
Iterative_Find_Interface(TV left,TV right,const int iterations) const
{
    T phi_left=Phi(left),phi_right=Phi(right);
    TV interface=LINEAR_INTERPOLATION<T,TV>::Linear(left,right,LEVELSET_UTILITIES<T>::Theta(phi_left,phi_right));
    int phi_left_sign=(phi_left<=0?-1:1),phi_right_sign=(phi_right<=0?-1:1);
    for(int i=0;i<iterations;i++){
        T phi=Phi(interface);
        int phi_sign=(phi<=0?-1:1);
        if(phi_left_sign*phi_sign<0){
            right=interface;phi_right=phi;phi_right_sign=phi_sign;
            interface=LINEAR_INTERPOLATION<T,TV>::Linear(left,interface,LEVELSET_UTILITIES<T>::Theta(phi_left,phi));}
        else if(phi_right_sign*phi_sign<0){
            left=interface;phi_left=phi;phi_left_sign=phi_sign;
            interface=LINEAR_INTERPOLATION<T,TV>::Linear(interface,right,LEVELSET_UTILITIES<T>::Theta(phi,phi_right));}
        else break;}
    return interface;
}
//#####################################################################
// Function Compute_Gradient
//#####################################################################
template<class TV> void LEVELSET_BASE<TV>::
Compute_Gradient(ARRAY<TV,TV_INT>& gradient,const T time) const
{
    TV one_over_two_dx=(T).5*grid.one_over_dX;
    int ghost_cells=3;
    T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    gradient.Resize(grid.Domain_Indices(ghost_cells-1));
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,ghost_cells-1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        for(int axis=0;axis<TV::dimension;axis++){
            TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            gradient(cell_index)(axis)=(phi_ghost(cell_index+axis_vector)-phi_ghost(cell_index-axis_vector))*one_over_two_dx(axis);}}
}
//#####################################################################
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV LEVELSET_BASE<TV>::
Normal(const TV& location) const
{
    if(normals) return normal_interpolation->Clamped_To_Array(grid,*normals,location).Normalized();
    TV N;
    for(int d=0;d<TV::m;d++){
        TV a(location),b(location);
        a(d)-=grid.dX(d);
        b(d)+=grid.dX(d);
        N(d)=(T).5*(Phi(b)-Phi(a))*grid.one_over_dX(d);}
    return N.Normalized();
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class TV> TV LEVELSET_BASE<TV>::
Extended_Normal(const TV& location) const
{
    if(normals) return normal_interpolation->Clamped_To_Array(grid,*normals,grid.Clamp(location)).Normalized();
    TV N;
    for(int d=0;d<TV::m;d++){
        TV a(location),b(location);
        a(d)-=grid.dX(d);
        b(d)+=grid.dX(d);
        N(d)=(T).5*(Extended_Phi(b)-Extended_Phi(a))*grid.one_over_dX(d);}
    return N.Normalized();
}
//#####################################################################
// Function Compute_Normals
//#####################################################################
// note that sqrt(phix^2+phiy^2+phiz^2)=1 if it's a distance function
template<class TV> void LEVELSET_BASE<TV>::
Compute_Normals(const T time)
{
    int ghost_cells=3;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells));
    boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    TV_INT offset;
    offset(TV::m-1)=1;
    for(int d=TV::m-1;d>=1;d--)
        offset(d-1)=offset(d)*phi_ghost.counts(d);

    if(!normals) normals=new ARRAY<TV,TV_INT>(grid.Domain_Indices(ghost_cells-1));
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,ghost_cells-1);iterator.Valid();iterator.Next()){
        const TV_INT& cell = iterator.Cell_Index();
        int index=phi_ghost.Standard_Index(iterator.Cell_Index());
        TV& N((*normals)(cell));
        for(int d=0;d<TV::m;d++)
            N(d)=(T).5*(phi_ghost.array(index+offset(d))-phi_ghost.array(index-offset(d)))*grid.one_over_dX(d);
        N.Normalize();}
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> LEVELSET_BASE<TV>::
Hessian(const TV& X) const
{
    SYMMETRIC_MATRIX<T,TV::m> H;
    for(int i=0;i<TV::m;i++){
        TV A(X),B(X);
        A(i)-=grid.dX(i);
        B(i)+=grid.dX(i);
        H(i,i)=(Phi(B)-2*Phi(X)+Phi(A))*sqr(grid.one_over_dX(i));}

    for(int i=0;i<TV::m;i++){
        TV A(X),B(X);
        A(i)-=grid.dX(i);
        B(i)+=grid.dX(i);
        for(int j=i+1;j<TV::m;j++){
            TV C(A),D(A),E(B),F(B);
            C(j)-=grid.dX(j);
            D(j)+=grid.dX(j);
            E(j)-=grid.dX(j);
            F(j)+=grid.dX(j);
            H(i,j)=(Phi(F)-Phi(E)-Phi(D)+Phi(A))*(T).25*grid.one_over_dX(i)*grid.one_over_dX(j);}}

    return H;
}
//#####################################################################
// Function Compute_Cell_Minimum_And_Maximum
//#####################################################################
template<class TV> void LEVELSET_BASE<TV>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    if(!recompute_if_exists && cell_range) return;
    TV_INT strides=phi.Strides();
    RANGE<TV_INT> pruned_range(phi.domain);
    pruned_range.max_corner-=1;
    if(!cell_range) cell_range=new ARRAY<INTERVAL<T>,TV_INT>(pruned_range);
    for(RANGE_ITERATOR<TV::m> it(pruned_range);it.Valid();it.Next()){
        int index=phi.Standard_Index(it.index);
        INTERVAL<T> interval(phi.array(index));
        for(int i=1;i<(1<<TV::m);i++){
            int ind=index;
            for(int d=0;d<TV::m;d++)
                ind+=((i>>d)&1)*strides(d);
            interval.Enlarge_To_Include_Point(phi.array(ind));}
        (*cell_range)(it.index)=interval;}
}
template class LEVELSET_BASE<VECTOR<float,1> >;
template class LEVELSET_BASE<VECTOR<float,2> >;
template class LEVELSET_BASE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_BASE<VECTOR<double,1> >;
template class LEVELSET_BASE<VECTOR<double,2> >;
template class LEVELSET_BASE<VECTOR<double,3> >;
#endif
