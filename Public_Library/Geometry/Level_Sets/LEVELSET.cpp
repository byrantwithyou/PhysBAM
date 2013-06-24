//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Sergey Koltakov, Neil Molino, Igor Neverov, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
template<class TV> typename LEVELSET<TV>::T_LINEAR_INTERPOLATION_SCALAR LEVELSET<TV>::interpolation_default;
template<class TV> typename LEVELSET<TV>::T_LINEAR_INTERPOLATION_VECTOR LEVELSET<TV>::normal_interpolation_default;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVELSET<TV>::
LEVELSET(GRID<TV>& grid_input,T_ARRAYS_SCALAR& phi_input,const int number_of_ghost_cells_input)
    :levelset_callbacks(0),collision_body_list(0),face_velocities_valid_mask_current(0),clamp_phi_with_collision_bodies(true),boundary_default(*new BOUNDARY<TV,T>),
    collision_aware_interpolation_plus(0),collision_aware_interpolation_minus(0),collision_unaware_interpolation(0),collidable_phi_replacement_value((T)1e-5),grid(grid_input),
    phi(phi_input),normals(0),curvature(0),cell_range(0),thread_queue(0),number_of_ghost_cells(number_of_ghost_cells_input)
{
    Set_Small_Number();
    Set_Max_Time_Step();
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
template<class TV> LEVELSET<TV>::
~LEVELSET()
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
template<class TV> void LEVELSET<TV>::
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
template<class TV> typename TV::SCALAR LEVELSET<TV>::
Collision_Aware_Phi(const TV& location) const
{
    assert(collision_body_list);
    return LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<GRID<TV>,T>(*collision_body_list,&valid_mask_current,collidable_phi_replacement_value).Clamped_To_Array(grid,phi,location);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET<TV>::
CFL(const T_FACE_ARRAYS_SCALAR& face_velocities) const
{
    T dt_convection=0;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        T local_V_norm=0;
        for(int axis=0;axis<TV::dimension;axis++)
            local_V_norm+=grid.one_over_dX[axis]*maxabs(face_velocities(axis,grid.First_Face_Index_In_Cell(axis,cell)),face_velocities(axis,grid.Second_Face_Index_In_Cell(axis,cell)));
        dt_convection=max(dt_convection,local_V_norm);}
    return 1/max(dt_convection,1/max_time_step);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET<TV>::
CFL(const ARRAY<TV,TV_INT>& velocity) const
{
    T dt_convection=0;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        dt_convection=max(dt_convection,TV::Dot_Product(velocity(cell),grid.one_over_dX));}
    return 1/max(dt_convection,1/max_time_step);
}
//#####################################################################
// Function Iterative_Find_Interface
//#####################################################################
template<class TV> TV LEVELSET<TV>::
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
template<class TV> void LEVELSET<TV>::
Compute_Gradient(ARRAY<TV,TV_INT>& gradient,const T time) const
{
    TV one_over_two_dx=(T).5*grid.one_over_dX;
    int ghost_cells=3;
    T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    gradient.Resize(grid.Domain_Indices(ghost_cells-1));
    for(CELL_ITERATOR<TV> iterator(grid,ghost_cells-1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        for(int axis=0;axis<TV::dimension;axis++){
            TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            gradient(cell_index)(axis)=(phi_ghost(cell_index+axis_vector)-phi_ghost(cell_index-axis_vector))*one_over_two_dx(axis);}}
}
//#####################################################################
// Function Gradient
//#####################################################################
template<class TV> TV LEVELSET<TV>::
Gradient(const TV& location) const
{
    TV G;
    for(int d=0;d<TV::m;d++){
        TV a(location),b(location);
        a(d)-=grid.dX(d);
        b(d)+=grid.dX(d);
        G(d)=(T).5*(Phi(b)-Phi(a))*grid.one_over_dX(d);}
    return G;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV LEVELSET<TV>::
Normal(const TV& location) const
{
    if(normals) return normal_interpolation->Clamped_To_Array(grid,*normals,location).Normalized();
    return Gradient(location).Normalized();
}
//#####################################################################
// Function Extended_Normal
//#####################################################################
template<class TV> TV LEVELSET<TV>::
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
template<class TV> void LEVELSET<TV>::
Compute_Normals(const T time)
{
    int ghost_cells=3;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells));
    boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);

    if(!normals) normals=new ARRAY<TV,TV_INT>(grid.Domain_Indices(ghost_cells-1));
    for(CELL_ITERATOR<TV> iterator(grid,ghost_cells-1);iterator.Valid();iterator.Next()){
        const TV_INT& cell = iterator.Cell_Index();
        int index=phi_ghost.Standard_Index(iterator.Cell_Index());
        TV& N((*normals)(cell));
        for(int d=0;d<TV::m;d++)
            N(d)=(T).5*(phi_ghost.array(index+phi_ghost.stride(d))-phi_ghost.array(index-phi_ghost.stride(d)))*grid.one_over_dX(d);
        N.Normalize();}
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> LEVELSET<TV>::
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
template<class TV> void LEVELSET<TV>::
Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists)
{
    if(!recompute_if_exists && cell_range) return;
    RANGE<TV_INT> pruned_range(phi.domain);
    pruned_range.max_corner-=1;
    if(!cell_range) cell_range=new ARRAY<INTERVAL<T>,TV_INT>(pruned_range);
    for(RANGE_ITERATOR<TV::m> it(pruned_range);it.Valid();it.Next()){
        int index=phi.Standard_Index(it.index);
        INTERVAL<T> interval(phi.array(index));
        for(int i=1;i<(1<<TV::m);i++){
            int ind=index;
            for(int d=0;d<TV::m;d++)
                ind+=((i>>d)&1)*phi.stride(d);
            interval.Enlarge_To_Include_Point(phi.array(ind));}
        (*cell_range)(it.index)=interval;}
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> LEVELSET<TV>::
Hessian(const ARRAY<T,TV_INT>& phi_input,const TV_INT& index) const
{
    int s_index=phi_input.Standard_Index(index);
    SYMMETRIC_MATRIX<T,TV::m> ddphi;
    for(int a=0;a<TV::m;a++)
        for(int b=0;b<a;b++)
            ddphi(a,b)=(T).25*(phi_input.array(s_index+phi_input.stride(a)+phi_input.stride(b))-phi_input.array(s_index+phi_input.stride(a)-phi_input.stride(b))
                -phi_input.array(s_index-phi_input.stride(a)+phi_input.stride(b))+phi_input.array(s_index-phi_input.stride(a)-phi_input.stride(b)))*grid.one_over_dX(a)*grid.one_over_dX(b);

    for(int a=0;a<TV::m;a++)
        ddphi(a,a)=(phi_input.array(s_index+phi_input.stride(a))-2*phi_input.array(s_index)+phi_input.array(s_index-phi_input.stride(a)))*sqr(grid.one_over_dX(a));
    return ddphi;
}
//#####################################################################
// Function Gradient
//#####################################################################
template<class TV> TV LEVELSET<TV>::
Gradient(const ARRAY<T,TV_INT>& phi_input,const TV_INT& index) const
{
    int s_index=phi_input.Standard_Index(index);
    TV dphi;
    for(int a=0;a<TV::m;a++)
        dphi(a)=(T).5*(phi_input.array(s_index+phi_input.stride(a))-phi_input.array(s_index-phi_input.stride(a)))*grid.one_over_dX(a);
    return dphi;
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET<TV>::
Compute_Curvature(const ARRAY<T,TV_INT>& phi_input,const TV_INT& index) const
{
    TV dphi=Gradient(phi_input,index);
    SYMMETRIC_MATRIX<T,TV::m> ddphi=Hessian(phi_input,index);

    T norm_squared=dphi.Magnitude_Squared(),norm=sqrt(norm_squared);
    if(norm<small_number) return LEVELSET_UTILITIES<T>::Sign(phi_input(index))/grid.dX.Min();
    T curvature=(dphi.Dot(ddphi*dphi)-norm_squared*ddphi.Trace())/(norm*norm_squared);
    return minmag(curvature,LEVELSET_UTILITIES<T>::Sign(curvature)/grid.dX.Min());
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
// kappa = - DIV(normal), negative for negative phi inside, positive for positive phi inside, sqrt(phix^2+phiy^2+phiy^2)=1 for distance functions
template<class TV> void LEVELSET<TV>::
Compute_Curvature(const T time)
{
    int ghost_cells=3;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);

    if(!curvature) curvature=new ARRAY<T,TV_INT>(grid.Domain_Indices(ghost_cells-1));
    for(CELL_ITERATOR<TV> iterator(grid,ghost_cells-1);iterator.Valid();iterator.Next())
        (*curvature)(iterator.Cell_Index())=Compute_Curvature(phi_ghost,iterator.Cell_Index());
}
//#####################################################################
// Function Compute_Curvature
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET<TV>::
Compute_Curvature(const TV& location) const
{
    TV l2=(location-(T).5*grid.dX-grid.domain.min_corner)*grid.one_over_dX;
    TV_INT cell(floor(l2));
    TV w(l2-TV(cell));

    T k[1<<TV::m];
    for(int i=0;i<(1<<TV::m);i++){
        TV_INT index(cell);
        for(int d=0;d<TV::m;d++)
            if(i&(1<<d))
                index(d)++;
        k[i]=Compute_Curvature(phi,index);}

    return LINEAR_INTERPOLATION<T,T>::Linear(k,w);
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
template<class TV> void LEVELSET<TV>::
Fast_Marching_Method(const T time,const T stopping_distance,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells,int process_sign)
{
    Get_Signed_Distance_Using_FMM(phi,time,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells,process_sign);
}
//#####################################################################
// Function Get_Signed_Distance_Using_FMM
//#####################################################################
template<class TV> void LEVELSET<TV>::
Get_Signed_Distance_Using_FMM(ARRAY<T,TV_INT>& signed_distance,const T time,const T stopping_distance,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells,int process_sign)
{
    const int ghost_cells=2*number_of_ghost_cells+1;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells),false);boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    FAST_MARCHING_METHOD_UNIFORM<GRID<TV> > fmm(static_cast<LEVELSET<TV>&>(*this),ghost_cells,thread_queue);
    fmm.Fast_Marching_Method(phi_ghost,stopping_distance,seed_indices,add_seed_indices_for_ghost_cells,process_sign);
    ARRAY<T,TV_INT>::Get(signed_distance,phi_ghost);
    boundary->Apply_Boundary_Condition(grid,signed_distance,time);
}
//#####################################################################
// Function Approximate_Surface_Size
//#####################################################################
// calculates the approximate perimeter using delta functions
template<class TV> typename TV::SCALAR LEVELSET<TV>::
Approximate_Surface_Size(const T interface_thickness,const T time) const
{
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(number_of_ghost_cells));
    boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,number_of_ghost_cells);
    T interface_half_width=interface_thickness*grid.dX.Max()/2,surface_size=0;
    for(RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next())
        surface_size+=LEVELSET_UTILITIES<T>::Delta(phi_ghost(it.index),interface_half_width)*Gradient(phi,it.index).Magnitude();
    return surface_size*grid.dX.Size();
}
template<class T> inline VECTOR<T,0> Principal_Curvatures_Helper(const VECTOR<T,1>& N,T norm,const SYMMETRIC_MATRIX<T,1>& H)
{
    return VECTOR<T,0>();
}
template<class T> inline VECTOR<T,1> Principal_Curvatures_Helper(const VECTOR<T,2>& N,T norm,const SYMMETRIC_MATRIX<T,2>& H)
{
    VECTOR<T,2> tangent=N.Perpendicular();
    return VECTOR<T,1>(tangent.Dot(H*tangent)/norm);
}
template<class T> inline VECTOR<T,2> Principal_Curvatures_Helper(const VECTOR<T,3>& N,T norm,const SYMMETRIC_MATRIX<T,3>& H)
{
    SYMMETRIC_MATRIX<T,3> P=(T)1-SYMMETRIC_MATRIX<T,3>::Outer_Product(N),M=SYMMETRIC_MATRIX<T,3>::Conjugate(P,H)/norm;
    T trace=M.Trace();
    QUADRATIC<T> quadratic(-1,trace,sqr(M(1,1))-M(0,1)*M(1,2)+sqr(M(2,1))-M(0,1)*M(2,3)+sqr(M(2,2))-M(1,2)*M(2,3));
    quadratic.Compute_Roots();
    if(quadratic.roots == 0) (T).5*VECTOR<T,2>(trace,trace);
    else if(quadratic.roots == 1) return VECTOR<T,2>(quadratic.root1,quadratic.root1);
    return VECTOR<T,2>(quadratic.root1,quadratic.root2);
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class TV> VECTOR<typename TV::SCALAR,TV::m-1> LEVELSET<TV>::
Principal_Curvatures(const TV& X) const
{
    TV N=Gradient(X);
    T grad_phi_magnitude=N.Normalize();
    return Principal_Curvatures_Helper(N,grad_phi_magnitude,Hessian(X));
}
namespace PhysBAM{
template class LEVELSET<VECTOR<float,1> >;
template class LEVELSET<VECTOR<float,2> >;
template class LEVELSET<VECTOR<float,3> >;
template class LEVELSET<VECTOR<double,1> >;
template class LEVELSET<VECTOR<double,2> >;
template class LEVELSET<VECTOR<double,3> >;
}
