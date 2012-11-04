//#####################################################################
// Copyright 2009, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_MACCORMACK_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Geometry/Advection_Collidable/ADVECTION_WRAPPER_COLLIDABLE_CELL.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVELSET_ADVECTION<TV>::
LEVELSET_ADVECTION(LEVELSET<TV>* levelset)
    :levelset(levelset),advection(0),nested_semi_lagrangian_collidable(0),semi_lagrangian_collidable(0),advection_maccormack(0)
{
    Set_Reinitialization_Runge_Kutta_Order();
    Set_Reinitialization_CFL();
    Use_WENO_For_Reinitialization();
    Use_Level_Set_Advection_Method();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_ADVECTION<TV>::
~LEVELSET_ADVECTION()
{
    delete nested_semi_lagrangian_collidable;
    delete semi_lagrangian_collidable;
}
//#####################################################################
// Function Use_Semi_Lagrangian_Collidable_Advection
//#####################################################################
template<class TV> void LEVELSET_ADVECTION<TV>::
Use_Semi_Lagrangian_Collidable_Advection(const T_GRID_BASED_COLLISION_GEOMETRY& body_list,const T phi_replacement_value,const T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input)
{
    assert(!nested_semi_lagrangian_collidable&&!semi_lagrangian_collidable);
    nested_semi_lagrangian_collidable=new T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL(body_list,levelset->valid_mask_current,levelset->valid_mask_next,phi_replacement_value,true);
    semi_lagrangian_collidable=new ADVECTION_WRAPPER_COLLIDABLE_CELL<GRID<TV>,T,T_FACE_LOOKUP,T_ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL,T_FACE_LOOKUP_COLLIDABLE>(*nested_semi_lagrangian_collidable,body_list,
        face_velocities_valid_mask_input);
    Set_Custom_Advection(*semi_lagrangian_collidable);
}
//#####################################################################
// Function HJ_WENO
// phi is (-2,m_3), phix_minus and phix_plus are (1,m)
//#####################################################################
template<class TV> void LEVELSET_ADVECTION<TV>::
HJ_WENO(const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,ARRAY<T,VECTOR<int,1> >& phix_minus,ARRAY<T,VECTOR<int,1> >& phix_plus) const
{
    T epsilon=(T)1e-6; // works only because phi is a distance function
    T one_over_dx=1/dx;
    ARRAY<T,VECTOR<int,1> > D1(-2,m+2);for(int i=-3;i<m+2;i++) D1(i)=(phi(i+1)-phi(i))*one_over_dx; // 1st divided difference
    for(int i=0;i<m;i++){
        phix_minus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::WENO(D1(i-3),D1(i-2),D1(i-1),D1(i),D1(i+1),epsilon);
        phix_plus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::WENO(D1(i+2),D1(i+1),D1(i),D1(i-1),D1(i-2),epsilon);}
}
//#####################################################################
// Function HJ_ENO
// order = 1, 2 or 3, phi is (-2,m_3), phix_minus and phix_plus are (1,m)
//#####################################################################
template<class TV> void LEVELSET_ADVECTION<TV>::
HJ_ENO(const int order,const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,ARRAY<T,VECTOR<int,1> >& phix_minus,ARRAY<T,VECTOR<int,1> >& phix_plus) const
{
    T one_over_dx=1/dx,one_over_two_dx=(T).5*one_over_dx,one_over_three_dx=(T)one_third*one_over_dx;
    ARRAY<T,VECTOR<int,1> > D1(-2,m+2),D2(-2,m+1),D3(-2,m); // divided differences
    for(int i=-3;i<m+2;i++) D1(i)=(phi(i+1)-phi(i))*one_over_dx;
    if(order >= 2) for(int i=-3;i<m+1;i++) D2(i)=(D1(i+1)-D1(i))*one_over_two_dx;
    if(order == 3) for(int i=-3;i<m;i++) D3(i)=(D2(i+1)-D2(i))*one_over_three_dx;

    if(order == 1) for(int i=0;i<m;i++){phix_minus(i)=D1(i-1);phix_plus(i)=D1(i);}
    else if(order == 2) for(int i=0;i<m;i++){
        phix_minus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,D1(i-1),D2(i-2),D2(i-1));
        phix_plus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,D1(i),-D2(i),-D2(i-1));}
    else if(order == 3) for(int i=0;i<m;i++){
        phix_minus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,D1(i-1),D2(i-2),D2(i-1),D3(i-3),D3(i-2),D3(i-1));
        phix_plus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,D1(i),-D2(i),-D2(i-1),D3(i),D3(i-1),D3(i-2));}
}
//#####################################################################
// Function Use_Maccormack_Advection
//#####################################################################
// calculates the approximate area using Heaviside functions
template<class TV> void LEVELSET_ADVECTION<TV>::
Use_Maccormack_Advection(const ARRAY<bool,TV_INT>& cell_mask)
{
    advection_maccormack=new ADVECTION_MACCORMACK_UNIFORM<GRID<TV>,T,ADVECTION<GRID<TV>,T> >(*advection,0,&cell_mask,0);
    Set_Custom_Advection(*advection_maccormack);
}
//#####################################################################
// Function Approximate_Negative_Material
//#####################################################################
// calculates the approximate area using Heaviside functions
template<class TV> typename TV::SCALAR LEVELSET_ADVECTION<TV>::
Approximate_Negative_Material(const T interface_thickness,const T time) const
{
    GRID<TV>& grid=levelset->grid;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    GRID<TV> node_grid=grid.Is_MAC_Grid()?grid.Get_Regular_Grid_At_MAC_Positions():grid;
    T interface_half_width=interface_thickness*grid.dX.Max()/2,volume=0;
    for(UNIFORM_GRID_ITERATOR_NODE<TV> iterator(node_grid);iterator.Valid();iterator.Next()) volume+=LEVELSET_UTILITIES<T>::Heaviside(-phi(iterator.Node_Index()),interface_half_width);
    return volume*grid.Cell_Size();
}
//#####################################################################
// Function Approximate_Positive_Material
//#####################################################################
// calculates the approximate area using Heaviside functions
template<class TV> typename TV::SCALAR LEVELSET_ADVECTION<TV>::
Approximate_Positive_Material(const T interface_thickness,const T time) const
{
    GRID<TV>& grid=levelset->grid;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    GRID<TV> node_grid=grid.Is_MAC_Grid()?grid.Get_Regular_Grid_At_MAC_Positions():grid;
    T interface_half_width=interface_thickness*grid.dX.Max()/2,volume=0;
    for(UNIFORM_GRID_ITERATOR_NODE<TV> iterator(node_grid);iterator.Valid();iterator.Next()) volume+=LEVELSET_UTILITIES<T>::Heaviside(phi(iterator.Node_Index()),interface_half_width);
    return volume*grid.Cell_Size();
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class TV> void LEVELSET_ADVECTION<TV>::
Euler_Step(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocity,const T dt,const T time,const int number_of_ghost_cells)
{
    GRID<TV>& grid=levelset->grid;
    BOUNDARY_UNIFORM<GRID<TV>,T>* boundary=levelset->boundary;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    assert(grid.Is_MAC_Grid());
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(number_of_ghost_cells));
    boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells);
    advection->Update_Advection_Equation_Cell(grid,phi,phi_ghost,face_velocity,*boundary,dt,time);
    boundary->Apply_Boundary_Condition(grid,phi,time+dt);
}
//#####################################################################
// Function Local_WENO_Advect
//#####################################################################
// phi is (-2,m+3), u, distance and u_phix are [0,m)
template<class T> static void
Local_WENO_Advect(const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& u,const ARRAY<T,VECTOR<int,1> >& distance,ARRAY<T,VECTOR<int,1> >& u_phix,const T half_band_width)
{
    T epsilon=(T)1e-6*sqr(dx); // 1e-6 works since phi is a distance function - sqr(dx) since undivided differences are used
    T one_over_dx=1/dx;
    for(int i=0;i<m;i++) if(abs(distance(i)) <= half_band_width){ // one_over_dx since undivided differences are used
        if(u(i) > 0) u_phix(i)=u(i)*ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<T,1> >,T>::WENO(phi(i-2)-phi(i-3),phi(i-1)-phi(i-2),phi(i)-phi(i-1),phi(i+1)-phi(i),phi(i+2)-phi(i+1),epsilon)*one_over_dx;
        else u_phix(i)=u(i)*ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<T,1> >,T>::WENO(phi(i+3)-phi(i+2),phi(i+2)-phi(i+1),phi(i+1)-phi(i),phi(i)-phi(i-1),phi(i-1)-phi(i-2),epsilon)*one_over_dx;}
}
//#####################################################################
// Function Local_ENO_Advect
//#####################################################################
// order = 1, 2 or 3, phi is (-2,m_3), u, distance and u_phix are [0,m)
template<class T> static void
Local_ENO_Advect(const int order,const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& u,const ARRAY<T,VECTOR<int,1> >& distance,ARRAY<T,VECTOR<int,1> >& u_phix,const T half_band_width)
{
    T one_over_dx=1/dx;
    if(order == 1){for(int i=0;i<m;i++) if(abs(distance(i)) <= half_band_width){
        if(u(i) > 0) u_phix(i)=u(i)*(phi(i)-phi(i-1))*one_over_dx;
        else u_phix(i)=u(i)*(phi(i+1)-phi(i))*one_over_dx;}}
    else if(order == 2){for(int i=0;i<m;i++) if(abs(distance(i)) <= half_band_width){
        if(u(i) > 0) u_phix(i)=u(i)*(phi(i)-phi(i-1)+(T).5*minmag(phi(i)-2*phi(i-1)+phi(i-2),phi(i+1)-2*phi(i)+phi(i-1)))*one_over_dx;
        else u_phix(i)=u(i)*(phi(i+1)-phi(i)-(T).5*minmag(phi(i+2)-2*phi(i+1)+phi(i),phi(i+1)-2*phi(i)+phi(i-1)))*one_over_dx;}}
    else if(order == 3){for(int i=0;i<m;i++) if(abs(distance(i)) <= half_band_width){
        if(u(i) > 0){
            T D2_left=phi(i)-2*phi(i-1)+phi(i-2),D2_right=phi(i+1)-2*phi(i)+phi(i-1); // both scaled (multiplied by) 2*dx*dx
            if(abs(D2_left) <= abs(D2_right))
                u_phix(i)=u(i)*(phi(i)-phi(i-1)+(T).5*D2_left+(T)one_third*minmag(phi(i)-3*(phi(i-1)-phi(i-2))-phi(i-3),phi(i+1)-3*(phi(i)-phi(i-1))-phi(i-2)))*one_over_dx;
            else u_phix(i)=u(i)*(phi(i)-phi(i-1)+(T).5*D2_right-(T)one_sixth*minmag(phi(i+2)-3*(phi(i+1)-phi(i))-phi(i-1),phi(i+1)-3*(phi(i)-phi(i-1))-phi(i-2)))*one_over_dx;}
        else{
            T D2_left=phi(i+2)-2*phi(i+1)+phi(i),D2_right=phi(i+1)-2*phi(i)+phi(i-1); // both scaled (multiplied by) -2*dx*dx
            if(abs(D2_left) <= abs(D2_right))
                u_phix(i)=u(i)*(phi(i+1)-phi(i)-(T).5*D2_left+(T)one_third*minmag(phi(i+3)-3*(phi(i+2)-phi(i+1))+phi(i),phi(i+2)-3*(phi(i+1)-phi(i))+phi(i-1)))*one_over_dx;
            else u_phix(i)=u(i)*(phi(i+1)-phi(i)-(T).5*D2_right-(T)one_sixth*minmag(phi(i+1)-3*(phi(i)-phi(i-1))+phi(i-2),phi(i+2)-3*(phi(i+1)-phi(i))+phi(i-1)))*one_over_dx;}}}
}
//#####################################################################
// Function Euler_Step_High_Order_Helper
//#####################################################################
template<class T> static void
Euler_Step_High_Order_Helper(const GRID<VECTOR<T,1> >& grid,const ARRAY<VECTOR<T,1> ,VECTOR<int,1> >& V,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& phi_ghost,ARRAY<T,VECTOR<int,1> >& rhs,const int spatial_order,
    const T half_band_width)
{
    int m=grid.counts.x;T dx=grid.dX.x;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_1d_x(-ghost_cells,m+ghost_cells),u_1d(0,m),distance_1d_x(0,m);
    for(int i=-ghost_cells;i<m+ghost_cells;i++) phi_1d_x(i)=phi_ghost(i);
    for(int i=0;i<m;i++){u_1d(i)=V(i).x;distance_1d_x(i)=phi(i);}
    if(spatial_order == 5) Local_WENO_Advect(m,dx,phi_1d_x,u_1d,distance_1d_x,rhs,half_band_width);
    else Local_ENO_Advect(spatial_order,m,dx,phi_1d_x,u_1d,distance_1d_x,rhs,half_band_width);
}
template<class T> static void
Euler_Step_High_Order_Helper(const GRID<VECTOR<T,2> >& grid,const ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& V,const ARRAY<T,VECTOR<int,2> >& phi,const ARRAY<T,VECTOR<int,2> >& phi_ghost,ARRAY<T,VECTOR<int,2> >& rhs,const int spatial_order,
    const T half_band_width)
{
    int m=grid.counts.x,n=grid.counts.y;T dx=grid.dX.x,dy=grid.dX.y;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_1d_x(-ghost_cells,m+ghost_cells),u_1d(0,m),distance_1d_x(0,m),u_phix_1d(0,m);
    for(int j=0;j<n;j++){
        for(int i=-ghost_cells;i<m+ghost_cells;i++) phi_1d_x(i)=phi_ghost(i,j);
        for(int i=0;i<m;i++){u_1d(i)=V(i,j).x;distance_1d_x(i)=phi(i,j);}
        if(spatial_order == 5) Local_WENO_Advect(m,dx,phi_1d_x,u_1d,distance_1d_x,u_phix_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,m,dx,phi_1d_x,u_1d,distance_1d_x,u_phix_1d,half_band_width);
        for(int i=0;i<m;i++) rhs(i,j)=u_phix_1d(i);}
    ARRAY<T,VECTOR<int,1> > phi_1d_y(-ghost_cells,n+ghost_cells),v_1d(0,n),distance_1d_y(0,n),v_phiy_1d(0,n);
    for(int i=0;i<m;i++){
        for(int j=-ghost_cells;j<n+ghost_cells;j++) phi_1d_y(j)=phi_ghost(i,j);
        for(int j=0;j<n;j++){v_1d(j)=V(i,j).y;distance_1d_y(j)=phi(i,j);}
        if(spatial_order == 5) Local_WENO_Advect(n,dy,phi_1d_y,v_1d,distance_1d_y,v_phiy_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,n,dy,phi_1d_y,v_1d,distance_1d_y,v_phiy_1d,half_band_width);
        for(int j=0;j<n;j++) rhs(i,j)+=v_phiy_1d(j);}
}
template<class T> static void
Euler_Step_High_Order_Helper(const GRID<VECTOR<T,3> >& grid,const ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& V,const ARRAY<T,VECTOR<int,3> >& phi,const ARRAY<T,VECTOR<int,3> >& phi_ghost,ARRAY<T,VECTOR<int,3> >& rhs,const int spatial_order,
    const T half_band_width)
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;T dx=grid.dX.x,dy=grid.dX.y,dz=grid.dX.z;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_1d_x(-ghost_cells,m+ghost_cells),u_1d(0,m),distance_1d_x(0,m),u_phix_1d(0,m);
    for(int j=0;j<n;j++) for(int ij=0;ij<mn;ij++){
        for(int i=-ghost_cells;i<m+ghost_cells;i++) phi_1d_x(i)=phi_ghost(i,j,ij);
        for(int i=0;i<m;i++){u_1d(i)=V(i,j,ij).x;distance_1d_x(i)=phi(i,j,ij);}
        if(spatial_order == 5) Local_WENO_Advect(m,dx,phi_1d_x,u_1d,distance_1d_x,u_phix_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,m,dx,phi_1d_x,u_1d,distance_1d_x,u_phix_1d,half_band_width);
        for(int i=0;i<m;i++) rhs(i,j,ij)=u_phix_1d(i);}
    ARRAY<T,VECTOR<int,1> > phi_1d_y(-ghost_cells,n+ghost_cells),v_1d(0,n),distance_1d_y(0,n),v_phiy_1d(0,n);
    for(int i=0;i<m;i++) for(int ij=0;ij<mn;ij++){
        for(int j=-ghost_cells;j<n+ghost_cells;j++) phi_1d_y(j)=phi_ghost(i,j,ij);
        for(int j=0;j<n;j++){v_1d(j)=V(i,j,ij).y;distance_1d_y(j)=phi(i,j,ij);}
        if(spatial_order == 5) Local_WENO_Advect(n,dy,phi_1d_y,v_1d,distance_1d_y,v_phiy_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,n,dy,phi_1d_y,v_1d,distance_1d_y,v_phiy_1d,half_band_width);
        for(int j=0;j<n;j++) rhs(i,j,ij)+=v_phiy_1d(j);}
    ARRAY<T,VECTOR<int,1> > phi_1d_z(-ghost_cells,mn+ghost_cells),w_1d(0,mn),distance_1d_z(0,mn),w_phiz_1d(0,mn);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){
        for(int ij=-ghost_cells;ij<mn+ghost_cells;ij++) phi_1d_z(ij)=phi_ghost(i,j,ij);
        for(int ij=0;ij<mn;ij++){w_1d(ij)=V(i,j,ij).z;distance_1d_z(ij)=phi(i,j,ij);}
        if(spatial_order == 5) Local_WENO_Advect(mn,dz,phi_1d_z,w_1d,distance_1d_z,w_phiz_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,mn,dz,phi_1d_z,w_1d,distance_1d_z,w_phiz_1d,half_band_width);
        for(int ij=0;ij<mn;ij++) rhs(i,j,ij)+=w_phiz_1d(ij);}
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class TV> void LEVELSET_ADVECTION<TV>::
Euler_Step(const ARRAY<TV,TV_INT>& V,const T dt,const T time,const int number_of_ghost_cells)
{
    GRID<TV>& grid=levelset->grid;
    T_ARRAYS_SCALAR& phi=levelset->phi;
    T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(number_of_ghost_cells));levelset->boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells);

    if(local_semi_lagrangian_advection){
        LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(abs(phi_ghost(cell)) <= levelset->half_band_width)
            phi(cell)=interpolation.Clamped_To_Array(grid,phi_ghost,iterator.Location()-dt*V(cell));}}
    else if(local_advection_spatial_order){
        T_ARRAYS_SCALAR rhs(grid.Domain_Indices());
        Euler_Step_High_Order_Helper(grid,V,phi,phi_ghost,rhs,local_advection_spatial_order,levelset->half_band_width);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(abs(phi(cell)) <= levelset->half_band_width) phi(cell)-=dt*rhs(cell);}}
    else // use the advection routine in the level set base class
        advection->Update_Advection_Equation_Node(grid,phi,phi_ghost,V,*levelset->boundary,dt,time);

    levelset->boundary->Apply_Boundary_Condition(grid,phi,time+dt);
}
//#####################################################################
template class LEVELSET_ADVECTION<VECTOR<float,1> >;
template class LEVELSET_ADVECTION<VECTOR<float,2> >;
template class LEVELSET_ADVECTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_ADVECTION<VECTOR<double,1> >;
template class LEVELSET_ADVECTION<VECTOR<double,2> >;
template class LEVELSET_ADVECTION<VECTOR<double,3> >;
#endif
