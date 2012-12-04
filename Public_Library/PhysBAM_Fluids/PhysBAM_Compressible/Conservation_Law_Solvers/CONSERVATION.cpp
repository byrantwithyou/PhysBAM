//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION  
//##################################################################### 
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UTILITIES.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/BOUNDARY_OBJECT.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/BOUNDARY_OBJECT_REFLECTION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,int d> CONSERVATION<T_GRID,d>::
CONSERVATION()
    :save_fluxes(1),object_boundary_default(*new BOUNDARY_OBJECT_REFLECTION<T_GRID,TV_DIMENSION>),use_exact_neumann_face_location(false),
    scale_outgoing_fluxes_to_clamp_variable(false),clamped_variable_index(0),clamped_value(0),clamp_rho(.5),clamp_e(.5),min_dt(0),adaptive_time_step(false)
{
    Set_Order();
    Use_Field_By_Field_Alpha();
    Amplify_Alpha();
    save_fluxes=save_fluxes||scale_outgoing_fluxes_to_clamp_variable;
    object_boundary=&object_boundary_default;
}
//#####################################################################
// Destructor
//##################################################################### 
template<class T_GRID,int d> CONSERVATION<T_GRID,d>::
~CONSERVATION()
{
    delete &object_boundary_default;
}
//#####################################################################
// Function Alpha
//#####################################################################
template<class T_GRID,int d> typename T_GRID::SCALAR CONSERVATION<T_GRID,d>::
Alpha(const ARRAY<T,VECTOR<int,1> >& lambda_left,const ARRAY<T,VECTOR<int,1> >& lambda_right,const int k,const int length)
{
    if(field_by_field_alpha) return amplification_factor*maxabs(lambda_left(k),lambda_right(k));
    else{
        T lambda_max=0;for(int kk=0;kk<length;kk++) lambda_max=maxabs(lambda_max,lambda_left(kk),lambda_right(kk));
        return amplification_factor*lambda_max;}
}
//#####################################################################
// Function Compute_Delta_Flux_For_Clamping_Variable
//#####################################################################
template<class T_GRID,int d> void CONSERVATION<T_GRID,d>::
Compute_Delta_Flux_For_Clamping_Variable(const T_GRID& grid,const int number_of_ghost_cells,T dt,const int clamped_variable_index,const T clamped_value,const T_FACE_ARRAYS_BOOL& psi_N,
    const T_ARRAYS_DIMENSION_SCALAR& U,const T_FACE_ARRAYS_DIMENSION_SCALAR& flux,T_FACE_ARRAYS_DIMENSION_SCALAR& delta_flux,T_ARRAYS_DIMENSION_SCALAR& rhs,T_ARRAYS_SCALAR& overshoot_percentages)
{
    TV one_over_dx=grid.one_over_dX;
    VECTOR<bool,T_GRID::dimension*2> clamp_flux;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        T outgoing_flux=0;overshoot_percentages(cell_index)=0;
        for(int i=0;i<T_GRID::dimension*2;i++) clamp_flux(i)=false;
        for(int axis=0;axis<T_GRID::dimension;axis++){TV_INT first_face_index=iterator.First_Face_Index(axis),second_face_index=iterator.Second_Face_Index(axis);
            if(flux.Component(axis)(first_face_index)(clamped_variable_index)<0){clamp_flux(2*axis)=true;
                outgoing_flux-=dt*flux.Component(axis)(first_face_index)(clamped_variable_index)*one_over_dx[axis];}
            if(flux.Component(axis)(second_face_index)(clamped_variable_index)>0){clamp_flux(2*axis+1)=true;
                outgoing_flux+=dt*flux.Component(axis)(second_face_index)(clamped_variable_index)*one_over_dx[axis];}}
        
        if((U(cell_index)(clamped_variable_index)-clamped_value)<=0) overshoot_percentages(cell_index)=1;
        else{
            T clamped_variable_value=U(cell_index)(clamped_variable_index)-outgoing_flux;
            if((clamped_variable_value-clamped_value)>=0)
                overshoot_percentages(cell_index)=0;
            else{
                T overshoot=clamped_value-clamped_variable_value;
                overshoot_percentages(cell_index)=overshoot/outgoing_flux;}}
        if(overshoot_percentages(cell_index))
            for(int axis=0;axis<T_GRID::dimension;axis++){TV_INT first_face_index=iterator.First_Face_Index(axis),second_face_index=iterator.Second_Face_Index(axis);
                if(clamp_flux(2*axis)){
                    if(!psi_N.Component(axis)(first_face_index))
                        delta_flux.Component(axis)(first_face_index)(clamped_variable_index)=(-overshoot_percentages(cell_index))*flux.Component(axis)(first_face_index)(clamped_variable_index);
                    else rhs(cell_index)(clamped_variable_index)+=overshoot_percentages(cell_index)*flux.Component(axis)(first_face_index)(clamped_variable_index)*one_over_dx[axis];}
                
                if(clamp_flux(2*axis+1)){
                    if(!psi_N.Component(axis)(second_face_index))
                        delta_flux.Component(axis)(second_face_index)(clamped_variable_index)=(-overshoot_percentages(cell_index))*flux.Component(axis)(second_face_index)(clamped_variable_index);
                    else rhs(cell_index)(clamped_variable_index)-=overshoot_percentages(cell_index)*flux.Component(axis)(second_face_index)(clamped_variable_index)*one_over_dx[axis];}}}
}
//#####################################################################
// Function Compute_Flux_Without_Clamping
//#####################################################################
template<class T_GRID,int d> void CONSERVATION<T_GRID,d>::
Compute_Flux_Without_Clamping(const T_GRID& grid,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
    VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems,VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems_explicit,const T_FACE_ARRAYS_BOOL& psi_N,
    const T_FACE_ARRAYS_SCALAR& face_velocities,const TV_BOOL& outflow_boundaries,T_ARRAYS_DIMENSION_SCALAR& rhs,const bool thinshell,const T_ARRAYS_DIMENSION_SCALAR* U_ghost_clamped)
{
    RANGE<TV_INT> U_domain_indices=U.Domain_Indices();
    RANGE<TV_INT> U_ghost_domain_indices=U_ghost.Domain_Indices();
    int U_start,U_end,U_ghost_start,U_ghost_end;
    TV dx=grid.dX;
    FLOOD_FILL<1> find_connected_components;

    for(int axis=0;axis<T_GRID::dimension;axis++){
        U_start=U_domain_indices.min_corner(axis);U_end=U_domain_indices.max_corner(axis);
        U_ghost_start=U_ghost_domain_indices.min_corner(axis);U_ghost_end=U_ghost_domain_indices.max_corner(axis);
        ARRAY<TV_DIMENSION,VECTOR<int,1> > U_1d_axis(U_ghost_start,U_ghost_end),flux_axis_1d(U_start,U_end);
        if(U_ghost_clamped) U_flux_1d_axis.Resize(U_ghost_start,U_ghost_end);
        ARRAY<bool,VECTOR<int,1> > psi_axis(U_start,U_end),psi_N_axis(U_start,U_end+1);
        VECTOR<bool,2> outflow_boundaries_axis(outflow_boundaries(2*axis),outflow_boundaries(2*axis+1));
        ARRAY<int,VECTOR<int,1> > filled_region_colors(U_start,U_end);filled_region_colors.Fill(-1);
        ARRAY<bool,VECTOR<int,1> > psi_axis_current_component(U_start,U_end);
        if(save_fluxes) flux_temp.Resize(U_start-1,U_end,true,false);
        GRID<TV_LOWER_DIM> lower_dimension_grid=grid.Remove_Dimension(axis);
        for(UNIFORM_GRID_ITERATOR_CELL<TV_LOWER_DIM> iterator(lower_dimension_grid);iterator.Valid();iterator.Next()){VECTOR<int,TV::m-1> cell_index=iterator.Cell_Index();
            VECTOR<int,3> slice_index;TV_INT cell_index_full_dimension=cell_index.Insert(0,axis);
            for(int axis_slice=0;axis_slice<T_GRID::dimension;axis_slice++){
                slice_index[axis_slice]=cell_index_full_dimension[axis_slice];}

            for(int i=U_start;i<U_end;i++){
                psi_axis(i)=psi(cell_index.Insert(i,axis));
                filled_region_colors(i)=psi_axis(i)?-1:-2;}
            for(int i=U_start;i<U_end+1;i++) psi_N_axis(i)=psi_N(axis,cell_index.Insert(i,axis));
            int number_of_regions=find_connected_components.Flood_Fill(filled_region_colors,psi_N_axis);
            for(int color=0;color<number_of_regions;color++){
                for(int i=U_ghost_start;i<U_ghost_end;i++) U_1d_axis(i)=U_ghost(cell_index.Insert(i,axis));
                if(U_ghost_clamped) for(int i=U_ghost_start;i<U_ghost_end;i++) U_flux_1d_axis(i)=(*U_ghost_clamped)(cell_index.Insert(i,axis));
                psi_axis_current_component.Fill(false);
                for(int i=U_start;i<U_end;i++) psi_axis_current_component(i)=(filled_region_colors(i)==color);
                VECTOR<int,2> region_boundary=find_connected_components.region_boundaries(color);
                VECTOR<bool,2> psi_N_boundary(psi_N_axis(region_boundary.x),psi_N_axis(region_boundary.y+1));
                if(thinshell) object_boundary->Fill_Ghost_Cells_Neumann(grid.Get_1D_Grid(axis),U_1d_axis,face_velocities,cell_index,axis,order,use_exact_neumann_face_location,
                    VECTOR<int,2>(U_start,U_end),find_connected_components.region_boundaries(color),psi_N_boundary,callbacks);
                if(U_ghost_clamped)
                    if(thinshell) object_boundary->Fill_Ghost_Cells_Neumann(grid.Get_1D_Grid(axis),U_flux_1d_axis,face_velocities,cell_index,axis,order,use_exact_neumann_face_location,
                        VECTOR<int,2>(U_start,U_end),find_connected_components.region_boundaries(color),psi_N_boundary,callbacks);
                VECTOR<bool,2> outflow_boundaries_current_component;
                outflow_boundaries_current_component(0)=outflow_boundaries_axis(0)&&(!psi_N_boundary(0));outflow_boundaries_current_component(0)=outflow_boundaries_axis(0)&&(!psi_N_boundary(0));
                (eigensystems[axis])->slice_index=slice_index;(eigensystems_explicit[axis])->slice_index=slice_index;
                ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_flux_pointer=U_ghost_clamped?(&U_flux_1d_axis):0;
                Conservation_Solver(U_end,dx[axis],psi_axis_current_component,U_1d_axis,flux_axis_1d,*eigensystems[axis],*eigensystems_explicit[axis],
                    outflow_boundaries_current_component,U_flux_pointer);}
            for(int i=U_start;i<U_end;i++) for(int k=0;k<d;k++)
                if(!scale_outgoing_fluxes_to_clamp_variable||(!U_ghost_clamped&&k==clamped_variable_index)||(U_ghost_clamped&&k!=clamped_variable_index))
                    rhs(cell_index.Insert(i,axis))(k)+=flux_axis_1d(i)(k);
            if(save_fluxes) for(int i=U_start-1;i<U_end;i++) for(int k=0;k<d;k++)
                if(!scale_outgoing_fluxes_to_clamp_variable||(!U_ghost_clamped&&k==clamped_variable_index)||(U_ghost_clamped&&k!=clamped_variable_index))
                    fluxes.Component(axis)(cell_index.Insert(i+1,axis))(k)=flux_temp(i)(k);}}
}
template<class T_GRID,int d> void CONSERVATION<T_GRID,d>::
Compute_Flux_With_Clamping(const T_GRID& grid,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
    VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems,VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems_explicit,const T_FACE_ARRAYS_BOOL& psi_N,
    const T_FACE_ARRAYS_SCALAR& face_velocities,const TV_BOOL& outflow_boundaries,T_ARRAYS_DIMENSION_SCALAR& rhs,const bool thinshell,
    VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>* eigensystems_auxiliary,T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary)
{
    RANGE<TV_INT> U_domain_indices=U.Domain_Indices();
    Compute_Flux_Without_Clamping(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs,thinshell); // After this call, we should have rhs for clamped variable
    // Use this rhs to update the fluxes for the clamped variable
    T_FACE_ARRAYS_DIMENSION_SCALAR delta_flux(grid);
    T_ARRAYS_DIMENSION_SCALAR delta_rhs(U_domain_indices,true);
    T_ARRAYS_SCALAR overshoot_percentages(grid.Domain_Indices());

    Compute_Delta_Flux_For_Clamping_Variable(grid,U_domain_indices.max_corner.x-grid.Domain_Indices().max_corner.x,dt,clamped_variable_index,clamped_value,psi_N,U,fluxes,delta_flux,rhs,
        overshoot_percentages);
    ARRAYS_UTILITIES<T_GRID,TV_DIMENSION>::Compute_Divergence_At_Cells_From_Face_Data(grid,delta_rhs,delta_flux,0);
    rhs+=delta_rhs;
    T_ARRAYS_DIMENSION_SCALAR U_ghost_clamped(U_ghost,true);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,U_domain_indices.max_corner.x-grid.Domain_Indices().max_corner.x);iterator.Valid();iterator.Next())
        U_ghost_clamped(iterator.Cell_Index())=(1-overshoot_percentages(iterator.Cell_Index()))*U_ghost(iterator.Cell_Index());
    // Compute the flux for the other variables
    if(fluxes_auxiliary){
        T_ARRAYS_DIMENSION_SCALAR rhs_auxiliary(U_domain_indices,true);
        Compute_Flux_Without_Clamping(grid,U,U_ghost,psi,dt,*eigensystems_auxiliary,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs_auxiliary,thinshell,&U_ghost_clamped);
        *fluxes_auxiliary=fluxes;}
    Compute_Flux_Without_Clamping(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs,thinshell,&U_ghost_clamped);
}
template<class T_GRID,int d> void CONSERVATION<T_GRID,d>::
Compute_Flux(const T_GRID& grid,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
    VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems,VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems_explicit,const T_FACE_ARRAYS_BOOL& psi_N,
    const T_FACE_ARRAYS_SCALAR& face_velocities,const TV_BOOL& outflow_boundaries,T_ARRAYS_DIMENSION_SCALAR& rhs,const bool thinshell,
    VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>* eigensystems_auxiliary,T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary)
{
    if(scale_outgoing_fluxes_to_clamp_variable) Compute_Flux_With_Clamping(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs,thinshell,
        eigensystems_auxiliary,fluxes_auxiliary);
    else{
        if(fluxes_auxiliary){
            RANGE<TV_INT> U_domain_indices=U.Domain_Indices();T_ARRAYS_DIMENSION_SCALAR rhs_auxiliary(U_domain_indices,true);
            Compute_Flux_Without_Clamping(grid,U,U_ghost,psi,dt,*eigensystems_auxiliary,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs_auxiliary,thinshell);
            *fluxes_auxiliary=fluxes;}
        Compute_Flux_Without_Clamping(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs,thinshell);}
}
//#####################################################################
// Function Update_Conservation_Law
//#####################################################################
template<class T_GRID,int d> void CONSERVATION<T_GRID,d>::
Update_Conservation_Law(T_GRID& grid,T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
    VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems,VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems_explicit,const T_FACE_ARRAYS_BOOL& psi_N,
    const T_FACE_ARRAYS_SCALAR& face_velocities,const bool thinshell,const TV_BOOL& outflow_boundaries,VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>* eigensystems_auxiliary,
    T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary)
{
    RANGE<TV_INT> U_domain_indices=U.Domain_Indices();
    T_ARRAYS_DIMENSION_SCALAR rhs(U_domain_indices,true);

    if(fluxes_auxiliary) save_fluxes=true;
    if(save_fluxes) fluxes.Resize(grid.Domain_Indices(),true,false);

    Compute_Flux(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs,thinshell,eigensystems_auxiliary,fluxes_auxiliary);

    T_ARRAYS_SCALAR rho_dt(grid.Domain_Indices()), e_dt(grid.Domain_Indices());
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        T clamp_rho_cell=clamp_rho*U_ghost(cell_index)(0);
        if(rhs(cell_index)(0)>0 && U_ghost(cell_index)(0)>=clamp_rho_cell) rho_dt(cell_index)=(U_ghost(cell_index)(0)-clamp_rho_cell)/rhs(cell_index)(0);
        else rho_dt(cell_index)=dt;

        T momentum_dt=dt;
        T clamp_e_cell=clamp_e*EULER<T_GRID>::e(U,cell_index);
        for(int axis=0;axis<TV::dimension;axis++){
            T tmp_dt=dt;
            T momentum_flux_sqr=rhs(cell_index)(axis+1)*rhs(cell_index)(axis+1);
            T a=2*rhs(cell_index)(0)*rhs(cell_index)(d-1)-momentum_flux_sqr-2*clamp_e_cell*rhs(cell_index)(0)*rhs(cell_index)(0);
            T c=2*U_ghost(cell_index)(d-1)*U_ghost(cell_index)(0)-2*clamp_e_cell*U_ghost(cell_index)(0)*U_ghost(cell_index)(0)-U_ghost(cell_index)(axis+1)*U_ghost(cell_index)(axis+1);
            if(dt*a>0){
                T b_over_two=U_ghost(cell_index)(axis+1)*rhs(cell_index)(axis+1)-U_ghost(cell_index)(0)*rhs(cell_index)(d-1)-U_ghost(cell_index)(d-1)*rhs(cell_index)(0)+2*clamp_e_cell*U_ghost(cell_index)(0)*rhs(cell_index)(0);
                T b_sqr_over_four=b_over_two*b_over_two;
                T ac=a*c;
                if(b_sqr_over_four>ac){
                    T tmp_dt1=(-1*b_over_two-sqrt(b_sqr_over_four-ac))/a;
                    T tmp_dt2=(-1*b_over_two+sqrt(b_sqr_over_four-ac))/a;
                    tmp_dt=(tmp_dt1>=0)?(tmp_dt2>=0)?min(tmp_dt1,tmp_dt2):tmp_dt1:(tmp_dt2>=0)?tmp_dt2:dt;}
                else tmp_dt=dt;}
            momentum_dt=min(momentum_dt,tmp_dt);}
        e_dt(cell_index)=dt;
        for(int axis=0;axis<TV::dimension;axis++) e_dt(cell_index)=min(e_dt(cell_index),momentum_dt);}

    min_dt=min(rho_dt.Min(),e_dt.Min());
    LOG::cout<<"dt: "<<dt<<" Min dt: "<<min_dt<<std::endl;

    if((min_dt==dt)||(!adaptive_time_step)){
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,U_domain_indices.max_corner.x-grid.Domain_Indices().max_corner.x);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();
            if(psi(cell_index)) U(cell_index)-=dt*rhs(cell_index);}}
}
//#####################################################################
// Function Update_Conservation_Law_For_Specialized_Shallow_Water_Equations
//#####################################################################
// A specialized function only used in SHALLOW_WATER_2D_SPECIALIZED
template<class T_GRID,int d> template<class T_ARRAYS> void CONSERVATION<T_GRID,d>::
Update_Conservation_Law_For_Specialized_Shallow_Water_Equations(GRID<TV>& grid,T_ARRAYS& U,const T_ARRAYS& U_ghost,const ARRAY<bool,VECTOR<int,2> >& psi,const T dt,
    EIGENSYSTEM<T,VECTOR<T,2> >& eigensystem_F,EIGENSYSTEM<T,VECTOR<T,2> >& eigensystem_G,CONSERVATION<GRID<TV>,2>& solver,const TV_BOOL& outflow_boundaries)
{
    STATIC_ASSERT((IS_SAME<T_ARRAYS,ARRAY<VECTOR<T,3> ,VECTOR<int,2> > >::value));
    if(save_fluxes!=solver.save_fluxes) PHYSBAM_FATAL_ERROR();

    int i,j;int m=grid.counts.x,n=grid.counts.y;T dx=grid.dX.x,dy=grid.dX.y;
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > rhs(0,m,0,n);

    if(save_fluxes) fluxes.Resize(grid);

    if(save_fluxes) solver.flux_temp.Resize(0,m,true,false);
    ARRAY<VECTOR<T,2> ,VECTOR<int,1> > U_1d_x(-3,m+3),Fx_1d(0,m);
    ARRAY<bool,VECTOR<int,1> > psi_x(0,m);
    for(j=0;j<n;j++){
        for(i=0;i<m;i++) psi_x(i)=psi(i,j);
        for(i=-3;i<m+3;i++){U_1d_x(i)(0)=U_ghost(i,j)(0);U_1d_x(i)(1)=U_ghost(i,j)(1);}
        eigensystem_F.slice_index=VECTOR<int,3>(0,j,0);
        solver.Conservation_Solver(m,dx,psi_x,U_1d_x,Fx_1d,eigensystem_F,eigensystem_F,VECTOR<bool,2>(outflow_boundaries(0),outflow_boundaries(1)));
        for(i=0;i<m;i++){rhs(i,j)(0)=Fx_1d(i)(0);rhs(i,j)(1)=Fx_1d(i)(1);}
        if(save_fluxes) 
            for(i=0;i<m;i++){fluxes.Component(0)(i+1,j)(0)=solver.flux_temp(i)(0);fluxes.Component(0)(i+1,j)(1)=solver.flux_temp(i)(1);fluxes.Component(0)(i+1,j)(2)=0;}}

    if(save_fluxes) solver.flux_temp.Resize(0,n,true,false);
    ARRAY<VECTOR<T,2> ,VECTOR<int,1> > U_1d_y(-3,n+3),Gy_1d(0,n);
    ARRAY<bool,VECTOR<int,1> > psi_y(0,n);
    for(i=0;i<m;i++){
        for(j=0;j<n;j++) psi_y(j)=psi(i,j);
        for(j=-3;j<n+3;j++){U_1d_y(j)(0)=U_ghost(i,j)(0);U_1d_y(j)(1)=U_ghost(i,j)(2);}
        eigensystem_G.slice_index=VECTOR<int,3>(i,0,0);
        solver.Conservation_Solver(n,dy,psi_y,U_1d_y,Gy_1d,eigensystem_G,eigensystem_G,VECTOR<bool,2>(outflow_boundaries(2),outflow_boundaries(3)));
        for(j=0;j<n;j++){rhs(i,j)(0)+=Gy_1d(j)(0);rhs(i,j)(2)=Gy_1d(j)(1);}
        if(save_fluxes)
            for(j=0;j<n;j++){fluxes.Component(1)(i,j+1)(0)=solver.flux_temp(j)(0);fluxes.Component(1)(i,j+1)(1)=0;fluxes.Component(1)(i,j+1)(2)=solver.flux_temp(j)(1);}}

    for(i=0;i<m;i++) for(j=0;j<n;j++) if(psi(i,j)) U(i,j)-=dt*rhs(i,j);
}
//#####################################################################
// Function Log_Parameters
//#####################################################################
template<class T_GRID,int d> void CONSERVATION<T_GRID,d>::
Log_Parameters() const
{
    LOG::SCOPE scope("CONSERVATION parameters");
    LOG::cout<<"order="<<order<<std::endl;
    LOG::cout<<"field_by_field_alpha="<<field_by_field_alpha<<std::endl;
    LOG::cout<<"amplification_factor="<<amplification_factor<<std::endl;
    LOG::cout<<"save_fluxes="<<save_fluxes<<std::endl;
    LOG::cout<<"use_exact_neumann_face_location="<<use_exact_neumann_face_location<<std::endl;
    LOG::cout<<"scale_outgoing_fluxes_to_clamp_variable="<<scale_outgoing_fluxes_to_clamp_variable<<std::endl;
    LOG::cout<<"clamped_variable_index="<<clamped_variable_index<<std::endl;
    LOG::cout<<"clamped_value="<<clamped_value<<std::endl;
}
//#####################################################################
namespace PhysBAM{
template class CONSERVATION<GRID<VECTOR<double,1> >,3>;
template class CONSERVATION<GRID<VECTOR<double,2> >,4>;
template class CONSERVATION<GRID<VECTOR<double,3> >,5>;
template class CONSERVATION<GRID<VECTOR<float,1> >,3>;
template class CONSERVATION<GRID<VECTOR<float,2> >,4>;
template class CONSERVATION<GRID<VECTOR<float,3> >,5>;
}
