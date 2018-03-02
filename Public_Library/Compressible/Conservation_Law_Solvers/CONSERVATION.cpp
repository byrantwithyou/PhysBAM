//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION  
//##################################################################### 
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Grid_Tools/Arrays/ARRAYS_UTILITIES.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Compressible/Conservation_Law_Solvers/BOUNDARY_OBJECT.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION.h>
#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
#include <Compressible/Euler_Equations/BOUNDARY_OBJECT_REFLECTION.h>
#include <Compressible/Euler_Equations/EULER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> CONSERVATION<TV,d>::
CONSERVATION()
    :save_fluxes(1),object_boundary_default(*new BOUNDARY_OBJECT_REFLECTION<TV,TV_DIMENSION>),use_exact_neumann_face_location(false),
    scale_outgoing_fluxes_to_clamp_variable(false),clamped_variable_index(0),clamped_value(0)
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
template<class TV,int d> CONSERVATION<TV,d>::
~CONSERVATION()
{
    delete &object_boundary_default;
}
//#####################################################################
// Function Alpha
//#####################################################################
template<class TV,int d> typename TV::SCALAR CONSERVATION<TV,d>::
Alpha(const VECTOR<T,d>& lambda_left,const VECTOR<T,d>& lambda_right,const int k,const int length)
{
    if(field_by_field_alpha) return amplification_factor*maxabs(lambda_left(k),lambda_right(k));
    else{
        T lambda_max=0;for(int kk=0;kk<length;kk++) lambda_max=maxabs(lambda_max,lambda_left(kk),lambda_right(kk));
        return amplification_factor*lambda_max;}
}
//#####################################################################
// Function Compute_Delta_Flux_For_Clamping_Variable
//#####################################################################
template<class TV,int d> void CONSERVATION<TV,d>::
Compute_Delta_Flux_For_Clamping_Variable(const GRID<TV>& grid,const int number_of_ghost_cells,T dt,const int clamped_variable_index,const T clamped_value,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
    const T_ARRAYS_DIMENSION_SCALAR& U,const T_FACE_ARRAYS_DIMENSION_SCALAR& flux,T_FACE_ARRAYS_DIMENSION_SCALAR& delta_flux,T_ARRAYS_DIMENSION_SCALAR& rhs,ARRAY<T,TV_INT>& overshoot_percentages)
{
    TV one_over_dx=grid.one_over_dX;
    VECTOR<bool,TV::m*2> clamp_flux;
    for(CELL_ITERATOR<TV> iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        T outgoing_flux=0;overshoot_percentages(cell_index)=0;
        for(int i=0;i<TV::m*2;i++) clamp_flux(i)=false;
        for(int axis=0;axis<TV::m;axis++){TV_INT first_face_index=iterator.First_Face_Index(axis),second_face_index=iterator.Second_Face_Index(axis);
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
            for(int axis=0;axis<TV::m;axis++){TV_INT first_face_index=iterator.First_Face_Index(axis),second_face_index=iterator.Second_Face_Index(axis);
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
template<class TV,int d> void CONSERVATION<TV,d>::
Compute_Flux_Without_Clamping(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
    VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems,VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems_explicit,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
    const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const VECTOR<bool,2*TV::m>& outflow_boundaries,T_ARRAYS_DIMENSION_SCALAR& rhs,const bool thinshell,const T_ARRAYS_DIMENSION_SCALAR* U_ghost_clamped)
{
    RANGE<TV_INT> U_domain_indices=U.Domain_Indices();
    RANGE<TV_INT> U_ghost_domain_indices=U_ghost.Domain_Indices();
    int U_start,U_end,U_ghost_start,U_ghost_end;
    TV dx=grid.dX;
    FLOOD_FILL<1> find_connected_components;

    for(int axis=0;axis<TV::m;axis++){
        RANGE<VECTOR<int,1> > U_range=U_domain_indices.Dimension_Range(axis);
        RANGE<VECTOR<int,1> > U_range_ext=U_range;
        RANGE<VECTOR<int,1> > U_range_less=U_range;
        U_range_ext.max_corner.x++;
        U_range_less.min_corner.x--;
        RANGE<VECTOR<int,1> > U_ghost_range=U_ghost_domain_indices.Dimension_Range(axis);
        U_start=U_domain_indices.min_corner(axis);U_end=U_domain_indices.max_corner(axis);
        U_ghost_start=U_ghost_domain_indices.min_corner(axis);U_ghost_end=U_ghost_domain_indices.max_corner(axis);
        ARRAY<TV_DIMENSION,VECTOR<int,1> > U_1d_axis(U_ghost_range),flux_axis_1d(U_range);
        if(U_ghost_clamped) U_flux_1d_axis.Resize(U_ghost_range);
        ARRAY<bool,VECTOR<int,1> > psi_axis(U_range),psi_N_axis(U_range_ext);
        VECTOR<bool,2> outflow_boundaries_axis(outflow_boundaries(2*axis),outflow_boundaries(2*axis+1));
        ARRAY<int,VECTOR<int,1> > filled_region_colors(U_range);filled_region_colors.Fill(-1);
        ARRAY<bool,VECTOR<int,1> > psi_axis_current_component(U_range);
        if(save_fluxes) flux_temp.Resize(U_range_less);
        GRID<TV_LOWER_DIM> lower_dimension_grid=grid.Remove_Dimension(axis);
        for(CELL_ITERATOR<TV_LOWER_DIM> iterator(lower_dimension_grid);iterator.Valid();iterator.Next()){VECTOR<int,TV::m-1> cell_index=iterator.Cell_Index();
            VECTOR<int,3> slice_index;TV_INT cell_index_full_dimension=cell_index.Insert(0,axis);
            for(int axis_slice=0;axis_slice<TV::m;axis_slice++){
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
                if(thinshell) object_boundary->Fill_Ghost_Cells_Neumann(grid.Get_1D_Grid(axis),U_1d_axis,face_velocities,cell_index,axis,order,use_exact_neumann_face_location,VECTOR<int,2>(U_range.min_corner.x,U_range.max_corner.x),find_connected_components.region_boundaries(color),psi_N_boundary,callbacks);
                if(U_ghost_clamped)
                    if(thinshell) object_boundary->Fill_Ghost_Cells_Neumann(grid.Get_1D_Grid(axis),U_flux_1d_axis,face_velocities,cell_index,axis,order,use_exact_neumann_face_location,
                        VECTOR<int,2>(U_range.min_corner.x,U_range.max_corner.x),find_connected_components.region_boundaries(color),psi_N_boundary,callbacks);
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
template<class TV,int d> void CONSERVATION<TV,d>::
Compute_Flux_With_Clamping(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
    VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems,VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems_explicit,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
    const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const VECTOR<bool,2*TV::m>& outflow_boundaries,T_ARRAYS_DIMENSION_SCALAR& rhs,const bool thinshell,
    VECTOR<EIGENSYSTEM<T,d>*,TV::m>* eigensystems_auxiliary,T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary)
{
    RANGE<TV_INT> U_domain_indices=U.Domain_Indices();
    Compute_Flux_Without_Clamping(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs,thinshell); // After this call, we should have rhs for clamped variable
    // Use this rhs to update the fluxes for the clamped variable
    T_FACE_ARRAYS_DIMENSION_SCALAR delta_flux(grid);
    T_ARRAYS_DIMENSION_SCALAR delta_rhs(U_domain_indices);
    ARRAY<T,TV_INT> overshoot_percentages(grid.Domain_Indices());

    Compute_Delta_Flux_For_Clamping_Variable(grid,U_domain_indices.max_corner.x-grid.Domain_Indices().max_corner.x,dt,clamped_variable_index,clamped_value,psi_N,U,fluxes,delta_flux,rhs,
        overshoot_percentages);
    ARRAYS_UTILITIES<TV,TV_DIMENSION>::Compute_Divergence_At_Cells_From_Face_Data(grid,delta_rhs,delta_flux,0);
    rhs+=delta_rhs;
    T_ARRAYS_DIMENSION_SCALAR U_ghost_clamped(U_ghost);
    for(CELL_ITERATOR<TV> iterator(grid,U_domain_indices.max_corner.x-grid.Domain_Indices().max_corner.x);iterator.Valid();iterator.Next())
        U_ghost_clamped(iterator.Cell_Index())=(1-overshoot_percentages(iterator.Cell_Index()))*U_ghost(iterator.Cell_Index());
    // Compute the flux for the other variables
    if(fluxes_auxiliary){
        T_ARRAYS_DIMENSION_SCALAR rhs_auxiliary(U_domain_indices);
        Compute_Flux_Without_Clamping(grid,U,U_ghost,psi,dt,*eigensystems_auxiliary,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs_auxiliary,thinshell,&U_ghost_clamped);
        *fluxes_auxiliary=fluxes;}
    Compute_Flux_Without_Clamping(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs,thinshell,&U_ghost_clamped);
}
template<class TV,int d> void CONSERVATION<TV,d>::
Compute_Flux(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
    VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems,VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems_explicit,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
    const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const VECTOR<bool,2*TV::m>& outflow_boundaries,T_ARRAYS_DIMENSION_SCALAR& rhs,const bool thinshell,
    VECTOR<EIGENSYSTEM<T,d>*,TV::m>* eigensystems_auxiliary,T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary)
{
    if(scale_outgoing_fluxes_to_clamp_variable) Compute_Flux_With_Clamping(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs,thinshell,
        eigensystems_auxiliary,fluxes_auxiliary);
    else{
        if(fluxes_auxiliary){
            RANGE<TV_INT> U_domain_indices=U.Domain_Indices();
            T_ARRAYS_DIMENSION_SCALAR rhs_auxiliary(U_domain_indices);
            Compute_Flux_Without_Clamping(grid,U,U_ghost,psi,dt,*eigensystems_auxiliary,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs_auxiliary,thinshell);
            *fluxes_auxiliary=fluxes;}
        Compute_Flux_Without_Clamping(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs,thinshell);}
}
//#####################################################################
// Function Update_Conservation_Law
//#####################################################################
template<class TV,int d> void CONSERVATION<TV,d>::
Update_Conservation_Law(GRID<TV>& grid,T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
    VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems,VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems_explicit,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
    const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const bool thinshell,const VECTOR<bool,2*TV::m>& outflow_boundaries,VECTOR<EIGENSYSTEM<T,d>*,TV::m>* eigensystems_auxiliary,
    T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary)
{
    RANGE<TV_INT> U_domain_indices=U.Domain_Indices();
    T_ARRAYS_DIMENSION_SCALAR rhs(U_domain_indices);

    if(fluxes_auxiliary) save_fluxes=true;
    if(save_fluxes) fluxes.Resize(grid.Domain_Indices(),true,false);

    Compute_Flux(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs,thinshell,eigensystems_auxiliary,fluxes_auxiliary);

    for(CELL_ITERATOR<TV> iterator(grid,U_domain_indices.max_corner.x-grid.Domain_Indices().max_corner.x);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(psi(cell_index)) U(cell_index)-=dt*rhs(cell_index);}
}
//#####################################################################
// Function Log_Parameters
//#####################################################################
template<class TV,int d> void CONSERVATION<TV,d>::
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
template class CONSERVATION<VECTOR<double,1>,1>;
template class CONSERVATION<VECTOR<double,1>,2>;
template class CONSERVATION<VECTOR<double,1>,3>;
template class CONSERVATION<VECTOR<double,1>,4>;
template class CONSERVATION<VECTOR<double,2>,2>;
template class CONSERVATION<VECTOR<double,2>,3>;
template class CONSERVATION<VECTOR<double,2>,4>;
template class CONSERVATION<VECTOR<double,2>,5>;
template class CONSERVATION<VECTOR<double,3>,5>;
template class CONSERVATION<VECTOR<double,3>,6>;
template class CONSERVATION<VECTOR<float,1>,1>;
template class CONSERVATION<VECTOR<float,1>,2>;
template class CONSERVATION<VECTOR<float,1>,3>;
template class CONSERVATION<VECTOR<float,1>,4>;
template class CONSERVATION<VECTOR<float,2>,2>;
template class CONSERVATION<VECTOR<float,2>,3>;
template class CONSERVATION<VECTOR<float,2>,4>;
template class CONSERVATION<VECTOR<float,2>,5>;
template class CONSERVATION<VECTOR<float,3>,5>;
template class CONSERVATION<VECTOR<float,3>,6>;
}
