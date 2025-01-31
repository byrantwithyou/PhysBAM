//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYBRID_SL_ENO_CONSERVATION  
//##################################################################### 
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/SCOPE.h>
#include <Grid_Tools/Arrays/ARRAYS_UTILITIES.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Compressible/Conservation_Law_Solvers/BOUNDARY_OBJECT.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION.h>
#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
#include <Compressible/Conservation_Law_Solvers/HYBRID_SL_ENO_CONSERVATION.h>
#include <Compressible/Euler_Equations/BOUNDARY_OBJECT_REFLECTION.h>
#include <Compressible/Euler_Equations/EULER.h>
using namespace PhysBAM;
//#####################################################################
// Function Update_Conservation_Law
//#####################################################################
template<class TV,int d> void HYBRID_SL_ENO_CONSERVATION<TV,d>::
Update_Conservation_Law(GRID<TV>& grid,T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const ARRAY<bool,TV_INT>& psi,const T dt,
    VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems,VECTOR<EIGENSYSTEM<T,d>*,TV::m>& eigensystems_explicit,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,
    const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const bool thinshell,const VECTOR<bool,2*TV::m>& outflow_boundaries,VECTOR<EIGENSYSTEM<T,d>*,TV::m>* eigensystems_auxiliary,
    T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary)
{
    const T cell_volume=grid.Cell_Size();
    const T one_over_cell_volume=(T)1/cell_volume;

    ARRAY<bool,TV_INT> regular_cell(psi);
    ARRAY<bool,TV_INT> cell_near_interface(psi),cell_near_interface_tmp(psi);
    T_FACE_ARRAYS_DIMENSION_SCALAR& face_fluxes(conservation->fluxes);
    T_ARRAYS_DIMENSION_SCALAR rhs(U.Domain_Indices());

    {
        LOG::SCOPE scope("Regular Update for Hybrid scheme.");
        for(CELL_ITERATOR<TV> iterator(grid,3);iterator.Valid();iterator.Next()){
            if(regular_cell(iterator.Cell_Index())){
                bool compute_self_weight=true;for(int dim=0;dim<TV::m;dim++) compute_self_weight &= flux_face(iterator.Full_First_Face_Index(dim)) & flux_face(iterator.Full_Second_Face_Index(dim));
                regular_cell(iterator.Cell_Index())=compute_self_weight;
                cell_near_interface_tmp(iterator.Cell_Index())=!compute_self_weight;}}
        for(CELL_ITERATOR<TV> iterator(grid,2);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            cell_near_interface(cell_index)=cell_near_interface_tmp(cell_index);
            if(cell_near_interface_tmp(cell_index)){
                for(int dim=0;dim<TV::m;dim++) {
                    if(psi(cell_index+TV_INT::Axis_Vector(dim))) cell_near_interface(cell_index+TV_INT::Axis_Vector(dim))=true;
                    if(psi(cell_index-TV_INT::Axis_Vector(dim))) cell_near_interface(cell_index-TV_INT::Axis_Vector(dim))=true;}
#if 0
                if(TV::m==2){
                    if(psi(cell_index+VECTOR<int,2>(-1,-1))) cell_near_interface(cell_index+VECTOR<int,2>(-1,-1))=true; if(psi(cell_index+VECTOR<int,2>( 1,-1))) cell_near_interface(cell_index+VECTOR<int,2>( 1,-1))=true;
                    if(psi(cell_index+VECTOR<int,2>(-1, 1))) cell_near_interface(cell_index+VECTOR<int,2>(-1, 1))=true; if(psi(cell_index+VECTOR<int,2>( 1, 1))) cell_near_interface(cell_index+VECTOR<int,2>( 1, 1))=true;}
                if(TV::m==3){
                    if(psi(cell_index+VECTOR<int,3>(-1,-1,-1))) cell_near_interface(cell_index+VECTOR<int,3>(-1,-1,-1))=true; if(psi(cell_index+VECTOR<int,3>( 1,-1,-1))) cell_near_interface(cell_index+VECTOR<int,3>( 1,-1,-1))=true;
                    if(psi(cell_index+VECTOR<int,3>(-1, 1,-1))) cell_near_interface(cell_index+VECTOR<int,3>(-1, 1,-1))=true; if(psi(cell_index+VECTOR<int,3>( 1, 1,-1))) cell_near_interface(cell_index+VECTOR<int,3>( 1, 1,-1))=true;
                    if(psi(cell_index+VECTOR<int,3>(-1,-1, 1))) cell_near_interface(cell_index+VECTOR<int,3>(-1,-1, 1))=true; if(psi(cell_index+VECTOR<int,3>( 1,-1, 1))) cell_near_interface(cell_index+VECTOR<int,3>( 1,-1, 1))=true;
                    if(psi(cell_index+VECTOR<int,3>(-1, 1, 1))) cell_near_interface(cell_index+VECTOR<int,3>(-1, 1, 1))=true; if(psi(cell_index+VECTOR<int,3>( 1, 1, 1))) cell_near_interface(cell_index+VECTOR<int,3>( 1, 1, 1))=true;}
#endif
            }
        }

        conservation->fluxes.Resize(grid.Domain_Indices(),init_all);
        conservation->Compute_Flux(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs,thinshell,eigensystems_auxiliary,fluxes_auxiliary);
    }

    {
        LOG::SCOPE scope("Irregular update for hybrid scheme (no collision bodies present).");
        U.Fill(TV_DIMENSION());
        for(int dim=0;dim<TV_DIMENSION::dimension;dim++){
            ARRAY<PAIR<T,INDEX> > weights;
            ARRAY<ARRAY<int>,INDEX> donors;donors.Resize(grid.Domain_Indices(3));
            ARRAY<ARRAY<int>,INDEX> receivers;receivers.Resize(grid.Domain_Indices(3));
            ARRAY<T,TV_INT> sigma;sigma.Resize(grid.Domain_Indices(3));

            for(FACE_ITERATOR<TV> iterator(grid,2);iterator.Valid();iterator.Next()){
                FACE_INDEX<TV::m> face_index=iterator.Full_Index();
                if(flux_face(face_index) && (cell_near_interface(iterator.First_Cell_Index()) || cell_near_interface(iterator.Second_Cell_Index()))){
                    T weight=dt*face_fluxes(face_index)(dim);
                    INDEX donor_cell    = (weight >= 0 ? iterator.First_Cell_Index() : iterator.Second_Cell_Index());
                    INDEX receiver_cell = (weight >= 0 ? iterator.Second_Cell_Index() : iterator.First_Cell_Index());
                    weight = abs(weight);
                    int index=weights.Append(PAIR<T, INDEX>(weight,receiver_cell));
                    sigma(donor_cell)+=weight; donors(donor_cell).Append(index);receivers(receiver_cell).Append(index);}}

            // "Finish" cells that are entirely updated via the high order method.
            for(CELL_ITERATOR<TV> iterator(grid,2);iterator.Valid();iterator.Next()){
                INDEX cell_index=iterator.Cell_Index();
                if(regular_cell(cell_index) && cell_near_interface(cell_index)){
                    T weight=U_ghost(cell_index)(dim)*cell_volume - sigma(cell_index);
                    int index=weights.Append(PAIR<T,INDEX>(weight,cell_index));
                    donors(cell_index).Append(index); receivers(cell_index).Append(index);
                    sigma(cell_index) = U_ghost(cell_index)(dim)*cell_volume;}}

            // Compute backward-cast weights
            LINEAR_INTERPOLATION_UNIFORM<TV,T> linear;
            for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){
                INDEX cell_index=iterator.Cell_Index();
                if(!regular_cell(cell_index)){
                    RANGE<TV> cell_preimage=iterator.Bounding_Box() - dt*linear.Clamped_To_Array_Face(grid,face_velocities,iterator.Location());
                    RANGE<TV_INT> affected_cells(grid.Index(cell_preimage.min_corner),grid.Index(cell_preimage.max_corner)+TV_INT::All_Ones_Vector());
                    for(CELL_ITERATOR<TV> intersecting_iter(grid,affected_cells);intersecting_iter.Valid();intersecting_iter.Next()){
                        INDEX donor_cell=intersecting_iter.Cell_Index();
                        if(regular_cell(donor_cell)) continue;
                        T weight=U_ghost(donor_cell)(dim) * cell_preimage.Intersection_Area(intersecting_iter.Bounding_Box());
                        if(weight > (T)1e-14){
                            int index=weights.Append(PAIR<T,INDEX>(weight,cell_index));
                            sigma(donor_cell) += weight; donors(donor_cell).Append(index); receivers(cell_index).Append(index);}}}}

            // Clamp and forward-cast weights
            for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){
                INDEX cell_index=iterator.Cell_Index();
                if(!regular_cell(cell_index)){
                    T cell_stuff = U_ghost(cell_index)(dim) * cell_volume;
                    if(abs(sigma(cell_index)) > abs(cell_stuff)) {
                        T one_over_sigma = cell_stuff / sigma(cell_index);
                        for(int i=0;i<donors(cell_index).Size();++i){
                            weights(donors(cell_index)(i)).x *= one_over_sigma;}}
                    else {
                        T remainder = (cell_stuff - sigma(cell_index))/cell_volume;
                        RANGE<TV> cell_postimage=iterator.Bounding_Box() + dt*linear.Clamped_To_Array_Face(grid,face_velocities,iterator.Location());
                        RANGE<TV_INT> affected_cells(grid.Index(cell_postimage.min_corner),grid.Index(cell_postimage.max_corner)+TV_INT::All_Ones_Vector());
                        for(CELL_ITERATOR<TV> intersecting_iter(grid,affected_cells);intersecting_iter.Valid();intersecting_iter.Next()){
                            INDEX receiver_cell=intersecting_iter.Cell_Index();
                            T weight= remainder * cell_postimage.Intersection_Area(intersecting_iter.Bounding_Box());
                            int index=weights.Append(PAIR<T,INDEX>(weight,receiver_cell));
                            donors(cell_index).Append(index); receivers(receiver_cell).Append(index);}}}}

            for(int i=0;i<weights.Size();i++){
                PAIR<T,INDEX>& weight_pair(weights(i));
                if(U.Valid_Index(weight_pair.y)) U(weight_pair.y)(dim) += weight_pair.x;}}
    }

    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){INDEX cell_index=iterator.Cell_Index();
        if(cell_near_interface(cell_index)) U(cell_index) *= one_over_cell_volume;
        else if(psi(cell_index)) U(cell_index) = U_ghost(cell_index) - dt*rhs(cell_index);}
}
namespace PhysBAM{
template class HYBRID_SL_ENO_CONSERVATION<VECTOR<float,1>,3>;
template class HYBRID_SL_ENO_CONSERVATION<VECTOR<float,2>,4>;
template class HYBRID_SL_ENO_CONSERVATION<VECTOR<float,3>,5>;
template class HYBRID_SL_ENO_CONSERVATION<VECTOR<double,1>,3>;
template class HYBRID_SL_ENO_CONSERVATION<VECTOR<double,2>,4>;
template class HYBRID_SL_ENO_CONSERVATION<VECTOR<double,3>,5>;
}
