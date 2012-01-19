//#####################################################################
// Copyright 2010, Mridul Aanjaneya, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM.h>

using namespace PhysBAM;

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
ADVECTION_CONSERVATIVE_UNIFORM()
    :use_second_order(false),clamp_weights(true),use_collision_object(false),total_mass_lost(T2()),total_mass_gained(T2()),total_potential_lost(T2()),total_potential_gained(T2()),find_closest_point(false),num_cells(1),number_of_ghost_cells(3),cfl(0),density(1000),num_iterations(1),num_diffusion_iterations(1),pls(0),mpi_grid(0),evenodd(0),evenodd_cell(0),smoke_density(0)
{}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Node(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
    const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
    const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
{
    PHYSBAM_FATAL_ERROR();
}
template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Cell(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
    const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
    const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
{
    PHYSBAM_FATAL_ERROR();
    PHYSBAM_ASSERT(!Z_min || !Z_max); //we don't support extrema yet
    T_ARRAYS_VECTOR V_ghost(grid.Domain_Indices(number_of_ghost_cells));
    BOUNDARY_UNIFORM<T_GRID,TV> boundary_vel;
    boundary_vel.Fill_Ghost_Cells(grid,V,V_ghost,dt,time,number_of_ghost_cells);
    ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT> weights; //for each cell an array of weights and cells they are associated with
    ARRAY<T,TV_INT> sum;
    weights.Resize(grid.Domain_Indices(number_of_ghost_cells)),sum.Resize(grid.Domain_Indices(number_of_ghost_cells));
    T_INTERPOLATION interpolation;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,TV> linear;
    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        TV X=iterator.Location()-dt*V_ghost(cell);
        if(use_second_order) X+=(iterator.Location()-(X+dt*linear.Clamped_To_Array(grid,V_ghost,X)))/2.;
        ARRAY<PAIR<TV_INT,T> > backwards_weights=interpolation.Clamped_To_Array_Weights(grid,Z_ghost,X);
        for(int i=0;i<backwards_weights.m;i++) weights(backwards_weights(i).x).Append(PAIR<TV_INT,T>(cell,backwards_weights(i).y));}
    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        ARRAY<PAIR<TV_INT,T> >& local_weights=weights(cell);
        sum(cell)=0;for(int i=0;i<local_weights.m;i++) sum(cell)+=local_weights(i).y;
        if(clamp_weights) if(sum(cell)>1) for(int i=0;i<local_weights.m;i++) local_weights(i).y/=sum(cell);}
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) Z(iterator.Cell_Index())=T2();
    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        ARRAY<PAIR<TV_INT,T> >& local_weights=weights(cell);
        for(int i=0;i<local_weights.m;i++) if(grid.Domain_Indices().Lazy_Inside(local_weights(i).x)) Z(local_weights(i).x)+=local_weights(i).y*Z_ghost(cell);}
    T_ARRAYS_T2 Z_ghost2(grid.Domain_Indices(2*number_of_ghost_cells+1));
    boundary.Fill_Ghost_Cells(grid,Z,Z_ghost2,dt,time,2*number_of_ghost_cells+1);
    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(clamp_weights) if(sum(cell)>=1) continue;
        T2 remaining=(1-sum(cell))*Z_ghost(cell);
        TV X=iterator.Location()+dt*V_ghost(cell);
        if(use_second_order) X+=(iterator.Location()-(X-dt*linear.Clamped_To_Array(grid,V_ghost,X)))/2.;
        ARRAY<PAIR<TV_INT,T> > forward_weights=interpolation.Clamped_To_Array_Weights(grid,Z_ghost2,X);
        for(int i=0;i<forward_weights.m;i++) if(grid.Domain_Indices().Lazy_Inside(forward_weights(i).x)) Z(forward_weights(i).x)+=forward_weights(i).y*remaining;}
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Clamp_Weights_To_Objects(const GRID<TV>& grid,ARRAY<PAIR<TV_INT,T> >& weights)
{
    T delta=(T)1e-5;
    for(int i=0;i<weights.m;i++) if(ghost_box.Lazy_Inside(grid.Center(weights(i).x))) weights(i).y=0;
    T sum=0;for(int i=0;i<weights.m;i++) sum+=weights(i).y;
    if(sum>delta && sum!=1) for(int i=0;i<weights.m;i++) weights(i).y/=sum;
    assert(sum==0 || sum>delta);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Clamp_Weights_To_Objects(const GRID<TV>& grid,ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& weights)
{
    T delta=(T)1e-5;
    for(int i=0;i<weights.m;i++) if(ghost_box.Lazy_Inside(grid.Axis_X_Face(weights(i).x))) weights(i).y=0;
    T sum=0;for(int i=0;i<weights.m;i++) sum+=weights(i).y;
    if(sum>delta && sum!=1) for(int i=0;i<weights.m;i++) weights(i).y/=sum;
    assert(sum==0 || sum>delta);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Clamp_Weights_To_Grid(const RANGE<TV_INT>& inside_domain,ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& weights)
{
    T delta=(T)1e-5;
    for(int i=0;i<weights.m;i++) if(!inside_domain.Lazy_Inside(weights(i).x.index)) for(int j=0;j<TV::dimension;j++){
        TV_INT index=TV_INT::All_Ones_Vector()*2;index(j)=weights(i).x.index(j);
        if(!inside_domain.Lazy_Inside(index)){
            int side=1;if(index(j)>inside_domain.max_corner(j)) side=2;else assert(index(j)<inside_domain.min_corner(j));
            if(solid_walls_hack_axis(j)(side)) weights(i).y=0;}}
    T sum=0;for(int i=0;i<weights.m;i++) sum+=weights(i).y;
    if(sum>delta && sum!=1) for(int i=0;i<weights.m;i++) weights(i).y/=sum;
    assert(sum==0 || sum>delta);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Clamp_X_To_Grid(const GRID<TV>& grid,TV& X)
{
    bool outside=false;RANGE<TV> range=grid.Domain();
    if(!range.Lazy_Inside(X)) for(int i=0;i<TV::dimension;i++){TV point=range.min_corner;point(i)=X(i);
        if(!range.Lazy_Inside(point)){
            int side=1;if(point(i)>range.max_corner(i)) side=2;else assert(point(i)<range.min_corner(i));
            if(solid_walls_hack_axis(i)(side)) outside=true;}}
    if(outside) X=grid.Clamp(X);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> bool ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Is_Outside(const RANGE<TV_INT>& inside_domain,const FACE_INDEX<TV::dimension>& face)
{
    if(!inside_domain.Lazy_Inside(face.index)) for(int i=0;i<TV::dimension;i++){
        TV_INT index=TV_INT::All_Ones_Vector()*2;index(i)=face.index(i);
        if(!inside_domain.Lazy_Inside(index)){
            int side=1;if(index(i)>inside_domain.max_corner(i)) side=2;else assert(index(i)<inside_domain.min_corner(i));
            if(solid_walls_hack_axis(i)(side)) return true;}}
    return false;
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> bool ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Is_MPI_Boundary(const RANGE<TV_INT>& inside_domain,const FACE_INDEX<TV::dimension>& face)
{
    return Is_MPI_Boundary(inside_domain,face.index);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> bool ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Is_MPI_Boundary(const RANGE<TV_INT>& inside_domain,const TV_INT& index)
{
    if(!inside_domain.Lazy_Inside(index)) for(int i=0;i<TV::dimension;i++){
        TV_INT tmp_index=TV_INT::All_Ones_Vector()*2;tmp_index(i)=index(i);
        if(!inside_domain.Lazy_Inside(tmp_index)){
            int side=1;if(tmp_index(i)>inside_domain.max_corner(i)) side=2;else assert(tmp_index(i)<inside_domain.min_corner(i));
            if(mpi_boundary(i)(side)) return true;}}
    return false;
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Clean_Weights(ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& weights)
{
    bool rescale=false;T delta=(T)1e-4;
    for(int i=0;i<weights.m;i++){
        assert(weights(i).y>-delta);
        if(weights(i).y<delta){weights(i).y=0;rescale=true;}}
    if(rescale){
        T sum=0;for(int i=0;i<weights.m;i++) sum+=weights(i).y;
        for(int i=0;i<weights.m;i++) weights(i).y/=sum;}
}

/*template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Cell_Diffusion_Helper(FACE_ITERATOR& iterator,ARRAY<T,TV_INT>& sum_jc_cell,T_ARRAYS_T2& Z,ARRAY<T,TV_INT>& give_jc,T_ARRAYS_T2& give_Z)
{
    T wjc_diff=(sum_jc_cell(iterator.First_Cell_Index())-sum_jc_cell(iterator.Second_Cell_Index()))/2.;
    T Z_diff=wjc_diff>0?wjc_diff/sum_jc_cell(iterator.First_Cell_Index()):wjc_diff/sum_jc_cell(iterator.Second_Cell_Index());
    sum_jc_cell(iterator.Second_Cell_Index())+=wjc_diff;sum_jc_cell(iterator.First_Cell_Index())-=wjc_diff;
    T2 local_Z=wjc_diff>0?Z(iterator.First_Cell_Index()):Z(iterator.Second_Cell_Index());
    Z(iterator.Second_Cell_Index())+=local_Z*(Z_diff);Z(iterator.First_Cell_Index())-=local_Z*(Z_diff);
    if(!grid.Domain_Indices().Lazy_Inside(iterator.First_Cell_Index())){give_jc(iterator.First_Cell_Index())-=wjc_diff;give_Z(iterator.First_Cell_Index())-=local_Z*(Z_diff)}
    if(!grid.Domain_Indices().Lazy_Inside(iterator.Second_Cell_Index())){give_jc(iterator.Second_Cell_Index())+=wjc_diff;give_Z(iterator.Second_Cell_Index())+=local_Z*(Z_diff)}
}*/

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Cell_Diffusion_Helper(FACE_ITERATOR& iterator,ARRAY<T,TV_INT>& sum_jc_cell,T_ARRAYS_T2& Z)
{
    if(ghost_box.Lazy_Inside(iterator.grid.Center(iterator.First_Cell_Index())) || ghost_box.Lazy_Inside(iterator.grid.Center(iterator.Second_Cell_Index()))) return;
    T wjc_diff=(T)((sum_jc_cell(iterator.First_Cell_Index())-sum_jc_cell(iterator.Second_Cell_Index()))/2.);
    T Z_diff=wjc_diff>0?wjc_diff/sum_jc_cell(iterator.First_Cell_Index()):wjc_diff/sum_jc_cell(iterator.Second_Cell_Index());
    sum_jc_cell(iterator.Second_Cell_Index())+=wjc_diff;sum_jc_cell(iterator.First_Cell_Index())-=wjc_diff;
    T2 local_Z=wjc_diff>0?Z(iterator.First_Cell_Index()):Z(iterator.Second_Cell_Index());
    Z(iterator.Second_Cell_Index())+=local_Z*(Z_diff);Z(iterator.First_Cell_Index())-=local_Z*(Z_diff);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Face_Diffusion_Helper(const GRID<TV>& grid,FACE_INDEX<TV::dimension>& first_face_index,FACE_INDEX<TV::dimension>& second_face_index,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside)
{ 
    if(inside && !((*inside)(first_face_index) && (*inside)(second_face_index))) return;
    if(ghost_box.Lazy_Inside(grid.Axis_X_Face(first_face_index)) || ghost_box.Lazy_Inside(grid.Axis_X_Face(second_face_index))) return;
    T wjc_diff=(T)((sum_jc(first_face_index)-sum_jc(second_face_index))/2.);
    T Z_diff=wjc_diff>0?wjc_diff/sum_jc(first_face_index):wjc_diff/sum_jc(second_face_index);
    sum_jc(second_face_index)+=wjc_diff;sum_jc(first_face_index)-=wjc_diff;
    T local_Z=wjc_diff>0?Z(first_face_index):Z(second_face_index);
    Z(second_face_index)+=local_Z*(Z_diff);Z(first_face_index)-=local_Z*(Z_diff);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Face_Diffusion_Helper(FACE_ITERATOR& iterator,int axis,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside)
{
    assert(axis!=iterator.Axis());
    FACE_INDEX<TV::dimension> first_face_index=FACE_INDEX<TV::dimension>(axis,iterator.First_Cell_Index()),second_face_index=FACE_INDEX<TV::dimension>(axis,iterator.Second_Cell_Index());
    Face_Diffusion_Helper(iterator.grid,first_face_index,second_face_index,sum_jc,Z,inside);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Face_Diffusion_Helper(CELL_ITERATOR& iterator,int axis,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside)
{
    FACE_INDEX<TV::dimension> first_face_index=FACE_INDEX<TV::dimension>(axis,iterator.First_Face_Index(axis)),second_face_index=FACE_INDEX<TV::dimension>(axis,iterator.Second_Face_Index(axis));
    Face_Diffusion_Helper(iterator.grid,first_face_index,second_face_index,sum_jc,Z,inside);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Cell_Diffusion(const T_GRID& grid,ARRAY<T,TV_INT>& sum_jc_cell,T_ARRAYS_T2& Z,T_BOUNDARY_T2& boundary,BOUNDARY_UNIFORM<T_GRID,T>* boundary_sum)
{
    for(int iter=0;iter<num_diffusion_iterations;iter++){
        for(int axis=0;axis<TV::dimension;axis++){
            if(evenodd_cell==0){
                RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner-=TV_INT::All_Ones_Vector();domain.min_corner+=TV_INT::All_Ones_Vector();domain.max_corner(axis)+=1;
                for(FACE_ITERATOR iterator(grid,domain,axis);iterator.Valid();iterator.Next()) Cell_Diffusion_Helper(iterator,sum_jc_cell,Z);}
            else{
                ARRAY<FACE_ITERATOR*> faces;
                RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner-=TV_INT::All_Ones_Vector();domain.min_corner+=TV_INT::All_Ones_Vector();domain.max_corner(axis)+=1;
                for(FACE_ITERATOR iterator(grid,domain,axis);iterator.Valid();iterator.Next()) faces.Append(new FACE_ITERATOR(iterator));
                for(int i=faces.m;i>=1;i--){FACE_ITERATOR& iterator=*faces(i);Cell_Diffusion_Helper(iterator,sum_jc_cell,Z);}
                for(int i=faces.m;i>=1;i--) delete faces(i);}
            if(mpi_grid){
                T_ARRAYS_T2 Z_ghost_local(grid.Domain_Indices(number_of_ghost_cells));
                ARRAY<T,TV_INT> sum_jc_cell_ghost(grid.Domain_Indices(number_of_ghost_cells));
                boundary.Fill_Ghost_Cells(grid,Z,Z_ghost_local,0,0,number_of_ghost_cells); //Sync for mpi boundaries
                boundary_sum->Fill_Ghost_Cells(grid,sum_jc_cell,sum_jc_cell_ghost,0,0,number_of_ghost_cells); //Sync for mpi boundaries
                for(int side=0;side<2;side++) if(mpi_boundary(axis)(side)){
                    if(evenodd_cell==0){
                        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner(axis)+=1;
                        if(side==1) domain.max_corner(axis)=domain.min_corner(axis);else domain.min_corner(axis)=domain.max_corner(axis);
                        for(FACE_ITERATOR iterator(grid,domain,axis);iterator.Valid();iterator.Next()) Cell_Diffusion_Helper(iterator,sum_jc_cell_ghost,Z_ghost_local);}
                    else{
                        ARRAY<FACE_ITERATOR*> faces;
                        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner(axis)+=1;
                        if(side==1) domain.max_corner(axis)=domain.min_corner(axis);else domain.min_corner(axis)=domain.max_corner(axis);
                        for(FACE_ITERATOR iterator(grid,domain,axis);iterator.Valid();iterator.Next()) faces.Append(new FACE_ITERATOR(iterator));
                        for(int i=faces.m;i>=1;i--){FACE_ITERATOR& iterator=*faces(i);Cell_Diffusion_Helper(iterator,sum_jc_cell_ghost,Z_ghost_local);}
                        for(int i=faces.m;i>=1;i--) delete faces(i);}}
                for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) Z(iterator.Cell_Index())=Z_ghost_local(iterator.Cell_Index());
                for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) sum_jc_cell(iterator.Cell_Index())=sum_jc_cell_ghost(iterator.Cell_Index());}}
        evenodd_cell++;evenodd_cell=evenodd_cell%2;}
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Face_Diffusion(const T_GRID& grid,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,T_BOUNDARY& boundary,BOUNDARY_UNIFORM<T_GRID,T>* boundary_sum,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Before diffusion",0,0);
    for(int iter=0;iter<num_diffusion_iterations;iter++){
        for(int axis=0;axis<TV::dimension;axis++){
            if(evenodd==0){
                for(int axis2=0;axis2<TV::dimension;axis2++){if(axis==axis2) continue;//handled above
                    RANGE<TV_INT> domain=grid.Domain_Indices();
                    domain.max_corner(axis)+=1;domain.min_corner(axis2)+=1;
                    if(solid_walls_hack_axis(axis)(1)) domain.min_corner(axis)+=1;
                    if(solid_walls_hack_axis(axis)(2)) domain.max_corner(axis)-=1;
                    for(FACE_ITERATOR iterator(grid,domain,axis2);iterator.Valid();iterator.Next()) Face_Diffusion_Helper(iterator,axis,sum_jc,Z,inside);}
                RANGE<TV_INT> domain=grid.Domain_Indices();
                if((mpi_grid && mpi_boundary(axis)(1)) || solid_walls_hack_axis(axis)(1)) domain.min_corner(axis)+=1; 
                if((mpi_grid && mpi_boundary(axis)(2)) || solid_walls_hack_axis(axis)(2)) domain.max_corner(axis)-=1; 
                for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()) Face_Diffusion_Helper(iterator,axis,sum_jc,Z,inside);}
            else{
                ARRAY<FACE_ITERATOR*> faces;ARRAY<CELL_ITERATOR*> cells;
                for(int axis2=0;axis2<TV::dimension;axis2++){if(axis==axis2) continue;//handled above
                    RANGE<TV_INT> domain=grid.Domain_Indices();
                    domain.max_corner(axis)+=1;domain.min_corner(axis2)+=1;
                    if(solid_walls_hack_axis(axis)(1)) domain.min_corner(axis)+=1;
                    if(solid_walls_hack_axis(axis)(2)) domain.max_corner(axis)-=1;
                    for(FACE_ITERATOR iterator(grid,domain,axis2);iterator.Valid();iterator.Next()) faces.Append(new FACE_ITERATOR(iterator));}
                RANGE<TV_INT> domain=grid.Domain_Indices();
                if((mpi_grid && mpi_boundary(axis)(1)) || solid_walls_hack_axis(axis)(1)) domain.min_corner(axis)+=1; 
                if((mpi_grid && mpi_boundary(axis)(2)) || solid_walls_hack_axis(axis)(2)) domain.max_corner(axis)-=1; 
                for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()) cells.Append(new CELL_ITERATOR(iterator));
                for(int i=faces.m;i>=1;i--){FACE_ITERATOR& iterator=*faces(i);Face_Diffusion_Helper(iterator,axis,sum_jc,Z,inside);}
                for(int i=cells.m;i>=1;i--){CELL_ITERATOR& iterator=*cells(i);Face_Diffusion_Helper(iterator,axis,sum_jc,Z,inside);}
                for(int i=faces.m;i>=1;i--) delete faces(i);
                for(int i=cells.m;i>=1;i--) delete cells(i);}
            if(mpi_grid){
                T_FACE_ARRAYS_SCALAR Z_ghost_local(grid,number_of_ghost_cells);
                ARRAY<T,FACE_INDEX<TV::dimension> > sum_jc_ghost(grid,number_of_ghost_cells);
                boundary.Apply_Boundary_Condition_Face(grid,Z,0);
                boundary.Fill_Ghost_Cells_Face(grid,Z,Z_ghost_local,0,number_of_ghost_cells); //Sync for mpi boundaries
                boundary_sum->Apply_Boundary_Condition_Face(grid,sum_jc,0);
                boundary_sum->Fill_Ghost_Cells_Face(grid,sum_jc,sum_jc_ghost,0,number_of_ghost_cells); //Sync for mpi boundaries
                for(int side=0;side<2;side++){
                    if(evenodd==0){
                        for(int axis2=0;axis2<TV::dimension;axis2++){if(axis==axis2 || !mpi_boundary(axis2)(side)) continue;//handled above
                            RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner(axis)++;
                            if(side==1) domain.max_corner(axis2)=domain.min_corner(axis2);
                            else{domain.max_corner(axis2)++;domain.min_corner(axis2)=domain.max_corner(axis2);}
                            for(FACE_ITERATOR iterator(grid,domain,axis2);iterator.Valid();iterator.Next()) Face_Diffusion_Helper(iterator,axis,sum_jc_ghost,Z_ghost_local,inside);}
                        RANGE<TV_INT> domain=grid.Domain_Indices(1);
                        if(side==1) domain.max_corner(axis)=domain.min_corner(axis)+1;
                        else domain.min_corner(axis)=domain.max_corner(axis)-1;
                        if(mpi_boundary(axis)(side)) for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()) Face_Diffusion_Helper(iterator,axis,sum_jc_ghost,Z_ghost_local,inside);}
                    else{
                        ARRAY<FACE_ITERATOR*> faces;ARRAY<CELL_ITERATOR*> cells;
                        for(int axis2=0;axis2<TV::dimension;axis2++){if(axis==axis2 || !mpi_boundary(axis2)(side)) continue;//handled above
                            RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner(axis)++;
                            if(side==1) domain.max_corner(axis2)=domain.min_corner(axis2);
                            else{domain.max_corner(axis2)++;domain.min_corner(axis2)=domain.max_corner(axis2);}
                            for(FACE_ITERATOR iterator(grid,domain,axis2);iterator.Valid();iterator.Next()) faces.Append(new FACE_ITERATOR(iterator));}
                        RANGE<TV_INT> domain=grid.Domain_Indices(1);
                        if(side==1) domain.max_corner(axis)=domain.min_corner(axis)+1;
                        else domain.min_corner(axis)=domain.max_corner(axis)-1;
                        if(mpi_boundary(axis)(side)) for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()) cells.Append(new CELL_ITERATOR(iterator));
                        for(int i=faces.m;i>=1;i--){FACE_ITERATOR& iterator=*faces(i);Face_Diffusion_Helper(iterator,axis,sum_jc_ghost,Z_ghost_local,inside);}
                        for(int i=cells.m;i>=1;i--){CELL_ITERATOR& iterator=*cells(i);Face_Diffusion_Helper(iterator,axis,sum_jc_ghost,Z_ghost_local,inside);}
                        for(int i=faces.m;i>=1;i--) delete faces(i);
                        for(int i=cells.m;i>=1;i--) delete cells(i);}}
                for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) Z(iterator.Full_Index())=Z_ghost_local(iterator.Full_Index());
                for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) sum_jc(iterator.Full_Index())=sum_jc_ghost(iterator.Full_Index());
                boundary.Apply_Boundary_Condition_Face(grid,Z,0);
                boundary_sum->Apply_Boundary_Condition_Face(grid,sum_jc,0);}}
        evenodd++;evenodd=evenodd%2;}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After diffusion",0,0);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
    const T_FACE_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
    const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
{
    static bool first_cell=true;
    PHYSBAM_ASSERT(!Z_min || !Z_max); //we don't support extrema yet
    T_FACE_ARRAYS_SCALAR face_velocities_ghost(grid,2*number_of_ghost_cells+1);
    BOUNDARY_UNIFORM<T_GRID,T> boundary_vel_scalar;boundary_vel_scalar.Set_Fixed_Boundary(true,0);
    BOUNDARY_UNIFORM<T_GRID,T>* boundary_vel=&boundary_vel_scalar;
    if(mpi_grid) boundary_vel=new BOUNDARY_MPI<GRID<TV> >(mpi_grid,boundary_vel_scalar);
    boundary_vel->Fill_Ghost_Cells_Face(grid,face_velocities.V_face,face_velocities_ghost,time,2*number_of_ghost_cells+1);
    T_FACE_LOOKUP face_velocities_ghost_lookup(face_velocities_ghost);
    BOUNDARY_UNIFORM<T_GRID,T> boundary_scalar;boundary_scalar.Set_Fixed_Boundary(true,1);
    BOUNDARY_UNIFORM<T_GRID,T>* boundary_sum=&boundary_scalar;
    if(mpi_grid) boundary_sum=new BOUNDARY_MPI<GRID<TV> >(mpi_grid,boundary_scalar);
    ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT> weights_to; //w ij
    ARRAY<ARRAY<PAIR<TV_INT,int> >,TV_INT> weights_from; //w ki - points to weights in weights to
    ARRAY<T,TV_INT> sum;
    weights_to.Resize(grid.Domain_Indices(number_of_ghost_cells)),weights_from.Resize(grid.Domain_Indices(number_of_ghost_cells)),sum.Resize(grid.Domain_Indices()),sum_jc_cell.Resize(grid.Domain_Indices());
    if(first_cell && time<=dt) for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) sum_jc_cell(iterator.Cell_Index())=1;
    first_cell=false;
    T_INTERPOLATION interpolation;T_AVERAGING averaging;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> linear;
    RANGE<TV_INT> ghost_domain=grid.Domain_Indices(number_of_ghost_cells),real_domain=grid.Domain_Indices();
    RANGE<TV> real_domain_x=grid.domain;real_domain_x.min_corner-=grid.dX;real_domain_x.max_corner+=grid.dX;
    for(CELL_ITERATOR iterator(grid,cfl);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(ghost_box.Lazy_Inside(iterator.Location())) continue;
        TV X=iterator.Location()-dt*averaging.Face_To_Cell_Vector(grid,cell,face_velocities_ghost_lookup);
        if(use_second_order) X+=(iterator.Location()-(X+dt*linear.Clamped_To_Array_Face(grid,face_velocities_ghost_lookup,X)))/2.;
        if(ghost_box.Lazy_Inside(X)) X=ghost_box.Surface(X);        
        if(!real_domain.Lazy_Inside(cell) && !real_domain_x.Lazy_Inside(X)) continue;
        ARRAY<PAIR<TV_INT,T> > backwards_weights=interpolation.Clamped_To_Array_Weights(grid,Z_ghost,X);
        Clamp_Weights_To_Objects(grid,backwards_weights);
        for(int i=0;i<backwards_weights.m;i++){assert(backwards_weights(i).y>-1e-5);if(backwards_weights(i).y<0) backwards_weights(i).y=0;}
        for(int i=0;i<backwards_weights.m;i++){if(!ghost_domain.Lazy_Inside(backwards_weights(i).x)){assert(!real_domain.Lazy_Inside(cell));continue;}
            weights_to(backwards_weights(i).x).Append(PAIR<TV_INT,T>(cell,backwards_weights(i).y));
            weights_from(cell).Append(PAIR<TV_INT,int>(backwards_weights(i).x,weights_to(backwards_weights(i).x).m));}}
    if(mpi_grid) mpi_grid->Sync_Common_Cell_Weights_To(weights_to,weights_from,number_of_ghost_cells);
    LOG::Time("After SL");

    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        ARRAY<PAIR<TV_INT,T> >& local_weights=weights_to(cell);
        sum(cell)=0;for(int i=0;i<local_weights.m;i++) sum(cell)+=local_weights(i).y;}

    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(clamp_weights) if(sum(cell)>=1) continue;
        T remaining=(1-sum(cell));
        TV X=iterator.Location()+dt*averaging.Face_To_Cell_Vector(grid,cell,face_velocities_ghost_lookup);
        if(use_second_order) X+=(iterator.Location()-(X-dt*linear.Clamped_To_Array_Face(grid,face_velocities_ghost_lookup,X)))/2.;
        if(!ghost_box.Lazy_Inside(iterator.Location()) && ghost_box.Lazy_Inside(X)) X=ghost_box.Surface(X);        
        ARRAY<PAIR<TV_INT,T> > forward_weights=interpolation.Clamped_To_Array_Weights(grid,Z_ghost,X);
        Clamp_Weights_To_Objects(grid,forward_weights);
        for(int i=0;i<forward_weights.m;i++){assert(forward_weights(i).y>-1e-5);if(forward_weights(i).y<0) forward_weights(i).y=0;}
        for(int i=0;i<forward_weights.m;i++){
            int index=0;for(int j=0;j<weights_to(cell).m;j++) if(weights_to(cell)(j).x==forward_weights(i).x) index=j;
            if(index) weights_to(cell)(index).y+=forward_weights(i).y*remaining;
            else{
                weights_to(cell).Append(PAIR<TV_INT,T>(forward_weights(i).x,forward_weights(i).y*remaining));
                weights_from(forward_weights(i).x).Append(PAIR<TV_INT,int>(cell,weights_to(cell).m));}}}
    if(mpi_grid) mpi_grid->Sync_Common_Cell_Weights_From(weights_to,weights_from,number_of_ghost_cells);
    LOG::Time("After 2");

    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(ghost_box.Lazy_Inside(iterator.Location())) continue;
        ARRAY<PAIR<TV_INT,T> >& local_weights=weights_to(cell);
        T sum=0;for(int i=0;i<local_weights.m;i++) sum+=local_weights(i).y;
        for(int i=0;i<local_weights.m;i++){assert(sum>1e-6);
            local_weights(i).y=local_weights(i).y/sum;}}
    if(mpi_grid) mpi_grid->Sync_Common_Cell_Weights_From(weights_to,weights_from,number_of_ghost_cells);
    LOG::Time("After Clamp 1");

    ARRAY<T,TV_INT> sum_ic(grid.Domain_Indices(number_of_ghost_cells));
    boundary_sum->Fill_Ghost_Cells(grid,sum_jc_cell,sum_ic,dt,time,number_of_ghost_cells);
    for(int i=0;i<num_iterations;i++){
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            if(ghost_box.Lazy_Inside(iterator.Location())) continue;            
            ARRAY<PAIR<TV_INT,int> >& local_weights=weights_from(cell);
            T sum=0;for(int i=0;i<local_weights.m;i++) sum+=sum_ic(local_weights(i).x)*weights_to(local_weights(i).x)(local_weights(i).y).y;
            for(int i=0;i<local_weights.m;i++){assert(sum_jc_cell(cell)>1e-6);
                weights_to(local_weights(i).x)(local_weights(i).y).y=weights_to(local_weights(i).x)(local_weights(i).y).y/sum;}}
        if(mpi_grid) mpi_grid->Sync_Common_Cell_Weights_To(weights_to,weights_from,number_of_ghost_cells);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            if(ghost_box.Lazy_Inside(iterator.Location())) continue;            
            ARRAY<PAIR<TV_INT,T> >& local_weights=weights_to(cell);
            T sum=0;for(int i=0;i<local_weights.m;i++) sum+=local_weights(i).y;
            for(int i=0;i<local_weights.m;i++){assert(sum>1e-5);
                local_weights(i).y=local_weights(i).y/sum;}}
        if(mpi_grid) mpi_grid->Sync_Common_Cell_Weights_From(weights_to,weights_from,number_of_ghost_cells);}
    LOG::Time("After Clamp final");

    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(ghost_box.Lazy_Inside(iterator.Location())){sum_jc_cell(cell)=1;continue;}
        ARRAY<PAIR<TV_INT,int> >& local_weights=weights_from(cell);
        sum_jc_cell(cell)=0;for(int i=0;i<local_weights.m;i++) sum_jc_cell(cell)+=sum_ic(local_weights(i).x)*weights_to(local_weights(i).x)(local_weights(i).y).y;}
    LOG::Time("After wjc calc");

    T_ARRAYS_T2 total_mass_lost_per_cell(grid.Domain_Indices(2*number_of_ghost_cells+1)),total_mass_gained_per_cell(grid.Domain_Indices(2*number_of_ghost_cells+1));
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(ghost_box.Lazy_Inside(iterator.Location())) continue;
        ARRAY<PAIR<TV_INT,T> >& local_weights=weights_to(cell);
        sum(cell)=0;for(int i=0;i<local_weights.m;i++) sum(cell)+=local_weights(i).y;}
    LOG::Time("After wj none");

    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){if(ghost_box.Lazy_Inside(iterator.Location())) continue;Z(iterator.Cell_Index())=T2();}
    for(CELL_ITERATOR iterator(grid,cfl);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        ARRAY<PAIR<TV_INT,T> >& local_weights=weights_to(cell);
        //if((!grid.Domain_Indices().Lazy_Inside(cell) && !Is_MPI_Boundary(grid.Domain_Indices(),cell)) || ghost_box.Lazy_Inside(iterator.Location())){
        //    for(int i=0;i<local_weights.m;i++) total_mass_gained_per_cell(local_weights(i).x)+=local_weights(i).y*Z_ghost(cell);
        //    T sum=0;for(int i=0;i<local_weights.m;i++) sum+=local_weights(i).y;total_mass_gained+=sum*Z_ghost(cell);}
        for(int i=0;i<local_weights.m;i++)
            if(grid.Domain_Indices().Lazy_Inside(local_weights(i).x) && !ghost_box.Lazy_Inside(grid.Center(local_weights(i).x))) Z(local_weights(i).x)+=local_weights(i).y*Z_ghost(cell);}
            //else if(!Is_MPI_Boundary(grid.Domain_Indices(),local_weights(i).x)){
            //    total_mass_lost_per_cell(local_weights(i).x)+=local_weights(i).y*Z_ghost(cell);
            //    total_mass_lost+=local_weights(i).y*Z_ghost(cell);}}
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) total_potential_gained+=total_mass_gained_per_cell(iterator.Cell_Index())*total_mass_gained_per_cell(iterator.Cell_Index())*(1-iterator.Location()(2));
    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next()) total_potential_lost+=total_mass_lost_per_cell(iterator.Cell_Index())*total_mass_lost_per_cell(iterator.Cell_Index())*(1-iterator.Location()(2));
    LOG::Time("After distrib calc");

    //diffusion
    Cell_Diffusion(grid,sum_jc_cell,Z,boundary,boundary_sum);
    if(mpi_grid) delete boundary_sum;
    LOG::Time("After diffusion");
}
 
template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Face_Lookup(T_GRID& grid,T_ARRAYS_SCALAR& phi1,T_ARRAYS_SCALAR& phi2,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
    const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
    const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
{
    static bool first=true;
    T delta=(T)1e-5;
    T_LEVELSET lsv1(grid,phi1),lsv2(grid,phi2);
    T_FACE_ARRAYS_BOOL inside1(grid,2*number_of_ghost_cells+1),inside2(grid,2*number_of_ghost_cells+1);
    BOUNDARY_UNIFORM<T_GRID,T> boundary_scalar;boundary_scalar.Set_Fixed_Boundary(true,1);
    BOUNDARY_UNIFORM<T_GRID,T>* boundary_sum=&boundary_scalar;
    if(mpi_grid) boundary_sum=new BOUNDARY_MPI<GRID<TV> >(mpi_grid,boundary_scalar);
    for(FACE_ITERATOR iterator(grid,2*number_of_ghost_cells+1);iterator.Valid();iterator.Next()){
        inside1(iterator.Full_Index())=false;
        inside2(iterator.Full_Index())=false;}
    for(CELL_ITERATOR iterator(grid,2*number_of_ghost_cells+1);iterator.Valid();iterator.Next()){
        for(int axis=0;axis<TV::dimension;axis++){
            if(!inside1(FACE_INDEX<TV::dimension>(axis,iterator.First_Face_Index(axis)))) inside1(FACE_INDEX<TV::dimension>(axis,iterator.First_Face_Index(axis)))=lsv1.Phi(iterator.Location())<=0;
            if(!inside1(FACE_INDEX<TV::dimension>(axis,iterator.Second_Face_Index(axis)))) inside1(FACE_INDEX<TV::dimension>(axis,iterator.Second_Face_Index(axis)))=lsv1.Phi(iterator.Location())<=0;
            if(!inside2(FACE_INDEX<TV::dimension>(axis,iterator.First_Face_Index(axis)))) inside2(FACE_INDEX<TV::dimension>(axis,iterator.First_Face_Index(axis)))=lsv2.Phi(iterator.Location())<=0;
            if(!inside2(FACE_INDEX<TV::dimension>(axis,iterator.Second_Face_Index(axis)))) inside2(FACE_INDEX<TV::dimension>(axis,iterator.Second_Face_Index(axis)))=lsv2.Phi(iterator.Location())<=0;}}
    weights_to.Clean_Memory();weights_from.Clean_Memory();
    PHYSBAM_ASSERT(!Z_min || !Z_max); //we don't support extrema yet
    ARRAY<T,FACE_INDEX<TV::dimension> > sum;
    weights_to.Resize(grid,number_of_ghost_cells),weights_from.Resize(grid,number_of_ghost_cells),sum.Resize(grid),sum_jc.Resize(grid),momentum_lost.Resize(grid.Domain_Indices());
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) momentum_lost(iterator.Cell_Index())=0;
    if(first && time<=dt) for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) sum_jc(iterator.Full_Index())=1;
    first=false;
    T_INTERPOLATION interpolation;T_AVERAGING averaging;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> linear;
    int ghost_cells=number_of_ghost_cells;
    RANGE<TV> real_domain=grid.domain;real_domain.min_corner-=grid.dX*(1+delta);real_domain.max_corner+=grid.dX*(1+delta);
    for(FACE_ITERATOR iterator(grid,cfl);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(!inside2(iterator.Full_Index())) continue; //Don't get velocity if not inside lsv
        if(ghost_box.Lazy_Inside(iterator.Location())) continue;
        RANGE<TV_INT> ghost_domain=grid.Domain_Indices(number_of_ghost_cells);ghost_domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        RANGE<TV_INT> inside_domain=grid.Domain_Indices();inside_domain.min_corner+=TV_INT::Axis_Vector(face.axis);
        if(Is_Outside(inside_domain,face)) continue;
        TV X=iterator.Location()-dt*averaging.Face_To_Face_Vector(grid,face.axis,face.index,face_velocities);
        if(!domain.Lazy_Inside(face.index) && !real_domain.Lazy_Inside(X)) continue;
        if(use_second_order) X+=(iterator.Location()-(X+dt*linear.Clamped_To_Array_Face(grid,face_velocities,X)))/2.;
        Clamp_X_To_Grid(grid,X);
        if(ghost_box.Lazy_Inside(X)) X=ghost_box.Surface(X);
        if(use_collision_object && collision_object.Lazy_Inside(X)) X=collision_object.Surface(X);
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > backwards_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X);
        for(int i=0;i<backwards_weights.m;i++) if(!inside1(backwards_weights(i).x)) backwards_weights(i).y=0;
        Clamp_Weights_To_Grid(inside_domain,backwards_weights);
        Clamp_Weights_To_Objects(grid,backwards_weights);   
        T sum=0;for(int i=0;i<backwards_weights.m;i++) sum+=backwards_weights(i).y;
        assert(sum<=1+delta);
        if(sum<1 && sum>=delta) for(int i=0;i<backwards_weights.m;i++) backwards_weights(i).y=backwards_weights(i).y/sum;
        else if(sum<delta){
            TV axis_vector=TV::Axis_Vector(iterator.Axis());
            TV_INT cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
            TV X_center1=grid.Center(cell1)-dt*averaging.Face_To_Cell_Vector(grid,cell1,face_velocities),X_center2=grid.Center(cell2)-dt*averaging.Face_To_Cell_Vector(grid,cell2,face_velocities);
            TV X_cell1=X_center1+grid.dX/2.*axis_vector,X_cell2=X_center2-grid.dX/2.*axis_vector;
            assert(lsv2.Phi(grid.Center(cell1))<=0||lsv2.Phi(grid.Center(cell2))<=0);
            assert(lsv1.Phi(X_center1)<=0||lsv1.Phi(X_center2)<=0);
            if(lsv1.Phi(X_cell1)>0){
                if(lsv1.Phi(X_center1)<0) X_cell1=X_center1+grid.dX/2.*axis_vector*lsv1.Phi(X_center1)/(lsv1.Phi(X_center1)-lsv1.Phi(X_cell1));
                else if(lsv1.Phi(X_center1)<lsv1.Phi(X_cell1)) X_cell1=X_center1;}
            if(lsv1.Phi(X_cell2)>0){
                if(lsv1.Phi(X_center2)<0) X_cell2=X_center2-grid.dX/2.*axis_vector*lsv1.Phi(X_center2)/(lsv1.Phi(X_center2)-lsv1.Phi(X_cell2));
                else if(lsv1.Phi(X_center2)<lsv1.Phi(X_cell2)) X_cell2=X_center2;}
            assert(!use_second_order);
            Clamp_X_To_Grid(grid,X_cell1);
            Clamp_X_To_Grid(grid,X_cell2);
            if(use_collision_object && collision_object.Lazy_Inside(X_cell1)) X_cell1=collision_object.Surface(X_cell1);
            if(use_collision_object && collision_object.Lazy_Inside(X_cell2)) X_cell2=collision_object.Surface(X_cell2);
            backwards_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X_cell1);
            ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > backwards_weights2=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X_cell2);
            for(int i=0;i<backwards_weights.m;i++) backwards_weights(i).y*=0.5;
            for(int i=0;i<backwards_weights2.m;i++) backwards_weights2(i).y*=0.5;
            for(int i=0;i<backwards_weights2.m;i++){
                int index=0;for(int j=0;j<backwards_weights.m;j++) if(backwards_weights(j).x==backwards_weights2(i).x) index=j;
                if(index) backwards_weights(index).y+=backwards_weights2(i).y;
                else backwards_weights.Append(backwards_weights2(i));}
            for(int i=0;i<backwards_weights.m;i++) if(!inside1(backwards_weights(i).x)) backwards_weights(i).y=0;
            Clamp_Weights_To_Grid(inside_domain,backwards_weights); 
            Clamp_Weights_To_Objects(grid,backwards_weights);   
            sum=0;for(int i=0;i<backwards_weights.m;i++) sum+=backwards_weights(i).y;
            assert(sum<=1+delta);
            if(sum<1 && sum>=delta) for(int i=0;i<backwards_weights.m;i++) backwards_weights(i).y=backwards_weights(i).y/sum;
            else if(sum<delta){
                if(domain.Lazy_Inside(face.index)) PHYSBAM_FATAL_ERROR("Cannot find backward weights for this cell");
                else for(int i=0;i<backwards_weights.m;i++) backwards_weights(i).y=0;}} //It's ok to not get anything if this is a ghost cell
        sum=0;for(int i=0;i<backwards_weights.m;i++) sum+=backwards_weights(i).y;
        assert(!domain.Lazy_Inside(face.index)||(sum>delta));
        for(int i=0;i<backwards_weights.m;i++){//assert(backwards_weights(i).y>-delta);
            if(backwards_weights(i).y<0) backwards_weights(i).y=0;}
        for(int i=0;i<backwards_weights.m;i++){if(backwards_weights(i).y==0) continue;
            if(!ghost_domain.Lazy_Inside(backwards_weights(i).x.index)){assert(!domain.Lazy_Inside(face.index));continue;}
            weights_to(backwards_weights(i).x).Append(PAIR<FACE_INDEX<TV::dimension>,T>(face,backwards_weights(i).y));
            weights_from(face).Append(PAIR<FACE_INDEX<TV::dimension>,int>(backwards_weights(i).x,weights_to(backwards_weights(i).x).m));}}
    if(mpi_grid) mpi_grid->Sync_Common_Face_Weights_To(weights_to,weights_from,ghost_cells);

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
        if(ghost_box.Lazy_Inside(iterator.Location())) continue;
        sum(face)=0;for(int i=0;i<local_weights.m;i++) sum(face)+=local_weights(i).y;}

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(ghost_box.Lazy_Inside(iterator.Location())) continue;
        if(clamp_weights) if(sum(face)>=1) continue;
        if(!inside1(iterator.Full_Index())) continue; //Don't get velocity if not inside lsv
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        RANGE<TV_INT> inside_domain=grid.Domain_Indices();inside_domain.min_corner+=TV_INT::Axis_Vector(face.axis);
        T remaining=(1-sum(face));
        if(Is_Outside(inside_domain,face)) continue;
        TV X=iterator.Location()+dt*averaging.Face_To_Face_Vector(grid,face.axis,face.index,face_velocities);
        if(use_second_order) X+=(iterator.Location()-(X-dt*linear.Clamped_To_Array_Face(grid,face_velocities,X)))/2.;
        Clamp_X_To_Grid(grid,X);
        if(use_collision_object && collision_object.Lazy_Inside(X)) X=collision_object.Surface(X);
        if(!ghost_box.Lazy_Inside(iterator.Location()) && ghost_box.Lazy_Inside(X)) X=ghost_box.Surface(X);
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > forward_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X);
        for(int i=0;i<forward_weights.m;i++) if(!inside2(forward_weights(i).x)) forward_weights(i).y=0;
        T sum=0;for(int i=0;i<forward_weights.m;i++) sum+=forward_weights(i).y;
        assert(sum<=1+delta);
        if(sum<1 && sum>=delta) for(int i=0;i<forward_weights.m;i++) forward_weights(i).y=forward_weights(i).y/sum;
        else if(sum<delta){int iterations=100,iter=0;
            for(int i=0;i<forward_weights.m;i++) forward_weights(i).y=0;
            TV X_start=X;
            if(!(lsv2.Phi(X)>num_cells*grid.dX(face.axis) && !find_closest_point && pls && pls->Fix_Momentum_With_Escaped_Particles(X,face_velocities.Starting_Point_Face(face.axis,face.index).V_face,density*grid.dX.Product()*remaining,1.5,1,time,false))){
                while(lsv2.Phi(X)>0 && iter<iterations){iter++;X=X-(T)(lsv2.Phi(X))*lsv2.Normal(X).Normalized();}
                forward_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X);
                for(int i=0;i<forward_weights.m;i++) if(!inside2(forward_weights(i).x)) forward_weights(i).y=0;
                sum=0;for(int i=0;i<forward_weights.m;i++) sum+=forward_weights(i).y;
                assert(sum<=1+delta);
                if(sum<1 && sum>=delta) for(int i=0;i<forward_weights.m;i++) forward_weights(i).y=forward_weights(i).y/sum;
                else if(sum<delta){X=X_start;iter=0;
                    while(lsv2.Phi(X)>0 && iter<iterations){iter++;X=X+(T)(lsv2.Phi(X))*lsv2.Normal(X).Normalized();}
                    forward_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X);
                    for(int i=0;i<forward_weights.m;i++) if(!inside2(forward_weights(i).x)) forward_weights(i).y=0;
                    sum=0;for(int i=0;i<forward_weights.m;i++) sum+=forward_weights(i).y;
                    assert(sum<=1+delta);
                    if(sum<1 && sum>=delta) for(int i=0;i<forward_weights.m;i++) forward_weights(i).y=forward_weights(i).y/sum;
                    if(sum<delta){
                        weights_to(face).Append(PAIR<FACE_INDEX<TV::dimension>,T>(FACE_INDEX<TV::dimension>(),remaining));
                        PHYSBAM_WARNING("Cannot find LSV or Particle for Momentum");
                        for(int i=0;i<forward_weights.m;i++) forward_weights(i).y=0;
                        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > forward_weights_local=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X_start);
                        for(int i=0;i<forward_weights_local.m;i++){
                            momentum_lost(forward_weights_local(i).x.index)+=(T)(density*grid.dX.Product()*forward_weights_local(i).y*remaining/2.);
                            momentum_lost(forward_weights_local(i).x.index+TV_INT::Axis_Vector(forward_weights_local(i).x.axis))+=(T)(density*grid.dX.Product()*forward_weights_local(i).y*remaining/2.);}}}}
            else weights_to(face).Append(PAIR<FACE_INDEX<TV::dimension>,T>(FACE_INDEX<TV::dimension>(),remaining));}
        Clamp_Weights_To_Grid(inside_domain,forward_weights);
        Clamp_Weights_To_Objects(grid,forward_weights);
        for(int i=0;i<forward_weights.m;i++){//assert(forward_weights(i).y>-delta);
            if(forward_weights(i).y<0) forward_weights(i).y=0;}
        for(int i=0;i<forward_weights.m;i++){if(forward_weights(i).y==0) continue;
            int index=0;for(int j=0;j<weights_to(face).m;j++) if(weights_to(face)(j).x==forward_weights(i).x) index=j;
            if(index) weights_to(face)(index).y+=forward_weights(i).y*remaining;
            else{
                weights_to(face).Append(PAIR<FACE_INDEX<TV::dimension>,T>(forward_weights(i).x,forward_weights(i).y*remaining));
                weights_from(forward_weights(i).x).Append(PAIR<FACE_INDEX<TV::dimension>,int>(face,weights_to(face).m));}}}
    if(mpi_grid) mpi_grid->Sync_Common_Face_Weights_From(weights_to,weights_from,ghost_cells);

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(!inside1(iterator.Full_Index())) continue; //Don't clamp weights for outside cells
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
        T sum=0;for(int i=0;i<local_weights.m;i++) sum+=local_weights(i).y;
        for(int i=0;i<local_weights.m;i++){assert(sum>delta);
            local_weights(i).y=local_weights(i).y/sum;}}
    if(mpi_grid){
        mpi_grid->ignore_boundary_faces=true;
        mpi_grid->Sync_Common_Face_Weights_From(weights_to,weights_from,ghost_cells);
        mpi_grid->ignore_boundary_faces=false;}

    ARRAY<T,FACE_INDEX<TV::dimension> > sum_ic(grid,number_of_ghost_cells,false);
    boundary_sum->Fill_Ghost_Cells_Face(grid,sum_jc,sum_ic,time,number_of_ghost_cells);
    for(int iter=0;iter<num_iterations;iter++){
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
            if(ghost_box.Lazy_Inside(iterator.Location())) continue;
            if(!inside2(face)) continue; //Don't clamp weights for outside cells
            ARRAY<PAIR<FACE_INDEX<TV::dimension>,int> >& local_weights=weights_from(face);
            T sum=0;for(int i=0;i<local_weights.m;i++) sum+=sum_ic(local_weights(i).x)*weights_to(local_weights(i).x)(local_weights(i).y).y;
            for(int i=0;i<local_weights.m;i++){assert(sum>delta);
                weights_to(local_weights(i).x)(local_weights(i).y).y=weights_to(local_weights(i).x)(local_weights(i).y).y/sum;}}
        if(mpi_grid){
            mpi_grid->ignore_boundary_faces=true;
            mpi_grid->Sync_Common_Face_Weights_To(weights_to,weights_from,ghost_cells);
            mpi_grid->ignore_boundary_faces=false;}
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
            if(ghost_box.Lazy_Inside(iterator.Location())) continue;
            if(!inside1(face)) continue; //Don't clamp weights for outside cells
            ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
            T sum=0;for(int i=0;i<local_weights.m;i++) sum+=local_weights(i).y;
            for(int i=0;i<local_weights.m;i++){assert(sum>delta);
                local_weights(i).y=local_weights(i).y/sum;}}
        if(mpi_grid){
            mpi_grid->ignore_boundary_faces=true;
            mpi_grid->Sync_Common_Face_Weights_From(weights_to,weights_from,ghost_cells);
            mpi_grid->ignore_boundary_faces=false;}}
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(ghost_box.Lazy_Inside(iterator.Location())){sum_jc(face)=1;continue;}
        if(!inside2(face)){sum_jc(face)=1;continue;} //Don't clamp weights for outside cells
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,int> >& local_weights=weights_from(face);
        sum_jc(face)=0;for(int i=0;i<local_weights.m;i++) sum_jc(face)+=sum_ic(local_weights(i).x)*weights_to(local_weights(i).x)(local_weights(i).y).y;}


    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){if(ghost_box.Lazy_Inside(iterator.Location()) || !inside2(iterator.Full_Index())) continue;Z(iterator.Full_Index())=T();}
    for(FACE_ITERATOR iterator(grid,cfl);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(!inside1(iterator.Full_Index())) continue; //Don't get velocity if not inside lsv
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        if(Is_Outside(domain,face)) continue;
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
        if(!domain.Lazy_Inside(face.index) || ghost_box.Lazy_Inside(grid.Face(face.axis,face.index))){
            if(!ghost_box.Lazy_Inside(grid.Face(face.axis,face.index))){T sum=0;for(int i=0;i<local_weights.m;i++) sum+=local_weights(i).y;total_mass_gained+=sum*Z_ghost(face);}}
        for(int i=0;i<local_weights.m;i++)
            if(domain.Lazy_Inside(local_weights(i).x.index) && !ghost_box.Lazy_Inside(grid.Face(local_weights(i).x.axis,local_weights(i).x.index))) 
                Z(local_weights(i).x)+=local_weights(i).y*Z_ghost(face); 
            else{
                total_mass_lost+=local_weights(i).y*Z_ghost(face);}}

    //diffusion
    Face_Diffusion(grid,sum_jc,Z,boundary,boundary_sum,&inside2);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
    const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
    const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
{
    static bool first=true;
    weights_to.Clean_Memory();weights_from.Clean_Memory();
    PHYSBAM_ASSERT(!Z_min || !Z_max); //we don't support extrema yet
    BOUNDARY_UNIFORM<T_GRID,T> boundary_scalar;boundary_scalar.Set_Fixed_Boundary(true,1);
    BOUNDARY_UNIFORM<T_GRID,T>* boundary_sum=&boundary_scalar;
    if(mpi_grid) boundary_sum=new BOUNDARY_MPI<GRID<TV> >(mpi_grid,boundary_scalar);
    ARRAY<T,FACE_INDEX<TV::dimension> > sum;
    weights_to.Resize(grid,number_of_ghost_cells),weights_from.Resize(grid,number_of_ghost_cells),sum.Resize(grid),sum_jc.Resize(grid);
    if(first && time<=dt) for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) sum_jc(iterator.Full_Index())=1;
    first=false;
    T_INTERPOLATION interpolation;T_AVERAGING averaging;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> linear;
    int ghost_cells=number_of_ghost_cells;
    RANGE<TV> domain_x=grid.domain;domain_x.min_corner-=grid.dX*1.5;domain_x.max_corner+=grid.dX*1.5;
    RANGE<TV> real_domain=grid.domain;real_domain.min_corner-=grid.dX*((T)(1+1e-5));real_domain.max_corner+=grid.dX*((T)(1+1e-5));
    for(FACE_ITERATOR iterator(grid,cfl);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(ghost_box.Lazy_Inside(iterator.Location())) continue;
        RANGE<TV_INT> ghost_domain=grid.Domain_Indices(number_of_ghost_cells);ghost_domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        RANGE<TV_INT> inside_domain=grid.Domain_Indices();inside_domain.min_corner+=TV_INT::Axis_Vector(face.axis);
        if(Is_Outside(inside_domain,face)) continue;
        TV X=iterator.Location()-dt*averaging.Face_To_Face_Vector(grid,face.axis,face.index,face_velocities);
        if(use_second_order) X+=(iterator.Location()-(X+dt*linear.Clamped_To_Array_Face(grid,face_velocities,X)))/2.;
        Clamp_X_To_Grid(grid,X);
        if(use_collision_object && collision_object.Lazy_Inside(X)) X=collision_object.Surface(X);
        if(ghost_box.Lazy_Inside(X)) X=ghost_box.Surface(X);
        if(!domain.Lazy_Inside(face.index) && !real_domain.Lazy_Inside(X)) continue;
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > backwards_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X);
        Clean_Weights(backwards_weights);
        Clamp_Weights_To_Grid(inside_domain,backwards_weights);
        Clamp_Weights_To_Objects(grid,backwards_weights);
        for(int i=0;i<backwards_weights.m;i++){if(!ghost_domain.Lazy_Inside(backwards_weights(i).x.index)){assert(!domain.Lazy_Inside(face.index));continue;}
            if(backwards_weights(i).y==0) continue;
            weights_to(backwards_weights(i).x).Append(PAIR<FACE_INDEX<TV::dimension>,T>(face,backwards_weights(i).y));
            weights_from(face).Append(PAIR<FACE_INDEX<TV::dimension>,int>(backwards_weights(i).x,weights_to(backwards_weights(i).x).m));}}
    if(mpi_grid) mpi_grid->Sync_Common_Face_Weights_To(weights_to,weights_from,ghost_cells);

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
        sum(face)=0;for(int i=0;i<local_weights.m;i++) sum(face)+=local_weights(i).y;}

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(clamp_weights) if(sum(face)>=1) continue;
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        RANGE<TV_INT> inside_domain=grid.Domain_Indices();inside_domain.min_corner+=TV_INT::Axis_Vector(face.axis);
        T remaining=(1-sum(face));
        if(Is_Outside(inside_domain,face)) continue;
        TV X=iterator.Location()+dt*averaging.Face_To_Face_Vector(grid,face.axis,face.index,face_velocities);
        if(use_second_order) X+=(iterator.Location()-(X-dt*linear.Clamped_To_Array_Face(grid,face_velocities,X)))/2.;
        Clamp_X_To_Grid(grid,X);
        if(use_collision_object && collision_object.Lazy_Inside(X)) X=collision_object.Surface(X);
        if(!ghost_box.Lazy_Inside(iterator.Location()) && ghost_box.Lazy_Inside(X)) X=ghost_box.Surface(X);
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > forward_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X);
        Clean_Weights(forward_weights);
        Clamp_Weights_To_Grid(inside_domain,forward_weights);
        Clamp_Weights_To_Objects(grid,forward_weights);
        for(int i=0;i<forward_weights.m;i++){
            int index=0;for(int j=0;j<weights_to(face).m;j++) if(weights_to(face)(j).x==forward_weights(i).x) index=j;
            if(index) weights_to(face)(index).y+=forward_weights(i).y*remaining;
            else{
                if(forward_weights(i).y==0) continue;
                weights_to(face).Append(PAIR<FACE_INDEX<TV::dimension>,T>(forward_weights(i).x,forward_weights(i).y*remaining));
                weights_from(forward_weights(i).x).Append(PAIR<FACE_INDEX<TV::dimension>,int>(face,weights_to(face).m));}}}
    if(mpi_grid) mpi_grid->Sync_Common_Face_Weights_From(weights_to,weights_from,ghost_cells);

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
        T sum=0;for(int i=0;i<local_weights.m;i++) sum+=local_weights(i).y;
        for(int i=0;i<local_weights.m;i++){assert(sum>1e-5);
            local_weights(i).y=local_weights(i).y/sum;}}

    if(mpi_grid){
        mpi_grid->ignore_boundary_faces=true;
        mpi_grid->Sync_Common_Face_Weights_From(weights_to,weights_from,ghost_cells);
        mpi_grid->ignore_boundary_faces=false;}

    ARRAY<T,FACE_INDEX<TV::dimension> > sum_ic(grid,number_of_ghost_cells,false);
    boundary_sum->Fill_Ghost_Cells_Face(grid,sum_jc,sum_ic,time,number_of_ghost_cells);
    for(int iter=0;iter<num_iterations;iter++){
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
            RANGE<TV_INT> inside_domain=grid.Domain_Indices();inside_domain.min_corner+=TV_INT::Axis_Vector(face.axis);
            if(Is_Outside(inside_domain,face)) continue;
            if(ghost_box.Lazy_Inside(iterator.Location())) continue;
            ARRAY<PAIR<FACE_INDEX<TV::dimension>,int> >& local_weights=weights_from(face);
            T sum=0;for(int i=0;i<local_weights.m;i++) sum+=sum_ic(local_weights(i).x)*weights_to(local_weights(i).x)(local_weights(i).y).y;
            for(int i=0;i<local_weights.m;i++){assert(sum>1e-5);
                weights_to(local_weights(i).x)(local_weights(i).y).y=weights_to(local_weights(i).x)(local_weights(i).y).y/sum;}}
        if(mpi_grid){
            mpi_grid->ignore_boundary_faces=true;
            mpi_grid->Sync_Common_Face_Weights_To(weights_to,weights_from,ghost_cells);
            mpi_grid->ignore_boundary_faces=false;}
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
            RANGE<TV_INT> inside_domain=grid.Domain_Indices();inside_domain.min_corner+=TV_INT::Axis_Vector(face.axis);
            if(Is_Outside(inside_domain,face)) continue;
            if(ghost_box.Lazy_Inside(iterator.Location())) continue;
            ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
            T sum=0;for(int i=0;i<local_weights.m;i++) sum+=local_weights(i).y;
            for(int i=0;i<local_weights.m;i++){assert(sum>1e-5);
                local_weights(i).y=local_weights(i).y/sum;}}
        if(mpi_grid){
            mpi_grid->ignore_boundary_faces=true;
            mpi_grid->Sync_Common_Face_Weights_From(weights_to,weights_from,ghost_cells);
            mpi_grid->ignore_boundary_faces=false;}}

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        RANGE<TV_INT> inside_domain=grid.Domain_Indices();inside_domain.min_corner+=TV_INT::Axis_Vector(face.axis);
        if(Is_Outside(inside_domain,face)) continue;
        if(ghost_box.Lazy_Inside(iterator.Location())){sum_jc(face)=1;continue;}
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,int> >& local_weights=weights_from(face);
        sum_jc(face)=0;for(int i=0;i<local_weights.m;i++) sum_jc(face)+=sum_ic(local_weights(i).x)*weights_to(local_weights(i).x)(local_weights(i).y).y;}

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){if(ghost_box.Lazy_Inside(iterator.Location())) continue;Z(iterator.Full_Index())=T();}
    for(FACE_ITERATOR iterator(grid,ghost_cells);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        RANGE<TV_INT> inside_domain=grid.Domain_Indices();inside_domain.min_corner+=TV_INT::Axis_Vector(face.axis);
        if(Is_Outside(inside_domain,face)) continue;
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
        if((!domain.Lazy_Inside(face.index) && !Is_MPI_Boundary(domain,face)) || ghost_box.Lazy_Inside(grid.Face(face.axis,face.index))){
            T sum=0;for(int i=0;i<local_weights.m;i++) sum+=local_weights(i).y;total_mass_gained+=sum*Z_ghost(face);}
        for(int i=0;i<local_weights.m;i++)
            if(domain.Lazy_Inside(local_weights(i).x.index) && !ghost_box.Lazy_Inside(grid.Face(local_weights(i).x.axis,local_weights(i).x.index))){
                TV_INT second_face=local_weights(i).x.index+TV_INT::Axis_Vector(local_weights(i).x.axis);
                if(smoke_density && ((*smoke_density)(iterator.First_Cell_Index())<1e-5 && (*smoke_density)(iterator.Second_Cell_Index())<1e-5)) total_mass_gained+=local_weights(i).y*Z_ghost(face);
                if(smoke_density && ((*smoke_density)(local_weights(i).x.index)<1e-5 && (*smoke_density)(second_face)<1e-5)) total_mass_lost+=local_weights(i).y*Z_ghost(face);
                Z(local_weights(i).x)+=local_weights(i).y*Z_ghost(face);}
            else if(!Is_MPI_Boundary(domain,local_weights(i).x)) total_mass_lost+=local_weights(i).y*Z_ghost(face);}

    //diffusion
    Face_Diffusion(grid,sum_jc,Z,boundary,boundary_sum);
}

template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,1> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,2> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,3> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,1> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >,CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,1> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,2> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >,CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,2> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,1> >,VECTOR<float,3>,AVERAGING_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,VECTOR<float,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,4>,AVERAGING_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,4>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,3> >,VECTOR<float,5>,AVERAGING_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,VECTOR<float,5>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,1> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,2> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,3> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,1> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >,CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,1> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,2> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >,CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,2> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,1> >,VECTOR<double,3>,AVERAGING_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,VECTOR<double,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,4>,AVERAGING_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,4>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,3> >,VECTOR<double,5>,AVERAGING_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,VECTOR<double,5>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > > >;
#endif


