//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_UNIFORM  
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Arrays/FLOOD_FILL.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Grid_PDE/Poisson/LAPLACE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LAPLACE_UNIFORM<TV>::
LAPLACE_UNIFORM(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input,const bool initialize_grid,const bool enforce_compatibility_input)
    :grid(grid_input),u(u_input),mpi_grid(0),psi_D_save_for_sph(0),psi_N_save_for_sph(0),enforce_compatibility(enforce_compatibility_input),solve_single_cell_neumann_regions(false),use_psi_R(false)
{
    if(initialize_grid) Initialize_Grid(grid);
    laplace_mpi=new LAPLACE_UNIFORM_MPI<TV>(*this);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LAPLACE_UNIFORM<TV>::
~LAPLACE_UNIFORM()
{
    delete laplace_mpi;
}
//#####################################################################
// Function Solve
//#####################################################################
template<class TV> void LAPLACE_UNIFORM<TV>::
Solve(const T time,const bool solution_regions_already_computed)
{
    if(!solution_regions_already_computed) Find_Solution_Regions();
    ARRAY<ARRAY<TV_INT> > matrix_index_to_cell_index_array(number_of_regions);
    T_ARRAYS_INT cell_index_to_matrix_index(grid.Domain_Indices(1));
    cell_index_to_matrix_index.Fill(-1);
    ARRAY<int,VECTOR<int,1> > filled_region_cell_count(-2,number_of_regions);
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > A_array(number_of_regions);ARRAY<ARRAY<T> > b_array(number_of_regions);
    for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()) filled_region_cell_count(filled_region_colors(iterator.Cell_Index()))++;
    for(int color=0;color<number_of_regions;color++) if(filled_region_touches_dirichlet(color)||solve_neumann_regions){
        matrix_index_to_cell_index_array(color).Resize(filled_region_cell_count(color));}
    filled_region_cell_count.Fill(0); // reusing this array in order to make the indirection arrays
    if(!mpi_grid) Compute_Matrix_Indices(grid.Domain_Indices(1),filled_region_cell_count,matrix_index_to_cell_index_array,cell_index_to_matrix_index);
    else laplace_mpi->Find_Matrix_Indices(filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);
    RANGE<TV_INT> domain=grid.Domain_Indices(1);
    Find_A(domain,A_array,b_array,filled_region_cell_count,cell_index_to_matrix_index);
    for(int color=0;color<number_of_regions;color++) if(filled_region_cell_count(color)>0 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){
        pcg.Enforce_Compatibility(!filled_region_touches_dirichlet(color)&&enforce_compatibility);
        Solve_Subregion(matrix_index_to_cell_index_array(color),A_array(color),b_array(color),color);}
    if(!solve_neumann_regions) for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){
        int filled_region_color=filled_region_colors(iterator.Cell_Index());if(filled_region_color>0 && !filled_region_touches_dirichlet(filled_region_color)) u(iterator.Cell_Index())=0;}
}
//#####################################################################
// Function Find_A
//#####################################################################
template<class TV> void LAPLACE_UNIFORM<TV>::
Find_A(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& A_array,ARRAY<ARRAY<T> >& b_array,const ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    ARRAY<ARRAY<int> > row_counts(A_array.m,false);
    for(int i=0;i<A_array.m;i++){
        row_counts(i).Resize(filled_region_cell_count(i),false,false);
        b_array(i).Resize(filled_region_cell_count(i));}
    // TODO: this should be rewritten in terms of faces cause this got really hacky with MPI
    for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        int color=filled_region_colors(cell_index);assert(color!=-1);
        if(color!=-2 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){
            int row_count=1;
            for(int axis=0;axis<TV::m;axis++){TV_INT offset;offset[axis]=1;
                if(((filled_region_colors.Valid_Index(cell_index-offset) && filled_region_colors(cell_index-offset)==color) ||
                    (grid.Domain_Indices().Lazy_Outside_Half_Open(cell_index-offset) && periodic_boundary[axis])) && !psi_N.Component(axis)(cell_index)) row_count++;
                if(((filled_region_colors.Valid_Index(cell_index+offset) && filled_region_colors(cell_index+offset)==color) ||
                    (grid.Domain_Indices().Lazy_Outside_Half_Open(cell_index+offset) && periodic_boundary[axis])) && !psi_N.Component(axis)(cell_index+offset)) row_count++;}
            row_counts(color)(cell_index_to_matrix_index(cell_index))=row_count;}}
    for(int i=0;i<A_array.m;i++) A_array(i).Set_Row_Lengths(row_counts(i));
    TV one_over_dx2=Inverse(grid.dX*grid.dX);
    T default_row_sum=-2*one_over_dx2.Sum_Abs(),r=0;
    TV_INT grid_counts=grid.counts;
    for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        int color=filled_region_colors(cell_index);
        if(color!=-2 && (filled_region_touches_dirichlet(color)||solve_neumann_regions)){
            int matrix_index=cell_index_to_matrix_index(cell_index);
            SPARSE_MATRIX_FLAT_MXN<T>& A=A_array(filled_region_colors(cell_index));ARRAY<T>& b=b_array(filled_region_colors(cell_index));b(matrix_index)=f(cell_index);
            T row_sum=default_row_sum;
            for(int axis=0;axis<TV::m;axis++){TV_INT offset;offset[axis]=1;
                if(filled_region_colors.Valid_Index(cell_index-offset)){
                    if(use_psi_R && (r=psi_R.Component(axis)(cell_index))) row_sum+=one_over_dx2[axis]*r;
                    else if(psi_N.Component(axis)(cell_index)) row_sum+=one_over_dx2[axis];
                    else if(grid.Domain_Indices().Lazy_Outside_Half_Open(cell_index-offset) && periodic_boundary[axis]){
                        TV_INT periodic_offset_cell=cell_index-offset;
                        int axis_periodic_cell=wrap(periodic_offset_cell[axis],grid_counts[axis]);
                        periodic_offset_cell[axis]=axis_periodic_cell;
                        A.Set_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),one_over_dx2[axis]);}
                    else if(psi_D(cell_index-offset)) b(matrix_index)-=one_over_dx2[axis]*u(cell_index-offset);
                    else{assert(filled_region_colors(cell_index-offset)==color);
                        A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index-offset),one_over_dx2[axis]);}}
                if(filled_region_colors.Valid_Index(cell_index+offset)){
                    if(use_psi_R && (r=psi_R.Component(axis)(cell_index+offset))) row_sum+=one_over_dx2[axis]*r;
                    else if(psi_N.Component(axis)(cell_index+offset)) row_sum+=one_over_dx2[axis];
                    else if(grid.Domain_Indices().Lazy_Outside_Half_Open(cell_index+offset) && periodic_boundary[axis]){
                        TV_INT periodic_offset_cell=cell_index+offset;
                        int axis_periodic_cell=wrap(periodic_offset_cell[axis],grid_counts[axis]);
                        periodic_offset_cell[axis]=axis_periodic_cell;
                        A.Set_Element(matrix_index,cell_index_to_matrix_index(periodic_offset_cell),one_over_dx2[axis]);}
                    else if(psi_D(cell_index+offset)) b(matrix_index)-=one_over_dx2[axis]*u(cell_index+offset);
                    else{assert(filled_region_colors(cell_index+offset)==color);
                        A.Set_Element(matrix_index,cell_index_to_matrix_index(cell_index+offset),one_over_dx2[axis]);}}}
            A.Set_Element(matrix_index,matrix_index,row_sum);}} // set diagonal and right hand side
}
//#####################################################################
// Function Solve_Subregion
//#####################################################################
template<class TV> void LAPLACE_UNIFORM<TV>::
Solve_Subregion(ARRAY<TV_INT>& matrix_index_to_cell_index,SPARSE_MATRIX_FLAT_MXN<T>& A,ARRAY<T>& b,const int color)
{
    int number_of_unknowns=matrix_index_to_cell_index.m;
    A.Negate();b*=(T)-1;
    ARRAY<T> x(number_of_unknowns),q;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;
    for(int i=0;i<number_of_unknowns;i++) x(i)=u(matrix_index_to_cell_index(i));
    Find_Tolerance(b); // needs to happen after b is completely set up
    if(pcg.show_results) LOG::cout << "solving " << number_of_unknowns << " cells to tolerance " << tolerance << std::endl;
    if(!mpi_grid) pcg.Solve(A,x,b,vectors,tolerance);
    else laplace_mpi->Solve(A,x,b,vectors,tolerance,color);
    for(int i=0;i<number_of_unknowns;i++){TV_INT cell_index=matrix_index_to_cell_index(i);u(cell_index)=x(i);}
    vectors.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Compute_Matrix_Indices
//#####################################################################
template<class TV> void LAPLACE_UNIFORM<TV>::
Compute_Matrix_Indices(ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    RANGE<TV_INT> domain(grid.Domain_Indices(1));
    Compute_Matrix_Indices(domain,filled_region_cell_count,matrix_index_to_cell_index_array,cell_index_to_matrix_index);
}
template<class TV> void LAPLACE_UNIFORM<TV>::
Compute_Matrix_Indices(const RANGE<TV_INT>& domain,ARRAY<int,VECTOR<int,1> >& filled_region_cell_count,ARRAY<ARRAY<TV_INT> >& matrix_index_to_cell_index_array,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    for(CELL_ITERATOR<TV> iterator(grid,domain);iterator.Valid();iterator.Next()){
        int color=filled_region_colors(iterator.Cell_Index());
        if(color>=0&&(filled_region_touches_dirichlet(color)||solve_neumann_regions)){
            matrix_index_to_cell_index_array(color)(filled_region_cell_count(color))=iterator.Cell_Index();
            cell_index_to_matrix_index(iterator.Cell_Index())=filled_region_cell_count(color)++;}}
}
//#####################################################################
// Function Build_Single_Solution_Region
//#####################################################################
template<class TV> void LAPLACE_UNIFORM<TV>::
Build_Single_Solution_Region(ARRAY<bool,TV_INT>& solve)
{
    number_of_regions=1;
    for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()) filled_region_colors(iterator.Cell_Index())=solve(iterator.Cell_Index())?0:-2;
    filled_region_touches_dirichlet.Resize(1);filled_region_touches_dirichlet(0)=true;
}
//#####################################################################
// Function Find_Solution_Regions
//#####################################################################
template<class TV> void LAPLACE_UNIFORM<TV>::
Find_Solution_Regions()
{
    FLOOD_FILL<TV::m> flood_fill;
    // set domain boundary cells and cells with objects to uncolorable
    for(CELL_ITERATOR<TV> iterator(grid,1,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next()) filled_region_colors(iterator.Cell_Index())=-2;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        if(psi_D(iterator.Cell_Index()) || (!solve_single_cell_neumann_regions && All_Cell_Faces_Neumann(iterator.Cell_Index()))) filled_region_colors(iterator.Cell_Index())=-2;
        else filled_region_colors(iterator.Cell_Index())=-1;}
    filled_region_touches_dirichlet.Remove_All();
    // do the fill
    if(mpi_grid){
        for(int axis=0;axis<TV::m;axis++) for(int side=0;side<2;side++) for(CELL_ITERATOR<TV> iterator(grid,1,GRID<TV>::GHOST_REGION,2*axis+side);iterator.Valid();iterator.Next()){
            for(int face=0;face<TV::m;face++)if(face!=axis){
                psi_N.Component(face)(iterator.Cell_Index())=true;psi_N.Component(face)(iterator.Cell_Index()+TV_INT::Axis_Vector(face))=true;}}
        for(int axis=0;axis<TV::m;axis++) for(int side=0;side<2;side++) for(CELL_ITERATOR<TV> iterator(grid,1,GRID<TV>::GHOST_REGION,2*axis+side);iterator.Valid();iterator.Next()){
            if(!psi_N.Component(axis)(iterator.Cell_Index()+(1-side)*TV_INT::Axis_Vector(axis))&&!psi_D(iterator.Cell_Index()))filled_region_colors(iterator.Cell_Index())=-1;}}
    number_of_regions=flood_fill.Flood_Fill(filled_region_colors,psi_N,&filled_region_touches_dirichlet);
    // correct flood fill for distributed grids
    if(mpi_grid)laplace_mpi->Synchronize_Solution_Regions();
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class TV> void LAPLACE_UNIFORM<TV>::
Initialize_Grid(const GRID<TV>& mac_grid_input)
{
    assert(mac_grid_input.DX()==TV() || mac_grid_input.Is_MAC_Grid());
    grid=mac_grid_input;
    f.Resize(grid.Domain_Indices(1));
    psi_N.Resize(grid,1);
    psi_D.Resize(grid.Domain_Indices(1));
    filled_region_colors.Resize(grid.Domain_Indices(1));
}
//#####################################################################
// Function Set_Neumann_Outer_Boundaries
//#####################################################################
template<class TV> void LAPLACE_UNIFORM<TV>::
Set_Neumann_Outer_Boundaries()
{
    for(FACE_ITERATOR<TV> iterator(grid,0,GRID<TV>::BOUNDARY_REGION);iterator.Valid();iterator.Next()) psi_N.Component(iterator.Axis())(iterator.Face_Index())=true;
    pcg.Enforce_Compatibility();
}
//#####################################################################
// Function Set_Dirichlet_Outer_Boundaries
//#####################################################################
template<class TV> void LAPLACE_UNIFORM<TV>::
Set_Dirichlet_Outer_Boundaries()
{
    for(CELL_ITERATOR<TV> iterator(grid,1,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next()) psi_D(iterator.Cell_Index())=true;
}
//#####################################################################
// Function Use_Psi_R
//#####################################################################
template<class TV> void LAPLACE_UNIFORM<TV>::
Use_Psi_R()
{
    use_psi_R=true;
    psi_R.Resize(psi_N.Domain_Indices());
}
//#####################################################################
namespace PhysBAM{
template class LAPLACE_UNIFORM<VECTOR<float,1> >;
template class LAPLACE_UNIFORM<VECTOR<float,2> >;
template class LAPLACE_UNIFORM<VECTOR<float,3> >;
template class LAPLACE_UNIFORM<VECTOR<double,1> >;
template class LAPLACE_UNIFORM<VECTOR<double,2> >;
template class LAPLACE_UNIFORM<VECTOR<double,3> >;
}
