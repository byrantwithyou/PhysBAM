#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MPI.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T2> BOUNDARY_MPI<TV,T2>::
BOUNDARY_MPI(T_MPI_GRID* mpi_grid_input,BOUNDARY<TV,T2>& boundary_input)
    :mpi_grid(mpi_grid_input),boundary(boundary_input)
{
    assert(&boundary);
}
//#####################################################################
// Function Set_Constant_Extrapolation
//#####################################################################
template<class TV,class T2> void BOUNDARY_MPI<TV,T2>::
Set_Constant_Extrapolation(const TV_SIDES& constant_extrapolation_input)
{
    boundary.Set_Constant_Extrapolation(constant_extrapolation_input);
}
//#####################################################################
// Function Constant_Extrapolation
//#####################################################################
template<class TV,class T2> bool BOUNDARY_MPI<TV,T2>::
Constant_Extrapolation(const int side) const
{
    return boundary.Constant_Extrapolation(side);
}
//#####################################################################
// Function Set_Fixed_Boundary
//#####################################################################
template<class TV,class T2> void BOUNDARY_MPI<TV,T2>::
Set_Fixed_Boundary(const bool use_fixed_boundary_input,const T2 fixed_boundary_value_input)
{
    boundary.Set_Fixed_Boundary(use_fixed_boundary_input,fixed_boundary_value_input);
}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY_MPI<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells_input) const
{
    boundary.Fill_Ghost_Cells(grid,u,u_ghost,dt,time,number_of_ghost_cells_input);
    mpi_grid->Exchange_Boundary_Cell_Data(u_ghost,number_of_ghost_cells_input);
}
//#####################################################################
// Function Fill_Ghost_Faces
//#####################################################################
template<class TV,class T2> void BOUNDARY_MPI<TV,T2>::
Fill_Ghost_Faces(const GRID<TV>& grid,const T_FACE_ARRAYS_T2& u,T_FACE_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells_input) const
{
    boundary.Fill_Ghost_Faces(grid,u,u_ghost,time,number_of_ghost_cells_input);
    mpi_grid->Exchange_Boundary_Face_Data(u_ghost,number_of_ghost_cells_input);
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV,class T2> void BOUNDARY_MPI<TV,T2>::
Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) const
{
    boundary.Apply_Boundary_Condition(grid,u,time);
}
//#####################################################################
// Function Apply_Boundary_Condition_Face
//#####################################################################
template<class TV,class T2> void BOUNDARY_MPI<TV,T2>::
Apply_Boundary_Condition_Face(const GRID<TV>& grid,T_FACE_ARRAYS_T2& u,const T time) const
{
    boundary.Apply_Boundary_Condition_Face(grid,u,time);
    mpi_grid->Average_Common_Face_Data(u);
}
namespace PhysBAM{
template class BOUNDARY_MPI<VECTOR<float,1>,VECTOR<float,1> >;
template class BOUNDARY_MPI<VECTOR<float,1>,VECTOR<float,3> >;
template class BOUNDARY_MPI<VECTOR<float,1>,float>;
template class BOUNDARY_MPI<VECTOR<float,2>,VECTOR<float,2> >;
template class BOUNDARY_MPI<VECTOR<float,2>,VECTOR<float,4> >;
template class BOUNDARY_MPI<VECTOR<float,2>,float>;
template class BOUNDARY_MPI<VECTOR<float,3>,VECTOR<float,3> >;
template class BOUNDARY_MPI<VECTOR<float,3>,VECTOR<float,5> >;
template class BOUNDARY_MPI<VECTOR<float,3>,float>;
template class BOUNDARY_MPI<VECTOR<float,1>,SYMMETRIC_MATRIX<float,1> >;
template class BOUNDARY_MPI<VECTOR<float,2>,SYMMETRIC_MATRIX<float,2> >;
template class BOUNDARY_MPI<VECTOR<float,3>,SYMMETRIC_MATRIX<float,3> >;
template class BOUNDARY_MPI<VECTOR<double,1>,VECTOR<double,1> >;
template class BOUNDARY_MPI<VECTOR<double,1>,VECTOR<double,3> >;
template class BOUNDARY_MPI<VECTOR<double,1>,double>;
template class BOUNDARY_MPI<VECTOR<double,2>,VECTOR<double,2> >;
template class BOUNDARY_MPI<VECTOR<double,2>,VECTOR<double,4> >;
template class BOUNDARY_MPI<VECTOR<double,2>,double>;
template class BOUNDARY_MPI<VECTOR<double,3>,VECTOR<double,3> >;
template class BOUNDARY_MPI<VECTOR<double,3>,VECTOR<double,5> >;
template class BOUNDARY_MPI<VECTOR<double,3>,double>;
template class BOUNDARY_MPI<VECTOR<double,1>,SYMMETRIC_MATRIX<double,1> >;
template class BOUNDARY_MPI<VECTOR<double,2>,SYMMETRIC_MATRIX<double,2> >;
template class BOUNDARY_MPI<VECTOR<double,3>,SYMMETRIC_MATRIX<double,3> >;
}
