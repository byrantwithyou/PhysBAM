//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
namespace PhysBAM{
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> EIGENSYSTEM<T,d>::
~EIGENSYSTEM()
{
}
//#####################################################################
// Function All_Eigenvalues_Same
//#####################################################################
template<class T,int d> bool EIGENSYSTEM<T,d>::
All_Eigenvalues_Same()
{
    return false;
}
//#####################################################################
// Function Flux_Divided_By_Velocity
//#####################################################################
template<class T,int d> void EIGENSYSTEM<T,d>::
Flux_Divided_By_Velocity(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,
    ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Get_Face_Velocity_Component
//#####################################################################
template<class T,int d> T EIGENSYSTEM<T,d>::
Get_Face_Velocity_Component(const int face_index,const bool use_standard_average,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U)
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Flux_Using_Face_Velocity
//#####################################################################
template<class T,int d> void EIGENSYSTEM<T,d>::
Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,
    ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped)
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
//#####################################################################
// Function Maximum_Magnitude_Eigenvalue
//#####################################################################
template<class T,int d> T EIGENSYSTEM<T,d>::
Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell)
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template class EIGENSYSTEM<float,1>;
template class EIGENSYSTEM<float,2>;
template class EIGENSYSTEM<float,3>;
template class EIGENSYSTEM<float,4>;
template class EIGENSYSTEM<float,5>;
template class EIGENSYSTEM<float,6>;
template class EIGENSYSTEM<double,1>;
template class EIGENSYSTEM<double,2>;
template class EIGENSYSTEM<double,3>;
template class EIGENSYSTEM<double,4>;
template class EIGENSYSTEM<double,5>;
template class EIGENSYSTEM<double,6>;
}

