//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_STEP_PARTICLES
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EULER_STEP_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Function Clamped_To_Array
//#####################################################################
template<class TV,class T_GRID> static TV Clamped_To_Array(LINEAR_INTERPOLATION_UNIFORM<T_GRID,TV>& interpolation,const T_GRID& grid,
    const ARRAY<TV,VECTOR<int,TV::m> >& U,const TV& X)
{
    return interpolation.Clamped_To_Array(grid,U,X);
}
//#####################################################################
// Function Euler_Step_Node
//#####################################################################
template<class T_GRID> void EULER_STEP_PARTICLES<T_GRID>::
Euler_Step_Node(ARRAY_VIEW<TV> X,const T_GRID& grid,const ARRAY<TV,TV_INT>& U,const T dt)
{
    typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR::template REBIND<TV>::TYPE interpolation;
    for(int k=0;k<X.Size();k++) X(k)+=dt*Clamped_To_Array(interpolation,grid,U,X(k));
}
//#####################################################################
// Function Euler_Step_Face
//#####################################################################
template<class T_GRID> void EULER_STEP_PARTICLES<T_GRID>::
Euler_Step_Face(ARRAY_VIEW<TV> X,const T_GRID& grid,const T_FACE_ARRAYS& face_velocities,const T dt)
{
    typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP lookup(face_velocities);
    for(int k=0;k<X.Size();k++){
        typename T_GRID::BLOCK block(grid,X(k),1);X(k)+=dt*INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face(block,lookup,X(k));}
}
//#####################################################################
// Function Second_Order_Runge_Kutta_Step
//#####################################################################
template<class T_GRID> void EULER_STEP_PARTICLES<T_GRID>::
Second_Order_Runge_Kutta_Step(ARRAY_VIEW<TV> X,const T_GRID& grid,const ARRAY<TV,TV_INT>& U,const T dt)
{
    typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR::template REBIND<TV>::TYPE interpolation;
    for(int k=0;k<X.Size();k++){
        TV velocity=Clamped_To_Array(interpolation,grid,U,X(k));
        TV X_new=X(k)+dt*velocity;
        X(k)+=(T).5*dt*(velocity+Clamped_To_Array(interpolation,grid,U,X_new));}
}
//#####################################################################

template class EULER_STEP_PARTICLES<GRID<VECTOR<float,2> > >;
template class EULER_STEP_PARTICLES<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EULER_STEP_PARTICLES<GRID<VECTOR<double,2> > >;
template class EULER_STEP_PARTICLES<GRID<VECTOR<double,3> > >;
#endif
}
