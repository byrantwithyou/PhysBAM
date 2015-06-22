//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_STEP_PARTICLES
//#####################################################################
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/EULER_STEP_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Function Clamped_To_Array
//#####################################################################
template<class TV> static TV Clamped_To_Array(LINEAR_INTERPOLATION_UNIFORM<TV,TV>& interpolation,const GRID<TV>& grid,
    const ARRAY<TV,VECTOR<int,TV::m> >& U,const TV& X)
{
    return interpolation.Clamped_To_Array(grid,U,X);
}
//#####################################################################
// Function Euler_Step_Node
//#####################################################################
template<class TV> void EULER_STEP_PARTICLES<TV>::
Euler_Step_Node(ARRAY_VIEW<TV> X,const GRID<TV>& grid,const ARRAY<TV,TV_INT>& U,const T dt)
{
    LINEAR_INTERPOLATION_UNIFORM<TV,TV> interpolation;
    for(int k=0;k<X.Size();k++) X(k)+=dt*Clamped_To_Array(interpolation,grid,U,X(k));
}
//#####################################################################
// Function Euler_Step_Face
//#####################################################################
template<class TV> void EULER_STEP_PARTICLES<TV>::
Euler_Step_Face(ARRAY_VIEW<TV> X,const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt)
{
    FACE_LOOKUP_UNIFORM<TV> lookup(face_velocities);
    for(int k=0;k<X.Size();k++){
        typename GRID<TV>::BLOCK block(grid,X(k),1);X(k)+=dt*LINEAR_INTERPOLATION_MAC_HELPER<TV>::Interpolate_Face(block,lookup,X(k));}
}
//#####################################################################
// Function Second_Order_Runge_Kutta_Step
//#####################################################################
template<class TV> void EULER_STEP_PARTICLES<TV>::
Second_Order_Runge_Kutta_Step(ARRAY_VIEW<TV> X,const GRID<TV>& grid,const ARRAY<TV,TV_INT>& U,const T dt)
{
    LINEAR_INTERPOLATION_UNIFORM<TV,TV> interpolation;
    for(int k=0;k<X.Size();k++){
        TV velocity=Clamped_To_Array(interpolation,grid,U,X(k));
        TV X_new=X(k)+dt*velocity;
        X(k)+=(T).5*dt*(velocity+Clamped_To_Array(interpolation,grid,U,X_new));}
}
//#####################################################################

template class EULER_STEP_PARTICLES<VECTOR<float,2> >;
template class EULER_STEP_PARTICLES<VECTOR<float,3> >;
template class EULER_STEP_PARTICLES<VECTOR<double,2> >;
template class EULER_STEP_PARTICLES<VECTOR<double,3> >;
}
