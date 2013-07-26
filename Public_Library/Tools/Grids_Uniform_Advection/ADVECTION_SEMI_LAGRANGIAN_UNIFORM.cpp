//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM_DEFINITION.h>
#include <Tools/Grids_Uniform_Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>

namespace PhysBAM{
#define INSTANTIATION_HELPER(T,TV,d) \
    template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,SYMMETRIC_MATRIX<T,d>,AVERAGING_UNIFORM<TV,FACE_LOOKUP_UNIFORM<TV> >,LINEAR_INTERPOLATION_UNIFORM<TV,SYMMETRIC_MATRIX<T,d>,FACE_LOOKUP_UNIFORM<TV> > >; \
    template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV,FACE_LOOKUP_UNIFORM<TV> >,LINEAR_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> > >; \
    template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV,FACE_LOOKUP_UNIFORM<TV> >,QUADRATIC_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> > >; \
    template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV,FACE_LOOKUP_UNIFORM<TV> >,CUBIC_MN_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> > >;
    
#define P(...) __VA_ARGS__
INSTANTIATION_HELPER(float,P(VECTOR<float,1>),2)
INSTANTIATION_HELPER(float,P(VECTOR<float,2>),2);
INSTANTIATION_HELPER(float,P(VECTOR<float,3>),3);
INSTANTIATION_HELPER(double,P(VECTOR<double,1>),2);
INSTANTIATION_HELPER(double,P(VECTOR<double,2>),2);
INSTANTIATION_HELPER(double,P(VECTOR<double,3>),3);
}
