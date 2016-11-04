//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/MATRIX.h>
#include <Grid_PDE/Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM_DEFINITION.h>
#include <Grid_PDE/Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
#include <Grid_PDE/Interpolation/CUBIC_MONOTONIC_INTERPOLATION_UNIFORM.h>
#include <Grid_PDE/Interpolation/CUBIC_SPLINE_INTERPOLATION_UNIFORM.h>
#include <Grid_PDE/Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>

namespace PhysBAM{
#define INSTANTIATION_HELPER(T,TV,d) \
    template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,SYMMETRIC_MATRIX<T,d>,AVERAGING_UNIFORM<TV,FACE_LOOKUP_UNIFORM<TV> >,LINEAR_INTERPOLATION_UNIFORM<TV,SYMMETRIC_MATRIX<T,d>,FACE_LOOKUP_UNIFORM<TV> > >; \
    template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,MATRIX<T,d>,AVERAGING_UNIFORM<TV,FACE_LOOKUP_UNIFORM<TV> >,LINEAR_INTERPOLATION_UNIFORM<TV,MATRIX<T,d>,FACE_LOOKUP_UNIFORM<TV> > >; \
    template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV,FACE_LOOKUP_UNIFORM<TV> >,LINEAR_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> > >; \
    template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV,FACE_LOOKUP_UNIFORM<TV> >,QUADRATIC_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> > >; \
    template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV,FACE_LOOKUP_UNIFORM<TV> >,CUBIC_MN_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> > >; \
    template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV,FACE_LOOKUP_UNIFORM<TV> >,CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> > >; \
    template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T,AVERAGING_UNIFORM<TV,FACE_LOOKUP_UNIFORM<TV> >,CUBIC_SPLINE_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> > >;
    
#define P(...) __VA_ARGS__
INSTANTIATION_HELPER(float,P(VECTOR<float,1>),2)
INSTANTIATION_HELPER(float,P(VECTOR<float,2>),2);
INSTANTIATION_HELPER(float,P(VECTOR<float,3>),3);
INSTANTIATION_HELPER(double,P(VECTOR<double,1>),2);
INSTANTIATION_HELPER(double,P(VECTOR<double,2>),2);
INSTANTIATION_HELPER(double,P(VECTOR<double,3>),3);
template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<VECTOR<float,2>,SYMMETRIC_MATRIX<float,2>,AVERAGING_UNIFORM<VECTOR<float,2>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > >,QUADRATIC_INTERPOLATION_UNIFORM<VECTOR<float,2>,SYMMETRIC_MATRIX<float,2>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > > >;
template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<VECTOR<double,2>,SYMMETRIC_MATRIX<double,2>,AVERAGING_UNIFORM<VECTOR<double,2>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > >,QUADRATIC_INTERPOLATION_UNIFORM<VECTOR<double,2>,SYMMETRIC_MATRIX<double,2>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > > >;
template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<VECTOR<float,3>,SYMMETRIC_MATRIX<float,3>,AVERAGING_UNIFORM<VECTOR<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >,QUADRATIC_INTERPOLATION_UNIFORM<VECTOR<float,3>,SYMMETRIC_MATRIX<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > > >;
template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<VECTOR<double,3>,SYMMETRIC_MATRIX<double,3>,AVERAGING_UNIFORM<VECTOR<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >,QUADRATIC_INTERPOLATION_UNIFORM<VECTOR<double,3>,SYMMETRIC_MATRIX<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > > >;
template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<VECTOR<float,2>,MATRIX<float,2>,AVERAGING_UNIFORM<VECTOR<float,2>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > >,QUADRATIC_INTERPOLATION_UNIFORM<VECTOR<float,2>,MATRIX<float,2>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > > >;
template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<VECTOR<double,2>,MATRIX<double,2>,AVERAGING_UNIFORM<VECTOR<double,2>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > >,QUADRATIC_INTERPOLATION_UNIFORM<VECTOR<double,2>,MATRIX<double,2>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > > >;
template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<VECTOR<float,3>,MATRIX<float,3>,AVERAGING_UNIFORM<VECTOR<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >,QUADRATIC_INTERPOLATION_UNIFORM<VECTOR<float,3>,MATRIX<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > > >;
template class ADVECTION_SEMI_LAGRANGIAN_UNIFORM<VECTOR<double,3>,MATRIX<double,3>,AVERAGING_UNIFORM<VECTOR<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >,QUADRATIC_INTERPOLATION_UNIFORM<VECTOR<double,3>,MATRIX<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > > >;
}