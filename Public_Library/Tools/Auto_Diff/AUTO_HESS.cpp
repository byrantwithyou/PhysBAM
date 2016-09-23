//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Auto_Diff/AUTO_HESS.h>
#include <Tools/Polynomials/QUADRATIC.h>
namespace PhysBAM{
template<class T> VECTOR<T,2>
Principal_Curvatures(const AUTO_HESS<T,VECTOR<T,3> >& h)
{
    SYMMETRIC_MATRIX<T,3> P=(T)1-SYMMETRIC_MATRIX<T,3>::Outer_Product(h.dx),M=SYMMETRIC_MATRIX<T,3>::Conjugate(P,h.ddx);
    T trace=M.Trace();
    QUADRATIC<T> quadratic(-1,trace,sqr(M(0,2))+sqr(M(0,1))+sqr(M(1,2))-M(1,1)*M(0,0)-M(2,2)*M(0,0)-M(2,2)*M(1,1));
    quadratic.Compute_Roots();
    if(quadratic.roots == 0) (T).5*VECTOR<T,2>(trace,trace);
    else if(quadratic.roots == 1) return VECTOR<T,2>(quadratic.root1,quadratic.root1);
    return VECTOR<T,2>(quadratic.root1,quadratic.root2);
}
template VECTOR<double,2> Principal_Curvatures<double>(AUTO_HESS<double,VECTOR<double,3> > const&);
template VECTOR<float,2> Principal_Curvatures<float>(AUTO_HESS<float,VECTOR<float,3> > const&);
}
