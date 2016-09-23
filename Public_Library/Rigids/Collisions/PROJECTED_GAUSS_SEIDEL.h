//#####################################################################
// Copyright 2010, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTED_GAUSS_SEIDEL
//##################################################################### 
#ifndef __PROJECTED_GAUSS_SEIDEL__
#define __PROJECTED_GAUSS_SEIDEL__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Log/SCOPE.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Vectors/VECTOR.h>
#include <Rigids/Collisions/SOLVE_CONTACT.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>

namespace PhysBAM{;

namespace PROJECTED_GAUSS_SEIDEL
{
template<class T> void Solve(SPARSE_MATRIX_FLAT_MXN<T>& A,ARRAY<T>& a,ARRAY<T>& x,T tolerance);

template<class T,int D> void Multiply(SPARSE_MATRIX_FLAT_MXN<VECTOR<T,D> >& A,ARRAY<VECTOR<T,D> >& x,ARRAY<T>& result);

template<class T,int D> void Solve(ARRAY<MATRIX<T,D> >& A_block_diagonal,
    ARRAY<ARRAY<PAIR<int,VECTOR<T,D> > > >& C_block,ARRAY<VECTOR<T,D> >& a_block,ARRAY<T>& c,
    ARRAY<VECTOR<T,D> >& x,ARRAY<T>& lambda,T tolerance);
template<class TV>
bool Solve(ARRAY<TWIST<TV> >& velocities,ARRAY<bool>& has_infinite_inertia,ARRAY<SOLVE_CONTACT::CONTACT<TV> >& contacts,
    ARRAY<typename TV::SCALAR>& lambda_normal,ARRAY<VECTOR<typename TV::SCALAR,TV::dimension-1> >& lambda_tangent,
    typename TV::SCALAR tolerance,int iteration_maximum,bool friction=false);
template<class TV> bool
Solve(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,ARRAY<SOLVE_CONTACT::CONTACT<TV> >& contacts,
    typename TV::SCALAR tolerance,int iteration_maximum);
}
}

#endif
