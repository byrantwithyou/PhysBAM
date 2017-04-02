//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_EIGENSYSTEM
//#####################################################################
#ifndef __SHALLOW_WATER_EIGENSYSTEM__
#define __SHALLOW_WATER_EIGENSYSTEM__

#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
namespace PhysBAM{

template<class TV>
class SHALLOW_WATER_EIGENSYSTEM:public EIGENSYSTEM<typename TV::SCALAR,TV::m+1>
{
    enum WORKAROUND1 {d=TV::m+1};typedef typename TV::SCALAR T;typedef VECTOR<T,d> TV_DIMENSION;
public:
    T gravity;
    int a;

    SHALLOW_WATER_EIGENSYSTEM(const T gravity_input,int a)
        :gravity(gravity_input),a(a)
    {}

//#####################################################################
    void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) override;
    T Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell) override;
    bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,VECTOR<T,d>& lambda,VECTOR<T,d>& lambda_left,VECTOR<T,d>& lambda_right) override;
    void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) override;
//#####################################################################
};
}
#endif
