//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REACTIVE_EULER_EIGENSYSTEM  
//#################################################################### 
#ifndef __REACTIVE_EULER_EIGENSYSTEM__
#define __REACTIVE_EULER_EIGENSYSTEM__   

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
#include <Compressible/Reactive_Euler_Equations/REACTIVE_EULER.h>
namespace PhysBAM{

template<class TV>
class REACTIVE_EULER_EIGENSYSTEM:public EIGENSYSTEM<typename TV::SCALAR,TV::m+3>
{
    enum WORKAROUND1 {d=TV::m+3};
    typedef typename TV::SCALAR T;typedef VECTOR<T,d> TV_DIMENSION;
public:

    REACTIVE_EOS<T>* eos;
    int a;

    REACTIVE_EULER_EIGENSYSTEM(REACTIVE_EOS<T>& eos_input,int a)
        :eos(&eos_input),a(a)
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
