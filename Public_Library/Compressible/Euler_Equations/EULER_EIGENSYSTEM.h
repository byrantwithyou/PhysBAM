//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_EIGENSYSTEM  
//##################################################################### 
#ifndef __EULER_EIGENSYSTEM__
#define __EULER_EIGENSYSTEM__   

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Matrices/MATRIX.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM_BASE.h>
namespace PhysBAM{

template<class TV>
class EULER_EIGENSYSTEM:public EULER_EIGENSYSTEM_BASE<TV>
{
    enum WORKAROUND1 {d=TV::m+2};
    typedef typename TV::SCALAR T;typedef VECTOR<T,d> TV_DIMENSION;
    typedef EULER_EIGENSYSTEM_BASE<TV> BASE;
public:
    using BASE::eos;

    bool only_pressure_flux;
    int a;

    EULER_EIGENSYSTEM(EOS<T>* eos_input,int a)
        :BASE(eos_input),only_pressure_flux(false),a(a)
    {}

//#####################################################################
    void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) override;
    void Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) override;
    T Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell) override;
    bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,VECTOR<T,d>& lambda,VECTOR<T,d>& lambda_left,VECTOR<T,d>& lambda_right) override;
    void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) override;
//#####################################################################
};
}
#endif

