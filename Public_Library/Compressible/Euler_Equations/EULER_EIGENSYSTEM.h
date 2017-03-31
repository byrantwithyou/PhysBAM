//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_EIGENSYSTEM  
//##################################################################### 
#ifndef __EULER_EIGENSYSTEM__
#define __EULER_EIGENSYSTEM__   

#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
namespace PhysBAM{

template<class T> class EOS;

template<class TV>
class EULER_EIGENSYSTEM:public EIGENSYSTEM<typename TV::SCALAR,VECTOR<typename TV::SCALAR,TV::m+2> >
{
    enum WORKAROUND1 {d=TV::m+2};
    typedef typename TV::SCALAR T;typedef VECTOR<T,d> TV_DIMENSION;
public:

    EOS<T>* eos;
    bool only_pressure_flux;
    bool only_advection;
    int a;

    EULER_EIGENSYSTEM(EOS<T>* eos_input,int a)
        :eos(eos_input),only_pressure_flux(false),only_advection(false),a(a)
    {}

//#####################################################################
    void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) override;
    void Flux_Divided_By_Velocity(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped) override;
    T Get_Face_Velocity_Component(const int face_index,const bool use_standard_average,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U) override;
    void Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) override;
    T Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell) override;
    bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,VECTOR<T,d>& lambda,VECTOR<T,d>& lambda_left,VECTOR<T,d>& lambda_right) override;
    void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) override;
    bool All_Eigenvalues_Same() override;
//#####################################################################
};
}
#endif

