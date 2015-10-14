//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_2D_EIGENSYSTEM_F  
//##################################################################### 
#ifndef __EULER_2D_EIGENSYSTEM_F__
#define __EULER_2D_EIGENSYSTEM_F__   

#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Matrices/MATRIX.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
namespace PhysBAM{

template<class T_input>
class EULER_2D_EIGENSYSTEM_F:public EULER_EIGENSYSTEM<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,4> TV_DIMENSION;
    enum WORKAROUND1 {d=TV_DIMENSION::m};
    typedef EULER_EIGENSYSTEM<VECTOR<T_input,2> > BASE;
public:
    using BASE::eos;using BASE::e;

    bool only_pressure_flux;

    EULER_2D_EIGENSYSTEM_F(bool only_pressure_flux_input=false):only_pressure_flux(only_pressure_flux_input)
    {}

//#####################################################################
    void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) override;
    void Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) override;
    T Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell) override;
    bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right) override; 
    void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) override; 
//#####################################################################
}; 
}  
#endif

