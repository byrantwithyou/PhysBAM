//#####################################################################
// Copyright 2007-2009, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY
//#####################################################################
#ifndef __EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY__
#define __EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Math_Tools/sqr.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
namespace PhysBAM{

template<class T_input>
class EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY:public EULER_EIGENSYSTEM<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,5> TV_DIMENSION;
    enum WORKAROUND1 {d=TV_DIMENSION::m};
public:
    EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY()
    {}

    bool All_Eigenvalues_Same() override {return true;}

//#####################################################################
    void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& H,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) override;
    void Flux_Divided_By_Velocity(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& H,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) override;
    T Get_Face_Velocity_Component(const int face_index,const bool use_standard_average,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U) override;
    void Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& H,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) override;
    T Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell) override;
    bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right) override;
    void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) override {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif
