//#####################################################################
// Copyright 2003-2007, Eran Guendelman, Nipun Kwatra, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_1D_EIGENSYSTEM_F  
//##################################################################### 
#ifndef __SHALLOW_WATER_1D_EIGENSYSTEM_F__
#define __SHALLOW_WATER_1D_EIGENSYSTEM_F__   

#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
namespace PhysBAM{

template<class T>
class SHALLOW_WATER_1D_EIGENSYSTEM_F:public EIGENSYSTEM<T,VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV_DIMENSION;
    enum WORKAROUND1 {d=TV_DIMENSION::m};
public:
    T gravity;

    SHALLOW_WATER_1D_EIGENSYSTEM_F(const T gravity_input)
        :gravity(gravity_input)
    {}
    
//#####################################################################
    void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) override;        
    T Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell) override;
    bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right) override; 
    void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) override; 
//#####################################################################
};
}    
#endif

