//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BURGERS_1D_EIGENSYSTEM_F  
//##################################################################### 
#ifndef __BURGERS_1D_EIGENSYSTEM_F__
#define __BURGERS_1D_EIGENSYSTEM_F__   

#include <Core/Matrices/MATRIX_FORWARD.h>
#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
namespace PhysBAM{

template<class T>
class BURGERS_1D_EIGENSYSTEM_F:public EIGENSYSTEM<T,1>
{
    typedef VECTOR<T,1> TV_DIMENSION;
    enum WORKAROUND1 {d=TV_DIMENSION::m};
    
public:
    BURGERS_1D_EIGENSYSTEM_F()
    {}
    
//#####################################################################
    void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) override;        
    bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,VECTOR<T,d>& lambda,VECTOR<T,d>& lambda_left,VECTOR<T,d>& lambda_right) override; 
    void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) override; 
//#####################################################################
};
}    
#endif
