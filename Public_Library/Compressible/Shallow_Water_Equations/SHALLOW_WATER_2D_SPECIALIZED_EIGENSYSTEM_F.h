//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Eran Guendelman, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F  
//##################################################################### 
#ifndef __SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F__
#define __SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F__   

#include <Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
namespace PhysBAM{

template<class T>
class SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F:public EIGENSYSTEM<T,2>
{
    enum WORKAROUND1 {d=2};
public:
    using EIGENSYSTEM<T,2>::slice_index;

    T gravity;
    ARRAY<T,VECTOR<int,2> >& eta_ghost;
    T min_height;

    SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F(const T gravity_input,ARRAY<T,VECTOR<int,2> >& eta_ghost_input,const T min_height_input)
        :gravity(gravity_input),eta_ghost(eta_ghost_input),min_height(min_height_input)
    {}
    
//#####################################################################
    void Flux(const int m,const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& F,ARRAY<VECTOR<T,2> ,VECTOR<int,1> >* U_clamped=0) override;        
    bool Eigenvalues(const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,const int i,VECTOR<T,d>& lambda,VECTOR<T,d>& lambda_left,VECTOR<T,d>& lambda_right) override; 
    void Eigenvectors(const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) override; 
//#####################################################################
};
}    
#endif

