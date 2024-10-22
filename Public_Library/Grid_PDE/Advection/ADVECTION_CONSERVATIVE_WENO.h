//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_CONSERVATIVE_WENO  
//##################################################################### 
#ifndef __ADVECTION_CONSERVATIVE_WENO__
#define __ADVECTION_CONSERVATIVE_WENO__

#include <Grid_PDE/Advection/ADVECTION_SEPARABLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input,class T2>
class ADVECTION_CONSERVATIVE_WENO:public ADVECTION_SEPARABLE_UNIFORM<T_input,T2>
{
    typedef T_input T;

    bool compute_epsilon;
    double epsilon;
public:

    ADVECTION_CONSERVATIVE_WENO()
    {
        Compute_Epsilon();
        Set_Epsilon();
    }

    void Set_Epsilon(const T epsilon_input=1e-6) override
    {compute_epsilon=false;epsilon=epsilon_input;}
    
    void Compute_Epsilon() override
    {compute_epsilon=true;}

//#####################################################################
    void Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx) override;
//#####################################################################
};   
}
#endif

