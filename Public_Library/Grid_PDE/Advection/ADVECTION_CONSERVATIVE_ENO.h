//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_CONSERVATIVE_ENO  
//##################################################################### 
#ifndef __ADVECTION_CONSERVATIVE_ENO__
#define __ADVECTION_CONSERVATIVE_ENO__

#include <Grid_PDE/Advection/ADVECTION_SEPARABLE_UNIFORM.h>
namespace PhysBAM{

template<class TV,class T2>
class ADVECTION_CONSERVATIVE_ENO:public ADVECTION_SEPARABLE_UNIFORM<TV,T2>
{
    typedef typename TV::SCALAR T;

    int order;
public:
    using ADVECTION_SEPARABLE_UNIFORM<TV,T2>::Advection_Solver;using ADVECTION_SEPARABLE_UNIFORM<TV,T2>::ENO;

    ADVECTION_CONSERVATIVE_ENO();
    virtual ~ADVECTION_CONSERVATIVE_ENO();

//#####################################################################
    void Set_Order(const int order_input=3) override;
    void Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx) override;
//#####################################################################
};   
}
#endif
