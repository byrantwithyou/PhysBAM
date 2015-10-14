//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_NONCONSERVATIVE_FLUX_FORM_ENO  
//##################################################################### 
#ifndef __ADVECTION_NONCONSERVATIVE_FLUX_FORM_ENO__
#define __ADVECTION_NONCONSERVATIVE_FLUX_FORM_ENO__

#include <Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input,class T2>
class ADVECTION_NONCONSERVATIVE_FLUX_FORM_ENO:public ADVECTION_SEPARABLE_UNIFORM<T_input,T2>
{
    typedef T_input T;

    int order;
public:

    ADVECTION_NONCONSERVATIVE_FLUX_FORM_ENO()
    {
        Set_Order();
    }

    void Set_Order(const int order_input=3)  override
    {order=order_input;assert(order >=1 && order <=3);}

//#####################################################################
    void Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx) override;
//#####################################################################
};  
} 
#endif

