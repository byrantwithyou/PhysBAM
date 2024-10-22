//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_HAMILTON_JACOBI_ENO  
//##################################################################### 
#ifndef __ADVECTION_HAMILTON_JACOBI_ENO__
#define __ADVECTION_HAMILTON_JACOBI_ENO__

#include <Grid_PDE/Advection/ADVECTION_SEPARABLE_UNIFORM.h>
namespace PhysBAM{

template<class TV,class T2>
class ADVECTION_HAMILTON_JACOBI_ENO:public ADVECTION_SEPARABLE_UNIFORM<TV,T2>
{
private:
    typedef typename TV::SCALAR T;

    int order;
public:
    using ADVECTION_SEPARABLE_UNIFORM<TV,T2>::ENO;

    ADVECTION_HAMILTON_JACOBI_ENO()
    {
        Set_Order();
    }

    void Set_Order(const int order_input=3)  override
    {order=order_input;assert(order >=1 && order <=3);}

//#####################################################################
    void Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx) override;
    void Advection_Solver(const int m_start,const int m_end,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx);
//#####################################################################
};
}    
#endif
