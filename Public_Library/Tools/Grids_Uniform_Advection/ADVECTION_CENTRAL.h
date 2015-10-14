//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_CENTRAL
//##################################################################### 
#ifndef __ADVECTION_CENTRAL__
#define __ADVECTION_CENTRAL__

#include <Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
namespace PhysBAM{

template<class TV,class T2>
class ADVECTION_CENTRAL:public ADVECTION_SEPARABLE_UNIFORM<TV,T2>
{
    typedef typename TV::SCALAR T;
public:
    ADVECTION_CENTRAL()
    {}

//#####################################################################
    void Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx) override;
//#####################################################################
};
}    
#endif
