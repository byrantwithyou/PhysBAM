//#####################################################################
// Copyright 2002-2004, Doug Enright, Ronald Fedkiw, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION_CENTRAL_ENO  
//##################################################################### 
#ifndef __CONSERVATION_CENTRAL_ENO__
#define __CONSERVATION_CENTRAL_ENO__   

#include <Compressible/Conservation_Law_Solvers/CONSERVATION.h>
namespace PhysBAM{

template<class TV,int d>
class CONSERVATION_CENTRAL_ENO:public CONSERVATION<TV,d>
{
    typedef typename TV::SCALAR T;typedef VECTOR<T,d> TV_DIMENSION;
public:
    typedef CONSERVATION<TV,d> BASE;
    using BASE::order;using BASE::Set_Order;using BASE::Use_Maximum_Alpha;using BASE::Alpha;

    CONSERVATION_CENTRAL_ENO()
    {
        Set_Order(2);
        Use_Maximum_Alpha();
    }

//#####################################################################
    void Conservation_Solver(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<VECTOR<T,d> ,VECTOR<int,1> >& U,ARRAY<VECTOR<T,d> ,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem_explicit,
        const VECTOR<bool,2>& outflow_boundaries,ARRAY<VECTOR<T,d> ,VECTOR<int,1> >* U_flux=0) override;
private:
    template<int order> void Conservation_Solver_Helper(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<VECTOR<T,d> ,VECTOR<int,1> >& U,ARRAY<VECTOR<T,d> ,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem,
        EIGENSYSTEM<T,TV_DIMENSION>& eigensystem_explicit,const VECTOR<bool,2>& outflow_boundaries);
//#####################################################################
};   
}
#endif
