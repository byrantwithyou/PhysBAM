//#####################################################################
// Copyright 2009, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FINITE_DIFFERENCE_DISCRETIZATION
//#####################################################################
#ifndef __FINITE_DIFFERENCE_DISCRETIZATION__
#define __FINITE_DIFFERENCE_DISCRETIZATION__

#include "BOX_ITERATOR.h"

namespace PhysBAM{

template<class T,int d> class FINITE_DIFFERENCE_DISCRETIZATION;
template<class T,int d> class STENCIL;

template<class T,int d>
class FINITE_DIFFERENCE_DISCRETIZATION
{
    typedef VECTOR<int,d> T_INDEX;
    typedef STENCIL<T,d> T_STENCIL;

public:
    const T h;

private:
    mutable T dummy;
public:

    FINITE_DIFFERENCE_DISCRETIZATION(const T h_input):h(h_input),dummy(0) {}
    
    // Stencil access helpers

    template<class ValidFunctor>
    T& Get_Or_Insert(T_STENCIL& stencil,const T_INDEX& index,ValidFunctor& valid) const
    {if(!valid(index)) PHYSBAM_FATAL_ERROR();return stencil.Get_Or_Insert(index);}

    template<class ValidFunctor,class ExistsFunctor>
    void Insert(T_STENCIL& stencil,const T data,const T_INDEX& index,ValidFunctor valid,ExistsFunctor exists) const
     {if(!valid(index)) PHYSBAM_FATAL_ERROR();if(exists(index)) stencil.Insert(index,data);}

    // Differentiation stencils

    template<class ValidFunctor>
    void Add_Forward_Difference_First_Derivative_Stencil(const int d1,const T_INDEX& index,T_STENCIL& stencil,const T scale,ValidFunctor valid) const
    {PHYSBAM_ASSERT((d1>=1 && d1<=d));
    Get_Or_Insert(stencil,index+T_INDEX::Axis_Vector(d1),valid)+=scale/h;
    Get_Or_Insert(stencil,index,valid)-=scale/h;}

    template<class ValidFunctor>
    void Add_Backward_Difference_First_Derivative_Stencil(const int d1,const T_INDEX& index,T_STENCIL& stencil,const T scale,ValidFunctor valid) const
    {PHYSBAM_ASSERT((d1>=1 && d1<=d));
    Get_Or_Insert(stencil,index,valid)+=scale/h;
    Get_Or_Insert(stencil,index-T_INDEX::Axis_Vector(d1),valid)-=scale/h;}

    // Interpolation stencils

    template<class ValidFunctor,class ExistsFunctor>
    void Set_Offset_Multilinear_Interpolation_Stencil(const T_INDEX& index,T_STENCIL& stencil,const T scale,ValidFunctor valid,ExistsFunctor exists) const
    {stencil.Remove_All();
    for(BOX_ITERATOR<d> iterator(RANGE<T_INDEX>(-T_INDEX::All_Ones_Vector(),2*T_INDEX::All_Ones_Vector()));iterator.Valid();iterator.Next()){
        const T_INDEX dindex=iterator.Index();
        T fscale=1;for(int v=1;v<=d;v++) switch(dindex(v)){case -1:case 2:fscale*=.25;break;;case 0:case 1:fscale*=.75;}
        Insert(stencil,scale*fscale,index+dindex,valid,exists);}}

//#####################################################################
};
}
#endif
