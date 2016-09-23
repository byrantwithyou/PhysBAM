//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARAMETER_SPACE
//#####################################################################
#ifndef __PARAMETER_SPACE__
#define __PARAMETER_SPACE__

#include <Core/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

template<class T>
struct PARAMETER_SPACE
{
    virtual ~PARAMETER_SPACE(){}

    virtual PARAMETER_SPACE<T>& Zero_Clone() const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Op(T a,const PARAMETER_SPACE& x,T b,const PARAMETER_SPACE& y)=0; // this=a*x+b*y
    virtual void Copy(const PARAMETER_SPACE& x) {Op(1,x,0,*this);} // this=x
    virtual void Zero() {Op(0,*this,0,*this);} // this=0
    virtual T Dot(const PARAMETER_SPACE& x) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();} // this dot x
};
}
#endif
