//#####################################################################
// Copyright 2003-2009, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUBIC_MN_INTERPOLATION 
//#####################################################################
#ifndef __CUBIC_MN_INTERPOLATION__
#define __CUBIC_MN_INTERPOLATION__

#include <Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
namespace PhysBAM{

template<class T,class T2>
class CUBIC_MN_INTERPOLATION
{
private:
    T m0,m1,m3,n0,n1,n2,n3; // used for acceleration
public:

    CUBIC_MN_INTERPOLATION();
    ~CUBIC_MN_INTERPOLATION();
    void Set_Parameters(const T b=1./3,const T c=1./3);
    T2 Cubic_MN(const T2& u_0,const T2& u_1,const T2& u_2,const T2& u_3,const T alpha) const;

    ARRAY<T> Cubic_MN_Weights(const T alpha) const;

//#####################################################################
};
}
#endif
