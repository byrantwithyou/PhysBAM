//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEXAHEDRON
//##################################################################### 
#ifndef __HEXAHEDRON__
#define __HEXAHEDRON__

#include <Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class HEXAHEDRON
{
public:
    VECTOR<T,3> x0,x1,x2,x3,x4,x5,x6,x7; 
    
    HEXAHEDRON()
        :x0(-1,-1,-1),x1(-1,-1,1),x2(-1,1,-1),x3(-1,1,1),x4(1,-1,-1),x5(1,-1,1),x6(1,1,-1),x7(1,1,1)
    {}

    HEXAHEDRON(const VECTOR<T,3>& x1_input,const VECTOR<T,3>& x2_input,const VECTOR<T,3>& x3_input,const VECTOR<T,3>& x4_input,
        const VECTOR<T,3>& x5_input,const VECTOR<T,3>& x6_input,const VECTOR<T,3>& x7_input,const VECTOR<T,3>& x8_input)
        :x0(x1_input),x1(x2_input),x2(x3_input),x3(x4_input),x4(x5_input),x5(x6_input),x6(x7_input),x7(x8_input)
    {}

    static T Volume(const VECTOR<T,3>& x0,const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,const VECTOR<T,3>& x3,
                    const VECTOR<T,3>& x4,const VECTOR<T,3>& x5,const VECTOR<T,3>& x6,const VECTOR<T,3>& x7)
    {return abs(Signed_Volume(x0,x1,x2,x3,x4,x5,x6,x7));}

    static T Signed_Volume(const VECTOR<T,3>& x0,const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,const VECTOR<T,3>& x3,
                           const VECTOR<T,3>& x4,const VECTOR<T,3>& x5,const VECTOR<T,3>& x6,const VECTOR<T,3>& x7)
    {VECTOR<T,3> x17=x7-x1,x12=x2-x1,x13=x3-x1,x24=x4-x2,x47=x7-x2,x46=x6-x4,x14=x4-x1,x15=x5-x1,x04=x4-x0;
    return (T)one_sixth*(VECTOR<T,3>::Triple_Product(x17,x12,x13)-VECTOR<T,3>::Triple_Product(x24,x47,x46)+VECTOR<T,3>::Triple_Product(x14,x17,x15)
        +VECTOR<T,3>::Triple_Product(x12,x17,x14)-VECTOR<T,3>::Triple_Product(x14,x24,x04));}

//#####################################################################
    T Volume() const;
    T Signed_Volume() const;
//#####################################################################
};   
}
#endif

