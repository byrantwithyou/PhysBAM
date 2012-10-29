//#####################################################################
// Copyright 2003-2007, Geoffrey Irving, Frank Losasso, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIPLE
//##################################################################### 
#ifndef __TRIPLE__
#define __TRIPLE__

#include <PhysBAM_Tools/Data_Structures/HASH_REDUCE.h>
#include <PhysBAM_Tools/Read_Write/READ_WRITE_FUNCTIONS.h>
namespace PhysBAM{
template<class T1,class T2,class T3>
class TRIPLE
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;

    T1 x;T2 y;T3 z;

    TRIPLE(int input=0) 
        :x(T1()),y(T2()),z(T3())
    {}

    TRIPLE(const T1& x_input,const T2& y_input,const T3& z_input) 
        :x(x_input),y(y_input),z(z_input)
    {}

    bool operator==(const TRIPLE& t) const
    {return x==t.x && y==t.y && z==t.z;}

    bool operator!=(const TRIPLE& t) const
    {return !(*this==t);}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,x,y,z);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,x,y,z);}

//#####################################################################
};  
template<class T1,class T2,class T3>
inline TRIPLE<T1,T2,T3> Tuple(const T1& x,const T2& y,const T3& z)
{return TRIPLE<T1,T2,T3>(x,y,z);}

template<class T1,class T2,class T3>
inline std::ostream& operator<<(std::ostream& output,const TRIPLE<T1,T2,T3>& object)
{output<<"("<<object.x<<" "<<object.y<<" "<<object.z<<")";return output;}
template<class T1,class T2,class T3> struct HASH_REDUCE<TRIPLE<T1,T2,T3> >
{static int H(const TRIPLE<T1,T2,T3>& key){return int_hash(HASH_REDUCE<T1>::H(key.x),HASH_REDUCE<T2>::H(key.y),HASH_REDUCE<T3>::H(key.z));}};
}
#endif
