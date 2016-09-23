//#####################################################################
// Copyright 2003-2008, Geoffrey Irving, Michael Lentine, Frank Losasso, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PAIR
//##################################################################### 
#ifndef __PAIR__
#define __PAIR__

#include <Core/Data_Structures/HASH_REDUCE.h>
#include <Core/Math_Tools/choice.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <iostream>
namespace PhysBAM{

template<class T1,class T2>
class PAIR
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;

    template<int i> struct GET{typedef typename conditional<i==0,T1,T2>::type TYPE;};

    T1 x;T2 y;

    PAIR() 
        :x(T1()),y(T2()) 
    {}

    PAIR(const T1& x_input,const T2& y_input) 
        :x(x_input),y(y_input) 
    {}

    ~PAIR()
    {}

    bool operator==(const PAIR& p) const
    {return x==p.x && y==p.y;}

    bool operator!=(const PAIR& p) const
    {return !(*this==p);}

    bool operator<(const PAIR& p) const
    {return x<p.x || (x==p.x && y<p.y);}
    
    bool operator>(const PAIR& p) const
    {return x>p.x || (x==p.x && y>p.y);}

    template<int i> typename GET<i>::TYPE& Get()
    {return choice<i>(x,y);}

    template<int i> const typename GET<i>::TYPE& Get() const
    {return choice<i>(x,y);}

    void Get(T1& a,T2& b) const
    {a=x;b=y;}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,x,y);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,x,y);}

//#####################################################################
};  
template<class T1,class T2>
inline PAIR<T1,T2> Tuple(const T1& x,const T2& y)
{return PAIR<T1,T2>(x,y);}

template<class T1,class T2>
inline std::istream& operator>>(std::istream& input,PAIR<T1,T2>& p)
{FILE_UTILITIES::Ignore(input,'(');input>>p.x>>p.y;FILE_UTILITIES::Ignore(input,')');return input;}

template<class T1,class T2>
inline std::ostream& operator<<(std::ostream& output,const PAIR<T1,T2>& p)
{output<<"("<<p.x<<" "<<p.y<<")";return output;}
template<class T1,class T2> struct HASH_REDUCE<PAIR<T1,T2> >
{static int H(const PAIR<T1,T2>& key){return int_hash(HASH_REDUCE<T1>::H(key.x),HASH_REDUCE<T2>::H(key.y));}};
}
#endif
