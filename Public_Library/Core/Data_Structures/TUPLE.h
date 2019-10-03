//#####################################################################
// Copyright 2003-2007, Geoffrey Irving, Frank Losasso, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TUPLE
//##################################################################### 
#ifndef __TUPLE__
#define __TUPLE__

#include <Core/Data_Structures/HASH_REDUCE.h>
#include <iostream>
#include <tuple>
namespace PhysBAM{

template<class... Args>
struct TUPLE_IO
{
    static const int n=std::tuple_size<std::tuple<Args...> >::value;
    template<int i> static enable_if_t<(i<n-1)>
    P(std::ostream& o,const std::tuple<Args...>& object)
    {
        o<<std::get<i>(object)<<" ";
        P<i+1>(o,object);
    }
    template<int i> static enable_if_t<(i==n-1)>
    P(std::ostream& o,const std::tuple<Args...>& object)
    {
        o<<std::get<i>(object);
    }
};
template<class... Args>
inline std::ostream& operator<<(std::ostream& o,const std::tuple<Args...>& object)
{
    o<<"(";
    TUPLE_IO<Args...>::template P<0>(o,object);
    o<<")";
    return o;
}

template<class... Args> struct HASH_REDUCE<std::tuple<Args...> >
{
    template<int n>
    static enable_if_t<(n>3),int>
    F(const std::tuple<Args...>& key)
    {
        return int_hash(F<n-2>(key),
            Hash_Reduce(std::get<n-2>(key)),
            Hash_Reduce(std::get<n-1>(key)));
    }

    template<int n>
    static enable_if_t<(n==3),int>
    F(const std::tuple<Args...>& key)
    {
        return int_hash(
            Hash_Reduce(std::get<0>(key)),
            Hash_Reduce(std::get<1>(key)),
            Hash_Reduce(std::get<2>(key)));
    }

    template<int n>
    static enable_if_t<(n==2),int>
    F(const std::tuple<Args...>& key)
    {
        return int_hash(
            Hash_Reduce(std::get<0>(key)),
            Hash_Reduce(std::get<1>(key)));
    }

    template<int n>
    static enable_if_t<(n==1),int>
    F(const std::tuple<Args...>& key)
    {
        return int_hash(Hash_Reduce(std::get<0>(key)));
    }

    template<int n>
    static enable_if_t<(n==0),int>
    F(const std::tuple<Args...>& key)
    {
        return missing_element_hash;
    }

    static int H(const std::tuple<Args...>& key)
    {
        return F<std::tuple_size<std::tuple<Args...> >::value>(key);
    }
};

}
#endif
