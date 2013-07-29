//#####################################################################
// Copyright 2013, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __MPLE_ITERATOR__
#define __MPLE_ITERATOR__

#include <Tools/Vectors/VECTOR.h>
#include "MPLE_POINT.h"

namespace PhysBAM{

template<class TV,int w>
class MPLE_ITERATOR
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

protected:

    TV_INT index;
public:

    const MPLE_POINT<TV,w>& point;

    MPLE_ITERATOR(const MPLE_POINT<TV,w>& point_input)
        :point(point_input){}

    void Next()
    {
        index(TV::m-1)++;
        int i=TV::m-1;
        while(i>0 && index(i)>=w){
            index(i)=0;i--;index(i)++;}
    }

    bool Valid() const
    {return index(0)<w;}

    TV_INT Node() const
    {return point.base_node+index;}

    T Weight() const
    {
        T weight=1;
        for(int i=0;i<TV::m;i++)
            weight*=point.weights(index(i))(i);
        return weight;
    }
};

}
#endif
