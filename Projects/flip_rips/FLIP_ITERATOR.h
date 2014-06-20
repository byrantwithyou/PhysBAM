//#####################################################################
// Copyright 2014, Alexey Stomakhin
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLIP_ITERATOR__
#define __FLIP_ITERATOR__

#include "FLIP_PARTICLE.h"

namespace PhysBAM{

template<class TV,int w>
class FLIP_ITERATOR
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;

protected:

    TV_INT index;

public:

    const FLIP_PARTICLE<TV,w>& particle;
    const int axis;

    FLIP_ITERATOR(const FLIP_PARTICLE<TV,w>& particle_input,const int axis_input)
        :particle(particle_input),axis(axis_input){}

    void Next()
    {
        index(TV::m-1)++;
        int i=TV::m-1;
        while(i>0 && index(i)>=w){
            index(i)=0;i--;index(i)++;}
    }

    bool Valid() const
    {return index(0)<w;}

    TV_INT Index() const
    {return particle.base_index(axis)+index;}

    T Weight() const
    {
        T weight=1;
        for(int i=0;i<TV::m;i++)
            weight*=particle.weights(axis)(index(i))(i);
        return weight;
    }

    inline TV Weight_Gradient() const
    {
        TV weight_gradient=TV()+1;
        for(int i=0;i<TV::m;i++)
        for(int j=0;j<TV::m;j++)
            weight_gradient(i)*=(i==j)?particle.dweight(axis)(index(j))(j):particle.weights(axis)(index(j))(j);
        return weight_gradient;
    }
};
}
#endif
