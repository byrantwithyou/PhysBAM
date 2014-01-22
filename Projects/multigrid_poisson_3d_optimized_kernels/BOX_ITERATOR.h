//#####################################################################
// Copyright 2009, Eftychios Sifakis, Yongning Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOX_ITERATOR
//#####################################################################
#ifndef __BOX_ITERATOR__
#define __BOX_ITERATOR__

#include <Tools/Math_Tools/RANGE.h>

namespace PhysBAM{

template<int d,int stride=1> class BOX_ITERATOR;

template<int d> 
class BOX_ITERATOR<d,1>
{
    typedef VECTOR<int,d> TV_INT;

    const RANGE<TV_INT>& box;
    TV_INT index;

public:
    BOX_ITERATOR(const RANGE<TV_INT>& box_input)
        :box(box_input)
    {
        Reset();
    }

    void Reset()
    {index=box.min_corner;}

    bool Valid() const
    {return index.x<=box.max_corner.x;}

    void Next()
    {for(int i=d;i>=1;i--) if(index(i)<box.max_corner(i) || i==1){index(i)++;return;} else index(i)=box.min_corner(i);}

    const TV_INT& Index()
    {return index;}
};

template<int d,int stride>
class BOX_ITERATOR
{
    STATIC_ASSERT((stride!=1));
    typedef VECTOR<int,d> TV_INT;

    const RANGE<TV_INT>& box;
    TV_INT index;

public:
    BOX_ITERATOR(const RANGE<TV_INT>& box_input)
        :box(box_input)
    {
        Reset();
    }

    void Reset()
    {index=box.min_corner;}

    bool Valid() const
    {return index.x<=box.max_corner.x;}

    void Next()
    {for(int i=d;i>=1;i--) if(index(i)+stride<=box.max_corner(i) || i==1){index(i)+=stride;return;} else index(i)=box.min_corner(i);}

    const TV_INT& Index()
    {return index;}
};

//#####################################################################
}
#endif
