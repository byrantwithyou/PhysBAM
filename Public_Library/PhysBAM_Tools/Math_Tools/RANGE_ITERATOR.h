//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Neil Molino, Igor Neverov, Duc Nguyen, Avi Robinson-Mosher,
//     Craig Schroeder, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANGE_ITERATOR
//#####################################################################
#ifndef __RANGE_ITERATOR__
#define __RANGE_ITERATOR__
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
namespace PhysBAM{
template<int d>
struct RANGE_ITERATOR
{
    typedef VECTOR<int,d> TV_INT;
    RANGE<TV_INT> domain;
    TV_INT index;

    RANGE_ITERATOR(const RANGE<TV_INT>& domain_input)
        :domain(domain_input)
    {Reset();}

    bool Valid() const PHYSBAM_ALWAYS_INLINE
    {return index(d-1)<domain.max_corner(d-1);}

    void Next() PHYSBAM_ALWAYS_INLINE
    {index(d-1)++;if(!Valid()) Next_Helper();}

    void Next_Helper()
    {
        index(d-1)=domain.min_corner(d-1);
        for(int i=d-2;i>=0;i--){
            if(++index(i)<domain.max_corner(i)) return;
            index(i)=domain.min_corner(i);}
        index(d-1)=domain.max_corner(d-1);
    }

    void Reset()
    {index=domain.min_corner;if(!index.All_Less(domain.max_corner)) index(d-1)=domain.max_corner(d-1);}
};

template<>
struct RANGE_ITERATOR<0>
{
    typedef VECTOR<int,0> TV_INT;
    RANGE<TV_INT> domain;
    TV_INT index;
    bool first;

    RANGE_ITERATOR(const RANGE<TV_INT>& domain_input)
    {first=true;}

    bool Valid() const PHYSBAM_ALWAYS_INLINE
    {return first;}

    void Next() PHYSBAM_ALWAYS_INLINE
    {first=false;}

    void Reset()
    {first=true;}
};
}
#endif
