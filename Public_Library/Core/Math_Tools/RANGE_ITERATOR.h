//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Neil Molino, Igor Neverov, Duc Nguyen, Avi Robinson-Mosher,
//     Craig Schroeder, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANGE_ITERATOR
//#####################################################################
#ifndef __RANGE_ITERATOR__
#define __RANGE_ITERATOR__
#include <Core/Math_Tools/RANGE.h>
#include <Core/Utilities/PHYSBAM_ATTRIBUTE.h>
namespace PhysBAM{
template<int d>
struct RANGE_ITERATOR
{
    enum INTERNAL_ENUM {END_ITERATOR};
    typedef VECTOR<int,d> TV_INT;
    RANGE<TV_INT> domain;
    TV_INT index;

    RANGE_ITERATOR(const RANGE<TV_INT>& domain_input) PHYSBAM_ALWAYS_INLINE
        :domain(domain_input)
    {Reset();}

    RANGE_ITERATOR(const RANGE<TV_INT>& domain_input,INTERNAL_ENUM) PHYSBAM_ALWAYS_INLINE
        :domain(domain_input)
    {index=domain.max_corner;}

    bool Valid() const PHYSBAM_ALWAYS_INLINE
    {return index(d-1)<domain.max_corner(d-1);}

    void Next() PHYSBAM_ALWAYS_INLINE
    {index(d-1)++;if(!Valid()) Next_Helper();}

    void Next_Helper() PHYSBAM_ALWAYS_INLINE
    {
        index(d-1)=domain.min_corner(d-1);
        for(int i=d-2;i>=0;i--){
            if(++index(i)<domain.max_corner(i)) return;
            index(i)=domain.min_corner(i);}
        index(d-1)=domain.max_corner(d-1);
    }

    void Prev() PHYSBAM_ALWAYS_INLINE
    {index(d-1)--;if(index(d-1)<domain.min_corner(d-1)) Prev_Helper();}

    void Prev_Helper() PHYSBAM_ALWAYS_INLINE
    {
        index(d-1)=domain.max_corner(d-1)-1;
        for(int i=d-2;i>=0;i--){
            if(--index(i)>=domain.max_corner(i)) return;
            index(i)=domain.max_corner(i);}
        index(d-1)=domain.min_corner(d-1)-1;
    }

    void Reset() PHYSBAM_ALWAYS_INLINE
    {index=domain.min_corner;if(!index.All_Less(domain.max_corner)) index(d-1)=domain.max_corner(d-1);}

    // stl
    RANGE_ITERATOR& operator++() PHYSBAM_ALWAYS_INLINE
    {Next();return *this;}

    RANGE_ITERATOR operator++(int) PHYSBAM_ALWAYS_INLINE
    {RANGE_ITERATOR it(*this);Next();return it;}

    RANGE_ITERATOR& operator--() PHYSBAM_ALWAYS_INLINE
    {Prev();return *this;}

    RANGE_ITERATOR operator--(int) PHYSBAM_ALWAYS_INLINE
    {RANGE_ITERATOR it(*this);Prev();return it;}

    TV_INT* operator->() PHYSBAM_ALWAYS_INLINE
    {return &index;}

    const TV_INT* operator->() const PHYSBAM_ALWAYS_INLINE
    {return &index;}

    TV_INT& operator*() PHYSBAM_ALWAYS_INLINE
    {return index;}

    const TV_INT& operator*() const PHYSBAM_ALWAYS_INLINE
    {return index;}

    bool operator==(const RANGE_ITERATOR& it) const PHYSBAM_ALWAYS_INLINE
    {for(int i=d-1;i>=0;i--) if(index(i)!=it.index(i)) return false;return true;}

    bool operator!=(const RANGE_ITERATOR& it) const PHYSBAM_ALWAYS_INLINE
    {return !(*this==it);}

    bool operator<(const RANGE_ITERATOR& it) const PHYSBAM_ALWAYS_INLINE
    {for(int i=d-1;i>=0;i--) if(index(i)!=it.index(i)) return index(i)<it.index(i);return true;}

    bool operator>(const RANGE_ITERATOR& it) const PHYSBAM_ALWAYS_INLINE
    {return it<*this;}

    bool operator<=(const RANGE_ITERATOR& it) const PHYSBAM_ALWAYS_INLINE
    {return !(*this>it);}

    bool operator>=(const RANGE_ITERATOR& it) const PHYSBAM_ALWAYS_INLINE
    {return !(*this<it);}
};

template<>
struct RANGE_ITERATOR<0>
{
    enum INTERNAL_ENUM {END_ITERATOR};
    typedef VECTOR<int,0> TV_INT;
    RANGE<TV_INT> domain;
    TV_INT index;
    bool first;

    RANGE_ITERATOR(const RANGE<TV_INT>& domain_input) PHYSBAM_ALWAYS_INLINE
    {first=true;}

    RANGE_ITERATOR(const RANGE<TV_INT>& domain_input,INTERNAL_ENUM) PHYSBAM_ALWAYS_INLINE
        :domain(domain_input)
    {first=false;}

    bool Valid() const PHYSBAM_ALWAYS_INLINE
    {return first;}

    void Next() PHYSBAM_ALWAYS_INLINE
    {first=false;}

    void Prev() PHYSBAM_ALWAYS_INLINE
    {first=true;}

    void Reset() PHYSBAM_ALWAYS_INLINE
    {first=true;}

    // stl
    RANGE_ITERATOR& operator++() PHYSBAM_ALWAYS_INLINE
    {Next();return *this;}

    RANGE_ITERATOR operator++(int) PHYSBAM_ALWAYS_INLINE
    {RANGE_ITERATOR it(*this);Next();return it;}

    RANGE_ITERATOR& operator--() PHYSBAM_ALWAYS_INLINE
    {Prev();return *this;}

    RANGE_ITERATOR operator--(int) PHYSBAM_ALWAYS_INLINE
    {RANGE_ITERATOR it(*this);Prev();return it;}

    TV_INT* operator->() PHYSBAM_ALWAYS_INLINE
    {return &index;}

    const TV_INT* operator->() const PHYSBAM_ALWAYS_INLINE
    {return &index;}

    TV_INT& operator*() PHYSBAM_ALWAYS_INLINE
    {return index;}

    const TV_INT& operator*() const PHYSBAM_ALWAYS_INLINE
    {return index;}

    bool operator==(const RANGE_ITERATOR& it) const PHYSBAM_ALWAYS_INLINE
    {return first==it.first;}

    bool operator!=(const RANGE_ITERATOR& it) const PHYSBAM_ALWAYS_INLINE
    {return !(*this==it);}

    bool operator<(const RANGE_ITERATOR& it) const PHYSBAM_ALWAYS_INLINE
    {return first>it.first;}

    bool operator>(const RANGE_ITERATOR& it) const PHYSBAM_ALWAYS_INLINE
    {return it<*this;}

    bool operator<=(const RANGE_ITERATOR& it) const PHYSBAM_ALWAYS_INLINE
    {return !(*this>it);}

    bool operator>=(const RANGE_ITERATOR& it) const PHYSBAM_ALWAYS_INLINE
    {return !(*this<it);}
};

template<int d> RANGE_ITERATOR<d> begin(const RANGE<VECTOR<int,d> >& range)
{return RANGE_ITERATOR<d>(range);}

template<int d> RANGE_ITERATOR<d> end(const RANGE<VECTOR<int,d> >& range)
{return RANGE_ITERATOR<d>(range,RANGE_ITERATOR<d>::END_ITERATOR);}

}
#endif
