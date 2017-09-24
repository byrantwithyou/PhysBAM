//#####################################################################
// Copyright 2017, Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANGE_ITERATOR
//#####################################################################
#ifndef __RANGE_ITERATOR__
#define __RANGE_ITERATOR__

#include <Core/Math_Tools/RANGE.h>
#include <Core/Utilities/PHYSBAM_ATTRIBUTE.h>

namespace PhysBAM{

// c c c b b b b b c c c
// c c c b b b b b c c c
// c c c b b b b b c c c
// b b b a a a a a b b b
// b b b a a a a a b b b
// b b b a a a a a b b b
// b b b a a a a a b b b
// b b b a a a a a b b b
// c c c b b b b b c c c
// c c c b b b b b c c c
// c c c b b b b b c c c
// 
// Flags: (side=-1)
//   (default)                                  a+b+c
//   ghost                                   *  b+c
//   duplicate_corners                          a+b+c+c
// 
// For combinations above with "*" these flags can be added:
//   delay_corners         Corners are normally covered with the first
//                         side that makes sense.  This flag delays the
//                         corners until the last side that makes sense.
//                         The cells covered does not change.
//   partial_single_side   This flag only affects the behavior when one
//                         side is being iterated.  Normally, when one
//                         side is being iterated, all corners are
//                         included.  This flag suppresses this.  Iterating
//                         one side covers the same cells as would be
//                         covered with that side when covering all sides.
// 
// Other flags:
//   reverse               Reverses the order in which cells are visited.
//                         Note that with this flag you must use Prev() and
//                         Prev_Valid().
//   end                   Like the reverse flag, but starts out one past
//                         the end.  Calling Prev() once will give you the
//                         iterator that the reverse flag returns.
//   side_mask             Allows you to iterate over a subset of axes.
//                         When this flag is set, the side input integer is
//                         a bitmask.  When set, the iterator always behaves
//                         as though multiple axes are being iterated.

enum class RI
{
    none=0,
    ghost=1,
    delay_corners=2,
    duplicate_corners=4,
    partial_single_side=0x20,
    reverse=0x40,
    end=0x80,
    side_mask=0x200
};
inline RI operator|(RI a,RI b){return RI((int)a|(int)b);}
inline RI operator&(RI a,RI b){return RI((int)a&(int)b);}
inline RI operator^(RI a,RI b){return RI((int)a^(int)b);}
inline RI operator~(RI a){return RI(~(int)a);}
inline bool operator!(RI a){return !(int)a;}
inline bool any(RI a){return (int)a;}

template<int d>
class RANGE_ITERATOR
{
    typedef VECTOR<int,d> TV_INT;

protected:
    TV_INT vecs[4];
    int side_adj;
    int indices;

    TV_INT current[2];
public:
    int side_mask;
    int side;
    TV_INT index;

protected:
    RANGE_ITERATOR()=default;
public:
    RANGE_ITERATOR(const RANGE<TV_INT>& outer,const RANGE<TV_INT>& inner,
        RI flags=RI::none,int side_input=-1);
    RANGE_ITERATOR(const RANGE<TV_INT>& range,int outer_ghost,
        int inner_ghost,RI flags=RI::none,int side_input=-1);
    RANGE_ITERATOR(const TV_INT& counts,int outer_ghost,int inner_ghost,
        RI flags=RI::none,int side_input=-1);
    explicit RANGE_ITERATOR(const RANGE<TV_INT>& range,RI flags=RI::none);
    explicit RANGE_ITERATOR(const TV_INT& counts,int outer_ghost=0,
        RI flags=RI::none);
    ~RANGE_ITERATOR()=default;

    void Next() PHYSBAM_ALWAYS_INLINE
    {index(d-1)++;if(!Valid()) Next_Helper();}

    bool Valid() const PHYSBAM_ALWAYS_INLINE
    {return index(d-1)<current[1](d-1);}

    void Prev() PHYSBAM_ALWAYS_INLINE
    {index(d-1)--;if(!Prev_Valid()) Prev_Helper();}

    bool Prev_Valid() const PHYSBAM_ALWAYS_INLINE
    {return index(d-1)>=current[0](d-1);}

    void Set_Range(const RANGE<TV_INT>& outer,const RANGE<TV_INT>& inner);
    void Set_Range(const RANGE<TV_INT>& range,int outer_ghost,int inner_ghost);
    void Set_Range(const TV_INT& counts,int outer_ghost,int inner_ghost);
    void Initialize(RI flags,int side_input);

private:
    void Next_Helper();
    void Prev_Helper();
    void Next_Side();
    void Prev_Side();
    void Fill_Current();
    void Encode(RI flags,int side_input);
    void Reset(RI flags);
//#####################################################################
};


template<>
class RANGE_ITERATOR<0>
{
    typedef VECTOR<int,0> TV_INT;
public:
    bool first;
    TV_INT index;
    
protected:
    RANGE_ITERATOR()=default;
public:
    RANGE_ITERATOR(const RANGE<TV_INT>& outer,const RANGE<TV_INT>& inner,
        RI flags=RI::none,int side_input=-1){first=!(flags&RI::end);}
    RANGE_ITERATOR(const RANGE<TV_INT>& range,int outer_ghost,
        int inner_ghost,RI flags=RI::none,int side_input=-1){first=!(flags&RI::end);}
    RANGE_ITERATOR(const TV_INT& counts,int outer_ghost,int inner_ghost,
        RI flags=RI::none,int side_input=-1){first=!(flags&RI::end);}
    explicit RANGE_ITERATOR(const RANGE<TV_INT>& range,
        RI flags=RI::none){first=!(flags&RI::end);}
    explicit RANGE_ITERATOR(const TV_INT& counts,int outer_ghost=0,
        RI flags=RI::none){first=!(flags&RI::end);}
    ~RANGE_ITERATOR()=default;

    void Next() PHYSBAM_ALWAYS_INLINE
    {first=false;}

    bool Valid() const PHYSBAM_ALWAYS_INLINE
    {return first;}

    void Prev() PHYSBAM_ALWAYS_INLINE
    {first=true;}

    bool Prev_Valid() const PHYSBAM_ALWAYS_INLINE
    {return !first;}

    void Set_Range(const RANGE<TV_INT>& outer,const RANGE<TV_INT>& inner){}
    void Set_Range(const RANGE<TV_INT>& range,int outer_ghost,int inner_ghost){}
    void Set_Range(const TV_INT& counts,int outer_ghost,int inner_ghost){}
    void Initialize(RI flags,int side_input){first=!(flags&RI::end);}
//#####################################################################
};

template<int d> RANGE_ITERATOR<d> begin(const RANGE<VECTOR<int,d> >& range)
{return RANGE_ITERATOR<d>(range);}

template<int d> RANGE_ITERATOR<d> end(const RANGE<VECTOR<int,d> >& range)
{return RANGE_ITERATOR<d>(range,RI::end);}

}
#endif
