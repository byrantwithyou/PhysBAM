//#####################################################################
// Copyright 2017, Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_RANGE_ITERATOR
//#####################################################################
#ifndef __FACE_RANGE_ITERATOR__
#define __FACE_RANGE_ITERATOR__

#include <Core/Math_Tools/RANGE.h>
#include <Core/Utilities/PHYSBAM_ATTRIBUTE.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>

namespace PhysBAM{

// 
//   h   h   h   g   f   f   f   g   h   h   h
// h X e X e X d X c X c X c X c X d X e X e X h
//   e   e   e   c   c   c   c   c   e   e   e
// h X e X e X d X c X c X c X c X d X e X e X h
//   e   e   e   c   c   c   c   c   e   e   e
// h X e X e X d X c X c X c X c X d X e X e X h
//   d   d   d   b   b   b   b   b   d   d   d
// g X c X c X b O a O a O a O a O b X c X c X g
//   c   c   c   a   a   a   a   a   c   c   c
// f X c X c X b O a O a O a O a O b X c X c X f
//   c   c   c   a   a   a   a   a   c   c   c
// f X c X c X b O a O a O a O a O b X c X c X f
//   c   c   c   a   a   a   a   a   c   c   c
// f X c X c X b O a O a O a O a O b X c X c X f
//   c   c   c   a   a   a   a   a   c   c   c
// g X c X c X b O a O a O a O a O b X c X c X g
//   d   d   d   b   b   b   b   b   d   d   d
// h X e X e X d X c X c X c X c X d X e X e X h
//   e   e   e   c   c   c   c   c   e   e   e
// h X e X e X d X c X c X c X c X d X e X e X h
//   e   e   e   c   c   c   c   c   e   e   e
// h X e X e X d X c X c X c X c X d X e X e X h
//   h   h   h   g   f   f   f   g   h   h   h
// 
// Flags: (side=-1)
//   (default)                               *  b+c+d+e+f+g+h
//   interior                                    a+b+c+d+e+f+g+h
// 
//   skip_inner                              *  c+d+e+f+g+h
//   skip_outer                              *  b+c+d+e
//   skip_inner|skip_outer                   *  c+d+e
//   interior|skip_outer                        a+b+c+d+e
// 
//   duplicate_corners                          b+c+d+d+e+e+f+g+h+h
//   duplicate_corners|skip_inner               c+d+e+e+f+g+h+h
//   duplicate_corners|skip_outer               b+c+d+d+e+e
//   duplicate_corners|skip_inner|skip_outer    c+d+e+e
// 
// For combinations above with "*" these flags can be added:
//   delay_corners         Corners are normally covered with the first
//                         side that makes sense.  This flag delays the
//                         corners until the last side that makes sense.
//                         The faces covered does not change.
//   partial_single_side   This flag only affects the behavior when one
//                         side is being iterated.  Normally, when one
//                         side is being iterated, all corners are
//                         included.  This flag suppresses this.  Iterating
//                         one side covers the same faces as would be
//                         covered with that side when covering all sides.
// 
// Other flags:
//   reverse               Reverses the order in which faces are visited.
//                         Note that with this flag you must use Prev() and
//                         Prev_Valid().
//   end                   Like the reverse flag, but starts out one past
//                         the end.  Calling Prev() once will give you the
//                         iterator that the reverse flag returns.
//   side_mask             Allows you to iterate over a subset of axes.
//                         When this flag is set, the side input integer is
//                         a bitmask.  When set, the iterator always behaves
//                         as though multiple axes are being iterated.
//   axis_mask             Allows you to iterate over a subset of axes.
//                         When this flag is set, the axis input integer is
//                         a bitmask.

enum class RF
{
    none=0,
    interior=1,
    delay_corners=2,
    duplicate_corners=4,
    skip_inner=8,
    skip_outer=0x10,
    partial_single_side=0x20,
    reverse=0x40,
    end=0x80,
    axis_mask=0x100,
    side_mask=0x200
};
inline RF operator|(RF a,RF b){return RF((int)a|(int)b);}
inline RF operator&(RF a,RF b){return RF((int)a&(int)b);}
inline RF operator^(RF a,RF b){return RF((int)a^(int)b);}
inline RF operator~(RF a){return RF(~(int)a);}
inline bool operator!(RF a){return !(int)a;}
inline bool any(RF a){return (int)a;}

template<int d>
class FACE_RANGE_ITERATOR
{
    typedef VECTOR<int,d> TV_INT;

    TV_INT vecs[4];
    int axis_adj;
    int side_adj;
    int indices;

    TV_INT current[2];
public:
    int axis_mask;
    int side_mask;
    int side;
    FACE_INDEX<d> face;

    FACE_RANGE_ITERATOR(const RANGE<TV_INT>& outer,const RANGE<TV_INT>& inner,
        RF flags=RF::none,int side_input=-1,int axis=-1);
    FACE_RANGE_ITERATOR(const RANGE<TV_INT>& range,int outer_ghost,
        int inner_ghost,RF flags=RF::none,int side_input=-1,int axis=-1);
    FACE_RANGE_ITERATOR(const TV_INT& counts,int outer_ghost,int inner_ghost,
        RF flags=RF::none,int side_input=-1,int axis=-1);
    FACE_RANGE_ITERATOR(const RANGE<TV_INT>& range,
        RF flags=RF::none,int side_input=-1,int axis=-1); // implict RF::interior
    FACE_RANGE_ITERATOR(const TV_INT& counts,int outer_ghost,
        RF flags=RF::none,int side_input=-1,int axis=-1); // implict RF::interior
    ~FACE_RANGE_ITERATOR()=default;

    void Next() PHYSBAM_ALWAYS_INLINE
    {face.index(d-1)++;if(!Valid()) Next_Helper();}

    bool Valid() const PHYSBAM_ALWAYS_INLINE
    {return face.index(d-1)<current[1](d-1);}

    void Prev() PHYSBAM_ALWAYS_INLINE
    {face.index(d-1)--;if(!Prev_Valid()) Prev_Helper();}

    bool Prev_Valid() const PHYSBAM_ALWAYS_INLINE
    {return face.index(d-1)>=current[0](d-1);}

private:
    void Next_Helper();
    void Prev_Helper();
    void Next_Axis_Side();
    void Prev_Axis_Side();
    void Fill_Current();
    void Encode(int side_input,int axis,RF flags);
    void Initialize(int side_input,int axis,RF flags);
    void Reset(RF flags);
//#####################################################################
};
}
#endif
