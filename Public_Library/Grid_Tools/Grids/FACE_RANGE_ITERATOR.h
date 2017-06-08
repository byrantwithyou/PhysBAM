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
        RF flags=RF::none,int axis=-1,int side_input=-1);
    FACE_RANGE_ITERATOR(const RANGE<TV_INT>& range,int outer_ghost,
        int inner_ghost,RF flags=RF::none,int axis=-1,int side_input=-1);
    FACE_RANGE_ITERATOR(const TV_INT& counts,int outer_ghost,int inner_ghost,
        RF flags=RF::none,int axis=-1,int side_input=-1);
    FACE_RANGE_ITERATOR(const RANGE<TV_INT>& range,
        RF flags=RF::none,int axis=-1,int side_input=-1); // implict RF::interior
    FACE_RANGE_ITERATOR(const TV_INT& counts,int outer_ghost,
        RF flags=RF::none,int axis=-1,int side_input=-1); // implict RF::interior
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
    void Encode(int axis,int side_input,RF flags);
    void Initialize(int axis,int side_input,RF flags);
    void Reset(RF flags);
//#####################################################################
};
}
#endif
