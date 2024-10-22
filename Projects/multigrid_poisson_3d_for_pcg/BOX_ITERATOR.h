//#####################################################################
// Copyright 2009, Eftychios Sifakis, Yongning Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOX_ITERATOR
//#####################################################################
#ifndef __BOX_ITERATOR__
#define __BOX_ITERATOR__

#include <Core/Math_Tools/RANGE.h>
namespace PhysBAM{

template<int d,int stride=1> class BOX_ITERATOR;

template<int d>
class BOUNDARY_ITERATOR
{
    typedef VECTOR<int,d> TV_INT;
    
    const RANGE<TV_INT>& box;
    TV_INT index;
    RANGE<TV_INT> regions[2*d];

    int current_region;
    int side; // side=1: xmin, side=2 xmax, side=3:ymin, side=4:ymax
    bool valid;
public:
    BOUNDARY_ITERATOR(const RANGE<TV_INT>& box_input,int side=0)
        :box(box_input),side(side)
    {
        assert(side<2*d);
        TV_INT min_corner=box.min_corner+1;
        TV_INT max_corner=box.max_corner-1;
        for(int v=0;v<d;v++){
            min_corner(v)=box.min_corner(v);
            max_corner(v)=box.min_corner(v);
            regions[v*2]=RANGE<TV_INT>(min_corner,max_corner);
            min_corner(v)=box.max_corner(v);
            max_corner(v)=box.max_corner(v);
            regions[v*2+1]=RANGE<TV_INT>(min_corner,max_corner);
            min_corner(v)=box.min_corner(v);
        }
        if(side)
            Reset(side-1);
        else
            Reset();
    }        
    
    void Reset(const int region_index=0)
    {
        valid=side?region_index<side:(region_index<2*d);
        if(valid){
            current_region=region_index;
            index=regions[current_region].min_corner;
        }
    }

    bool Valid() const
    {
        return valid;
    }

    void Next()
    {
        if(!valid) return;
        for(int i=d-1;i>=0;i--){
            if(index(i)<regions[current_region].max_corner(i)){
                index(i)++; return;
            }
            else
                index(i)=regions[current_region].min_corner(i);
        }
        Reset(current_region+1);
    }
        
    const TV_INT&  Index()
    {return index;}
        

};

//#####################################################################
}
#endif
