//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CELL_MAPPING
//#####################################################################
#ifndef __CELL_MAPPING__
#define __CELL_MAPPING__

#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Math_Tools/RANGE.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class CELL_MAPPING:public NONCOPYABLE
{
public:
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    const GRID<TV>& grid;

    ARRAY<int,TV_INT> cell_index[2];
    int next_index;

    VECTOR<bool,TV::m> periodic;

    int Get_Index(TV_INT index, bool is_neg)
    {
        ARRAY<int,TV_INT>& ar=cell_index[is_neg];
        for(int j=0;j<TV::m;j++)
            if(periodic(j)){
                if(index(j)<ar.domain.min_corner(j)) index(j)+=ar.domain.Edge_Lengths()(j);
                if(index(j)>=ar.domain.max_corner(j)) index(j)-=ar.domain.Edge_Lengths()(j);}

        int& ci=ar(index);
        if(ci<0) ci=next_index++;
        return ci;
    }

    int Get_Index_Fixed(TV_INT index, bool is_neg) const
    {
        const ARRAY<int,TV_INT>& ar=cell_index[is_neg];
        for(int j=0;j<TV::m;j++)
            if(periodic(j)){
                if(index(j)<ar.domain.min_corner(j)) index(j)+=ar.domain.Edge_Lengths()(j);
                if(index(j)>=ar.domain.max_corner(j)) index(j)-=ar.domain.Edge_Lengths()(j);}

        return ar(index);
    }

    void Shift(const int shift)
    {
        for(int s=0;s<2;s++)
            for(RANGE_ITERATOR<TV::m> it(cell_index[s].domain);it.Valid();it.Next())
                if(cell_index[s](it.index)>=0) cell_index[s](it.index)+=shift;
    }

    CELL_MAPPING(const GRID<TV>& grid_input)
        :grid(grid_input),next_index(0)
    {for(int i=0;i<2;i++){cell_index[i].Resize(grid.Domain_Indices());cell_index[i].Fill(-1);}}

    ~CELL_MAPPING(){}
};
}
#endif
