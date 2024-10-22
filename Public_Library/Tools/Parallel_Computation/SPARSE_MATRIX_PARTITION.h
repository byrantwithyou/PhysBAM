//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_MATRIX_PARTITION
//#####################################################################
#ifndef __SPARSE_MATRIX_PARTITION__
#define __SPARSE_MATRIX_PARTITION__

#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/INTERVAL.h>
namespace PhysBAM{
template<class T> class SPARSE_MATRIX_FLAT_MXN;

class SPARSE_MATRIX_PARTITION
{
public:
    int number_of_sides;
    INTERVAL<int> interior_indices;
    ARRAY<INTERVAL<int> > ghost_indices;
    ARRAY<ARRAY<int> > boundary_indices;
    int interior_offset;
    ARRAY<int> neighbor_ranks;
    ARRAY<SPARSE_MATRIX_PARTITION*> neighbors;

    SPARSE_MATRIX_PARTITION()
        :number_of_sides(0)
    {}

    SPARSE_MATRIX_PARTITION(const int number_of_sides_input)
    {
        Set_Number_Of_Sides(number_of_sides_input);
    }

    void Set_Number_Of_Sides(const int number_of_sides_input)
    {number_of_sides=number_of_sides_input;
    ghost_indices.Resize(number_of_sides);boundary_indices.Resize(number_of_sides);
    neighbor_ranks.Resize(number_of_sides);neighbors.Resize(number_of_sides);}

    int Interior_Rows() const
    {return interior_indices.Size();}

    template<class T>
    int Interior_Entries(const SPARSE_MATRIX_FLAT_MXN<T>& A) const
    {return A.offsets(interior_indices.max_corner)-A.offsets(interior_indices.min_corner);}

    void Set_Interior_Offset(const int previous_rows)
    {interior_offset=previous_rows-interior_indices.min_corner;}

    int Translate_Index(const int j) const
    {if(interior_indices.Lazy_Inside_Half_Open(j)) return Translate_Interior_Index(j);
    for(int r=0;r<number_of_sides-1;r++) if(ghost_indices(r).Lazy_Inside_Half_Open(j)) return Translate_Ghost_Index(j,r);
    return Translate_Ghost_Index(j,number_of_sides);}

    int Translate_Interior_Index(const int j) const
    {assert(interior_indices.Lazy_Inside_Half_Open(j));return j+interior_offset;}

    int Translate_Ghost_Index(const int j,const int region) const
    {assert(ghost_indices(region).Lazy_Inside_Half_Open(j));
    assert(ghost_indices(region).Size()==neighbors(region)->boundary_indices(region^1).m);
    return neighbors(region)->boundary_indices(region^1)(j-ghost_indices(region).min_corner+1)+neighbors(region)->interior_offset;}

//#####################################################################
};
}
#endif
