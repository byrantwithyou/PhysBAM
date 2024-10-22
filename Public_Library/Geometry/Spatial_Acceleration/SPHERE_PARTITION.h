//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPHERE_PARTITION
//##################################################################### 
#ifndef __SPHERE_PARTITION__
#define __SPHERE_PARTITION__

#include <Core/Arrays/ARRAY.h>
#include <Core/Math_Tools/RANGE.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
namespace PhysBAM{

template<class T>
class SPHERE_PARTITION
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    ARRAY<SPHERE<TV> > spheres; // list of all the spheres
    GRID<TV> grid;
    RANGE<TV> box; // box containing the spheres
    ARRAY<ARRAY<int>*,TV_INT> voxel_sphere_list;
    
    SPHERE_PARTITION(const int number_input)
        :spheres(number_input)
    {
        Set_Up_Grid(TV_INT()+1);
    }
    SPHERE_PARTITION(const SPHERE_PARTITION&) = delete;
    void operator=(const SPHERE_PARTITION&) = delete;

    virtual ~SPHERE_PARTITION()
    {voxel_sphere_list.Delete_Pointers_And_Clean_Memory();}

//#####################################################################
    TV Normal(const TV& location,const int aggregate=0) const;
    bool Inside(const TV& location,const T thickness_over_two=0) const;
    bool Outside(const TV& location,const T thickness_over_two=0) const;
    bool Boundary(const TV& location,const T thickness_over_two) const;
    void Set_Up_Grid(const TV_INT& size); // use this to initialize!
    void Find_Voxel(const TV& location,int& i_left,int& j_bottom,int& ij_front) const;
//#####################################################################
};   
}
#endif

