//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DUALCONTOUR_3D
//##################################################################### 
#ifndef __DUALCONTOUR_3D__
#define __DUALCONTOUR_3D__

#include <Core/Arrays/ARRAY.h>
#include <Geometry/Level_Sets/LEVELSET.h>
namespace PhysBAM{

template<class T> class TRIANGULATED_SURFACE;

template<class T>
class DUALCONTOUR_3D
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    LEVELSET<TV>& levelset;
    ARRAY<VECTOR<int,4> > topology;
    ARRAY<int,TV_INT> vertices;
    ARRAY<TV> geometry;
    ARRAY<TV> normals;
    GRID<TV> grid; 

    DUALCONTOUR_3D(LEVELSET<TV>& levelset_input)
        :levelset(levelset_input)
    {
        if(levelset_input.grid.Is_MAC_Grid()) grid=levelset_input.grid.Get_Regular_Grid_At_MAC_Positions();else grid=levelset_input.grid;
    }
    DUALCONTOUR_3D(const DUALCONTOUR_3D&) = delete;
    void operator=(const DUALCONTOUR_3D&) = delete;

    void Dualcontour()
    {Generate_Topology();Generate_Vertices();Ensure_Vertices_In_Correct_Cells();}
    
    static TRIANGULATED_SURFACE<T>* Create_Triangulated_Surface_From_Levelset(LEVELSET<TV>& levelset)
    {DUALCONTOUR_3D<T> dualcontour(levelset);dualcontour.Dualcontour();return dualcontour.Get_Triangulated_Surface();} 

//#####################################################################
    void Generate_Topology();
    void Generate_Vertices();
    void Ensure_Vertices_In_Correct_Cells();
    TRIANGULATED_SURFACE<T>* Get_Triangulated_Surface();
//#####################################################################
};
}
#endif
