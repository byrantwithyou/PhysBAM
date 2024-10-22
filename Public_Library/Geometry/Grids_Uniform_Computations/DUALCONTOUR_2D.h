//#####################################################################
// Copyright 2004-2008, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DUALCONTOUR_2D
//#####################################################################
#ifndef __DUALCONTOUR_2D__
#define __DUALCONTOUR_2D__

#include <Core/Arrays/ARRAY.h>
#include <Geometry/Level_Sets/LEVELSET.h>
namespace PhysBAM{

template<class T> class SEGMENTED_CURVE_2D;
template<class T> class TRIANGULATED_AREA;

template<class T>
class DUALCONTOUR_2D
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
private:
    LEVELSET<TV>& levelset;
    T contour_value;
    ARRAY<VECTOR<int,2> > topology;
    ARRAY<int,TV_INT> vertices;
    ARRAY<TV> geometry;
    ARRAY<TV> normals;
    GRID<TV> grid;
    bool is_distance_field;
public:

    DUALCONTOUR_2D(LEVELSET<TV>& levelset_input,const T contour_value_input=0,const bool is_distance_field_input=true)
        :levelset(levelset_input),contour_value(contour_value_input),is_distance_field(is_distance_field_input)
    {if(levelset_input.grid.Is_MAC_Grid()) grid=levelset_input.grid.Get_Regular_Grid_At_MAC_Positions();else grid=levelset_input.grid;}

    DUALCONTOUR_2D(const DUALCONTOUR_2D&) = delete;
    void operator=(const DUALCONTOUR_2D&) = delete;

    void Dualcontour()
    {Generate_Topology();Generate_Vertices();Ensure_Vertices_In_Correct_Cells();}
    
    static SEGMENTED_CURVE_2D<T>* Create_Segmented_Curve_From_Levelset(LEVELSET<TV>& levelset,const T contour_value=0,const bool is_distance_field_input=true)
    {DUALCONTOUR_2D<T> dualcontour(levelset,contour_value,is_distance_field_input);dualcontour.Dualcontour();return dualcontour.Get_Segmented_Curve();} 

    static TRIANGULATED_AREA<T>* Create_Triangulated_Area_From_Levelset(LEVELSET<TV>& levelset,const int sign=-1,const T contour_value=0,const bool is_distance_field_input=true)
    {DUALCONTOUR_2D<T> dualcontour(levelset,contour_value,is_distance_field_input);dualcontour.Dualcontour();return dualcontour.Get_Triangulated_Area(sign);} 

//#####################################################################
    void Generate_Topology();
    void Generate_Vertices();
    void Ensure_Vertices_In_Correct_Cells();
    SEGMENTED_CURVE_2D<T>* Get_Segmented_Curve();
    TRIANGULATED_AREA<T>* Get_Triangulated_Area(const int sign=-1);
//#####################################################################
};
}
#endif
