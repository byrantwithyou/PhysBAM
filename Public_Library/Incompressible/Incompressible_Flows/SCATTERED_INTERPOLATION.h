//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SCATTERED_INTERPOLATION
//#####################################################################
#ifndef __SCATTERED_INTERPOLATION__
#define __SCATTERED_INTERPOLATION__

#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Math_Tools/RANGE.h>
namespace PhysBAM{

template<class TV>
class SCATTERED_INTERPOLATION
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<ARRAY<int>,TV_INT> T_ARRAYS_ARRAY_INT;
public:
    T radius_of_influence;
    bool use_tent_weights;
    bool use_distance_averaged_weights;
private:
    T radius_of_influence_squared;
public:

    SCATTERED_INTERPOLATION();
    ~SCATTERED_INTERPOLATION();

    void Set_Radius_Of_Influence(const T radius_of_influence_input)
    {radius_of_influence=radius_of_influence_input;radius_of_influence_squared=sqr(radius_of_influence);}

    void Use_Tent_Weights(const bool use=true)
    {use_tent_weights=use;if(use_tent_weights) use_distance_averaged_weights=false;}

    void Use_Distance_Averaged_Weights(const bool use=true)
    {use_distance_averaged_weights=use;if(use_distance_averaged_weights) use_tent_weights=false;}

private:
    RANGE<TV_INT> Grid_Influence_Bounds(const GRID<TV>& grid,const TV& X) const
    {return RANGE<TV_INT>(grid.Clamp_To_Cell(RANGE<TV>(X).Thickened(radius_of_influence)));}
public:

    template<class T_ARRAY_T2,class T_ARRAYS_T2> void Transfer_To_Grid(ARRAY_VIEW<const TV> domain,const T_ARRAY_T2& range,const GRID<TV>& grid,T_ARRAYS_T2& grid_data) const
    {Transfer_To_Grid_Helper<T_ARRAYS_T2>(domain,range,grid,grid_data);}

//#####################################################################
    template<class T_ARRAYS_T2> void Transfer_To_Grid_Helper(ARRAY_VIEW<const TV> domain,ARRAY_VIEW<const typename T_ARRAYS_T2::ELEMENT> range,const GRID<TV>& grid,T_ARRAYS_T2& grid_data) const;
private:
    static void Bin_Domain_Values(ARRAY_VIEW<const TV> domain_values,const GRID<TV>& grid,T_ARRAYS_ARRAY_INT& points_in_cell);
    template<class T_ARRAYS_T2> void Transfer_With_Distance_Averaged_Weights(ARRAY_VIEW<const TV> domain,ARRAY_VIEW<const typename T_ARRAYS_T2::ELEMENT> range,const GRID<TV>& grid,
        T_ARRAYS_T2& grid_data,const T_ARRAYS_ARRAY_INT& points_in_cell) const;
    template<class T_ARRAYS_T2> void Transfer_With_Tent_Weights(ARRAY_VIEW<const TV> domain,ARRAY_VIEW<const typename T_ARRAYS_T2::ELEMENT> range,const GRID<TV>& grid,T_ARRAYS_T2& grid_data,
        const T_ARRAYS_ARRAY_INT& points_in_cell) const;
    void Instantiate();
//#####################################################################
};
}
#endif
