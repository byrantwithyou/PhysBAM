//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CELL_DOMAIN_INTERFACE_COLOR
//#####################################################################
#ifndef __CELL_DOMAIN_INTERFACE_COLOR__
#define __CELL_DOMAIN_INTERFACE_COLOR__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/LINE_2D.h>
#include <Geometry/Basic_Geometry/PLANE.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>

namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class CELL_DOMAIN_INTERFACE_COLOR:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    TV_INT a;
    int b;

    ARRAY<int> cell_location; // [1] - outside, [0] - boundary (inside, within padding from boundary), [-1] - inside (strictly)

public:

    typedef typename MARCHING_CUBES_COLOR<TV>::CELL_ELEMENTS CELL_ELEMENTS;
    typedef typename MARCHING_CUBES_COLOR<TV>::INTERFACE_ELEMENT INTERFACE_ELEMENT;
    typedef typename MARCHING_CUBES_COLOR<TV>::BOUNDARY_ELEMENT BOUNDARY_ELEMENT;
    typedef typename MARCHING_CUBES_COLOR<TV>::T_SURFACE T_SURFACE;
    typedef typename MARCHING_CUBES_COLOR<TV>::T_FACE T_FACE;

    struct BOUNDARY_CONDITIONS
    {
        enum WORKAROUND{SLIP=-3,DIRICHLET=-2,NEUMANN=-1};
    };

    const GRID<TV>& grid;
    int padding;
    TV_INT size;
    int flat_size;
    int colors;

    bool nc_present;
    bool dc_present;
    bool sc_present;

    ARRAY<int> remap; // maps ghost cells inside in case of wrapping, identity otherwise

    int constraint_base_n;
    int constraint_base_t;
    VECTOR<int*,TV::m> constraint_base;
    int total_number_of_surface_constraints;
    ARRAY<int> flat_base_n;
    ARRAY<int> flat_base_t;
    VECTOR<ARRAY<int>*,TV::m> flat_base;

    int constraint_base_scalar;
    ARRAY<int> flat_base_scalar;

    HASHTABLE<TV_INT,CELL_ELEMENTS> index_to_cell_elements;
    
    CELL_DOMAIN_INTERFACE_COLOR(const GRID<TV>& grid_input,int padding_input,int colors_input);

    int Flatten(const TV_INT& index) const
    {return index.Dot(a)+b;}

    int Flatten_Diff(const TV_INT& index) const
    {return index.Dot(a);}

    bool Is_Outside_Cell(int i) const
    {return cell_location(i)==1;}

    bool Is_Boundary_Cell(int i) const
    {return cell_location(i)==0;}

    bool Is_Boundary_Constraint(int i,int orientation) const
    {return Is_Boundary_Cell((*flat_base(orientation))(i));}

    bool Is_Boundary_Constraint_Scalar(int i) const
    {return Is_Boundary_Cell(flat_base_scalar(i));}

    void Set_Flat_Base_And_Resize(int extra_constraints_n,int extra_constraints_t,const TV_INT& index);
    void Set_Flat_Base_And_Resize_Scalar(int extra_constraints_scalar,const TV_INT& index);
    void Update_Constraint_Count();
    void Update_Total_Constraint_Count();
    void Construct_Surface_Meshes(const GRID<TV>& phi_grid,const ARRAY<T,TV_INT>& phi_value,const ARRAY<int,TV_INT>& phi_color);

    static void Interpolate_Level_Set_To_Double_Fine_Grid(const RANGE<TV_INT>& range_input,const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input,
        const RANGE<TV_INT>& range,ARRAY<T,TV_INT>& phi_value,ARRAY<int,TV_INT>& phi_color,T threshold);
    static void Interpolate_Level_Set_To_Double_Fine_Grid(const GRID<TV>& phi_grid_input,const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input,
        const GRID<TV>& phi_grid,ARRAY<T,TV_INT>& phi_value,ARRAY<int,TV_INT>& phi_color,T tol);
    static void Interpolate_Mac_Level_Set_To_Double_Fine_Grid(const GRID<TV>& phi_grid_input,const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input,
        const GRID<TV>& phi_grid,ARRAY<T,TV_INT>& phi_value,ARRAY<int,TV_INT>& phi_color,T tol);
    static void Interpolate_Level_Set_To_Double_Fine_Grid(const RANGE<TV_INT>& range_input,const ARRAY<T,TV_INT>& phi_value_input,
        const RANGE<TV_INT>& range,ARRAY<T,TV_INT>& phi_value,T threshold);
    static void Interpolate_Level_Set_To_Double_Fine_Grid(const GRID<TV>& phi_grid_input,const ARRAY<T,TV_INT>& phi_value_input,
        const GRID<TV>& phi_grid,ARRAY<T,TV_INT>& phi_value,T tol);
    static void Interpolate_Mac_Level_Set_To_Double_Fine_Grid(const GRID<TV>& phi_grid_input,const ARRAY<T,TV_INT>& phi_value_input,
        const GRID<TV>& phi_grid,ARRAY<T,TV_INT>& phi_value,T tol);

};
}
#endif
