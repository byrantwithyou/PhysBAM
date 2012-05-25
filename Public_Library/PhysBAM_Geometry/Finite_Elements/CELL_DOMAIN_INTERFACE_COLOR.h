//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CELL_DOMAIN_INTERFACE_COLOR
//#####################################################################
#ifndef __CELL_DOMAIN_INTERFACE_COLOR__
#define __CELL_DOMAIN_INTERFACE_COLOR__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

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

    struct BOUNDARY_CONDITIONS
    {
        enum WORKAROUND{SLIP=-3,DIRICHLET=-2,NEUMANN=-1};
    };

    const GRID<TV>& grid;
    const int padding;
    const TV_INT size;
    const int flat_size;
    const int colors;
    const bool wrap;

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
    
    CELL_DOMAIN_INTERFACE_COLOR(const GRID<TV>& grid_input,int padding_input,int colors_input,bool wrap_input);

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

    static void Interpolate_Level_Set_To_Double_Fine_Grid(const GRID<TV>& phi_grid_input,
        const ARRAY<T,TV_INT>& phi_value_input,const ARRAY<int,TV_INT>& phi_color_input,
        const GRID<TV>& phi_grid,ARRAY<T,TV_INT>& phi_value,ARRAY<int,TV_INT>& phi_color,T tol=1e-2);
};
}
#endif
