//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIPLE_JUNCTION_CORRECTION
//#####################################################################
#ifndef __TRIPLE_JUNCTION_CORRECTION__
#define __TRIPLE_JUNCTION_CORRECTION__
#include <Core/Arrays/ARRAY.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
namespace PhysBAM{
template<class TV>
class TRIPLE_JUNCTION_CORRECTION
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename MARCHING_CUBES_COLOR<TV>::HASH_CELL_DATA HASH_CELL_DATA;
    typedef typename MARCHING_CUBES_COLOR<TV>::CELL_ELEMENTS CELL_ELEMENTS;
    typedef typename MARCHING_CUBES_COLOR<TV>::BOUNDARY_ELEMENT BOUNDARY_ELEMENT;
    typedef typename MARCHING_CUBES_COLOR<TV>::INTERFACE_ELEMENT INTERFACE_ELEMENT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    typedef VECTOR<T,(1<<TV::m)> PHI;
public:
    enum WORKAROUND {num_stencils=TV::m==2?2:5};
    const GRID<TV>& grid;
    ARRAY<ARRAY<T,TV_INT> >& phi;
    int ghost;
    T default_phi;
    T trust_buffer; // distance from sonic points to trusted region
    T valid_width; // distance from interface we care about
    T extent; // compute pairwise level set far enough so we can exprapolate
    int extrap_width; // how far to extrapolate
    int bc_colors;

    struct PAIRWISE_LEVEL_SET_DATA
    {
        T phi;
        int valid_flags;
        VECTOR<short,2> trust;

        PAIRWISE_LEVEL_SET_DATA();
    };

    ARRAY<PAIRWISE_LEVEL_SET_DATA,TV_INT> pairwise_data;
    ARRAY<ARRAY<ARRAY<T,TV_INT> > > pairwise_phi;
    ARRAY<T,TV_INT> combined_phi;
    ARRAY<int,TV_INT> combined_color;

    TRIPLE_JUNCTION_CORRECTION(const GRID<TV>& grid,ARRAY<ARRAY<T,TV_INT> >& phi,int ghost);
    void Compute_Pairwise_Data();
    void Initialize_Pairwise_Level_Set();
    void Fill_Valid_Region_With_Exprapolation();
    void One_Step_Triple_Junction_Correction();
    void Update_Color_Level_Sets();
    void Cut_Interface(HASHTABLE<TV_INT,CELL_ELEMENTS>& index_to_cell_data);
    void Cut_Cell_With_Pairwise_Phi(HASHTABLE<TV_INT,CELL_ELEMENTS>& index_to_cell_data,const TV_INT& cell);
    void Compute_Pairwise_Level_Set_Data();
    int Fill_Combined_Level_Set_At_Index(const TV_INT& node);
    int Fill_Phi_From_Pairwise_Level_Set_At_Index(const TV_INT& node,int color);
};
}
#endif
