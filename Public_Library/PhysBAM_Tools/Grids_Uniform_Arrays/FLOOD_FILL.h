//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Eran Guendelman, Frank Losasso, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOOD_FILL
//#####################################################################
#ifndef __FLOOD_FILL__
#define __FLOOD_FILL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <climits>
namespace PhysBAM{

template<int d>
class FLOOD_FILL
{
    typedef VECTOR<int,d> TV_INT;
private:
    STACK<TV_INT> flood_fill_stack;
    bool optimize_fill_for_single_cell_regions;
    TV_INT last_uncolored_node;

public:
    FLOOD_FILL();

    void Optimize_Fill_For_Single_Cell_Regions(const bool optimize=true)
    {optimize_fill_for_single_cell_regions=optimize;}

    int Flood_Fill(ARRAY<int,TV_INT>& colors,const ARRAY<bool,FACE_INDEX<TV_INT::m> >& edge_is_blocked,ARRAY<bool>* color_touches_uncolorable_node=0)
    {return Flood_Fill(colors,edge_is_blocked.data,color_touches_uncolorable_node);}

    // colors should be initialized by the user with -1's where colors will be filled and negative values for nodes which will not be colored.
    // -2 is the distinguished uncolorable node (which color_touches_uncolorable_node refers to)
    int Flood_Fill(ARRAY<int,TV_INT>& colors,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,ARRAY<bool>* color_touches_uncolorable_node=0);
    void Fill_Single_Cell_Regions(ARRAY<int,TV_INT>& colors,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,int& fill_color);
    bool Find_Uncolored_Node(const ARRAY<int,TV_INT>& colors,TV_INT& node_index);
    void Flood_Fill_From_Seed_Node(ARRAY<int,TV_INT>& colors,const int fill_color,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,
        bool& touches_uncolorable_node,const TV_INT& seed_node);
    void Flood_Fill_Node(ARRAY<int,TV_INT>& colors,const int fill_color,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,bool& touches_uncolorable_node,
        STACK<TV_INT>& flood_fill_stack,const TV_INT& node);
    void Identify_Colors_Touching_Boundary(const int number_of_colors,const ARRAY<int,TV_INT>& colors,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,
        ARRAY<bool>& color_touches_boundary);
    void Identify_Colors_Touching_Color(const int color,const int number_of_colors,const ARRAY<int,TV_INT>& colors,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,
        ARRAY<bool>& color_touches_color);
    bool Path_Between_Nodes(const RANGE<TV_INT>& domain,const TV_INT& start_node,const TV_INT& end_node,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,ARRAY<TV_INT>* path=0);
    void Explore_Path(ARRAY<TV_INT,TV_INT>& parents,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,const TV_INT& node);
//#####################################################################
};

template<>
class FLOOD_FILL<1>
{
    typedef VECTOR<int,1> TV_INT;
public:
    ARRAY<VECTOR<int,2> > region_boundaries;

    int Flood_Fill(ARRAYS_ND_BASE<int,TV_INT>& colors,const ARRAY<bool,FACE_INDEX<1> >& edge_is_blocked,ARRAY<bool>* color_touches_uncolorable_node=0)
    {return Flood_Fill(colors,edge_is_blocked.Component(0),color_touches_uncolorable_node);}

    // colors should be initialized by the user with 0's where colors will be filled and negative values for nodes which will not be colored.
    // -1 is the distinguished uncolorable node (which color_touches_uncolorable_node refers to)
    int Flood_Fill(ARRAYS_ND_BASE<int,TV_INT>& colors,const ARRAYS_ND_BASE<bool,TV_INT>& edge_is_blocked_x,ARRAY<bool>* color_touches_uncolorable_node=0);
    void Identify_Colors_Touching_Boundary(const int number_of_colors,const ARRAYS_ND_BASE<int,TV_INT>& colors,const ARRAYS_ND_BASE<bool,TV_INT>& edge_is_blocked_x,ARRAY<bool>& color_touches_boundary);
//#####################################################################
};
}
#endif
