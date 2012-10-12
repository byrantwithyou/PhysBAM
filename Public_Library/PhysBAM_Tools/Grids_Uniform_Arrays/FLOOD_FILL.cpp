//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Eran Guendelman, Frank Losasso, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<int d> FLOOD_FILL<d>::
FLOOD_FILL()
{
    Optimize_Fill_For_Single_Cell_Regions(false);
}
//#####################################################################
// Constructor
//#####################################################################
template<int d> int FLOOD_FILL<d>::
Flood_Fill(ARRAY<int,TV_INT>& colors,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,ARRAY<bool>* color_touches_uncolorable_node)
{
    TV_INT seed_node;
    int fill_color=0;
    if(color_touches_uncolorable_node) color_touches_uncolorable_node->Remove_All();
    if(optimize_fill_for_single_cell_regions){
        Fill_Single_Cell_Regions(colors,edge_is_blocked,fill_color);
        if(color_touches_uncolorable_node) color_touches_uncolorable_node->Resize(fill_color);} // Resize sets new elements to false for us
    flood_fill_stack.Preallocate(colors.counts.Product());
    last_uncolored_node=colors.domain.min_corner;
    while(Find_Uncolored_Node(colors,seed_node)){
        bool touches_uncolorable_node;
        Flood_Fill_From_Seed_Node(colors,fill_color,edge_is_blocked,touches_uncolorable_node,seed_node);
        if(color_touches_uncolorable_node)color_touches_uncolorable_node->Append(touches_uncolorable_node);
        fill_color++;}
    flood_fill_stack.Clean_Memory();
    return fill_color;
}
//#####################################################################
// Function Fill_Single_Cell_Regions
//#####################################################################
template<int d> void FLOOD_FILL<d>::
Fill_Single_Cell_Regions(ARRAY<int,TV_INT>& colors,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,int& fill_color)
{
    for(RANGE_ITERATOR<TV_INT::m> it(colors.domain);it.Valid();it.Next())
        if(colors(it.index)==-1){
            bool ok=true;
            for(int a=0;a<TV_INT::m;a++){
                TV_INT index1(it.index);
                index1(a)++;
                if(!((it.index(a)==colors.domain.min_corner(a) || edge_is_blocked(a)(it.index)) && (it.index(a)==colors.domain.max_corner(a) || edge_is_blocked(a)(index1)))) ok=false;}
            if(ok) colors(it.index)=fill_color++;}
}
//#####################################################################
// Function Find_Uncolored_Node
//#####################################################################
template<int d> bool FLOOD_FILL<d>::
Find_Uncolored_Node(const ARRAY<int,TV_INT>& colors,TV_INT& node_index)
{
    bool first_time_through_loop=true;
    for(RANGE_ITERATOR<TV_INT::m> it(colors.domain);it.Valid();it.Next()){
        if(first_time_through_loop){first_time_through_loop=false;it.index=last_uncolored_node;}
        if(colors(it.index)==-1){last_uncolored_node=node_index=it.index;return true;}}
    return false;
}
//#####################################################################
// Function Flood_Fill_From_Seed_Node
//#####################################################################
template<int d> void FLOOD_FILL<d>::
Flood_Fill_From_Seed_Node(ARRAY<int,TV_INT>& colors,const int fill_color,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,
    bool& touches_uncolorable_node,const TV_INT& seed_node)
{
    assert(colors(seed_node)==-1);
    touches_uncolorable_node=false;
    flood_fill_stack.Remove_All();
    flood_fill_stack.Push(seed_node);
    while(!flood_fill_stack.Empty()){
        TV_INT node=flood_fill_stack.Pop();
        Flood_Fill_Node(colors,fill_color,edge_is_blocked,touches_uncolorable_node,flood_fill_stack,node);}
}
//#####################################################################
// Function Flood_Fill_Node
//#####################################################################
template<int d> void FLOOD_FILL<d>::
Flood_Fill_Node(ARRAY<int,TV_INT>& colors,const int fill_color,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,bool& touches_uncolorable_node,
    STACK<TV_INT>& flood_fill_stack,const TV_INT& node)
{
    if(colors(node)==-2){touches_uncolorable_node=true;return;}
    else if(colors(node)!=-1)return;colors(node)=fill_color;

    for(int a=0;a<TV_INT::m;a++){
        TV_INT node0(node),node1(node);
        node0(a)--;
        node1(a)++;
        if(node(a)>colors.domain.min_corner(a) && !edge_is_blocked(a)(node) && colors(node0)<0) flood_fill_stack.Push(node0);
        if(node(a)<colors.domain.max_corner(a)-1 && !edge_is_blocked(a)(node1) && colors(node1)<0) flood_fill_stack.Push(node1);}
}
//#####################################################################
// Function Identify_Colors_Touching_Boundary
//#####################################################################
template<int d> void FLOOD_FILL<d>::
Identify_Colors_Touching_Boundary(const int number_of_colors,const ARRAY<int,TV_INT>& colors,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,
    ARRAY<bool>& color_touches_boundary)
{
    color_touches_boundary.Resize(number_of_colors);
    color_touches_boundary.Fill(false);
    for(int a=0;a<TV_INT::m;a++)
        for(RANGE_ITERATOR<TV_INT::m-1> it(colors.domain.Remove_Dimension(a));it.Valid();it.Next()){
            int lo=colors(it.index.Insert(colors.domain.min_corner(a),a)),hi=colors(it.index.Insert(colors.domain.max_corner(a),a));
            if(lo>=0) color_touches_boundary(lo)=true;
            if(hi>=0) color_touches_boundary(hi)=true;}
}
//#####################################################################
// Function Identify_Colors_Touching_Color
//#####################################################################
template<int d> void FLOOD_FILL<d>::
Identify_Colors_Touching_Color(const int color,const int number_of_colors,const ARRAY<int,TV_INT>& colors,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,
    ARRAY<bool>& color_touches_color)
{
    color_touches_color.Resize(number_of_colors);
    color_touches_color.Fill(false);
    for(RANGE_ITERATOR<TV_INT::m> it(colors.domain);it.Valid();it.Next())
        if(colors(it.index)==color)
            for(int a=0;a<TV_INT::m;a++){
                TV_INT index0(it.index),index1(it.index);
                index0(a)--;
                index1(a)++;
                if(it.index(a)>colors.domain.min_corner(a) && !edge_is_blocked(a)(it.index) && colors(index0)>=0) color_touches_color(colors(index0))=true;
                if(it.index(a)<colors.domain.max_corner(a)-1 && !edge_is_blocked(a)(index1) && colors(index1)>=0) color_touches_color(colors(index1))=true;}
}
//#####################################################################
// Function Path_Between_Nodes
//#####################################################################
template<int d> bool FLOOD_FILL<d>::
Path_Between_Nodes(const RANGE<TV_INT>& domain,const TV_INT& start_node,const TV_INT& end_node,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,ARRAY<TV_INT>* path)
{
    ARRAY<TV_INT,TV_INT> parents(domain,false);
    parents.Fill(TV_INT()+INT_MAX);
    flood_fill_stack.Remove_All();
    flood_fill_stack.Preallocate(parents.counts.Product());
    flood_fill_stack.Push(end_node);
    bool success=false;
    parents(end_node)=end_node;
    while(!flood_fill_stack.Empty()){
        TV_INT node=flood_fill_stack.Pop();
        if(node==start_node){success=true;break;}
        Explore_Path(parents,edge_is_blocked,node);}
    if(success){
        if(path){path->Remove_All();path->Preallocate(20);
            TV_INT node=start_node;
            for(;;){path->Append(node);node=parents(node);if(node==end_node){path->Append(node);break;}}}
        return true;}
    else return false;
}
//#####################################################################
// Function Explore_Path
//#####################################################################
template<int d> void FLOOD_FILL<d>::
Explore_Path(ARRAY<TV_INT,TV_INT>& parents,const VECTOR<ARRAY_VIEW<bool,TV_INT>,TV_INT::m>& edge_is_blocked,const TV_INT& node)
{
    for(int a=0;a<TV_INT::m;a++){
        TV_INT node0(node),node1(node);
        node0(a)--;
        node1(a)++;
        if(node(a)>parents.domain.min_corner(a) && !edge_is_blocked(a)(node) && parents(node0)(a)==INT_MAX){
            parents(node0)=node;
            flood_fill_stack.Push(node0);}
        if(node(a)<parents.domain.max_corner(a)-1 && !edge_is_blocked(a)(node1) && parents(node0)(a)==INT_MAX){
            parents(node1)=node;
            flood_fill_stack.Push(node1);}}
}
//#####################################################################
// Constructor
//#####################################################################
int FLOOD_FILL<1>::
Flood_Fill(ARRAYS_ND_BASE<int,TV_INT>& colors,const ARRAYS_ND_BASE<bool,TV_INT>& edge_is_blocked_x,ARRAY<bool>* color_touches_uncolorable_node)
{
    int fill_color=0,number_of_regions=0;if(color_touches_uncolorable_node) color_touches_uncolorable_node->Remove_All();region_boundaries.Clean_Memory();
    bool touches_uncolorable_node=false;
    for(int i=colors.domain.min_corner.x;i<colors.domain.max_corner.x;i++){
        if(colors(i)!=-2){
            colors(i)=fill_color;
            if(number_of_regions!=fill_color+1){
                if(color_touches_uncolorable_node) touches_uncolorable_node=i>colors.domain.min_corner.x&&colors(i-1)==-2&&!edge_is_blocked_x(i);
                region_boundaries.Append(VECTOR<int,2>(i,i));number_of_regions=fill_color+1;}
            if(edge_is_blocked_x(i+1)||i==colors.domain.max_corner.x-1||colors(i+1)==-2){
                region_boundaries(fill_color).y=i;fill_color++;
                if(color_touches_uncolorable_node)
                    color_touches_uncolorable_node->Append(touches_uncolorable_node||(i+1<colors.domain.max_corner.x&&colors(i+1)==-2&&!edge_is_blocked_x(i+1)));}}}
    return region_boundaries.m;
}
//#####################################################################
// Function Identify_Colors_Touching_Boundary
//#####################################################################
void FLOOD_FILL<1>::
Identify_Colors_Touching_Boundary(const int number_of_colors,const ARRAYS_ND_BASE<int,TV_INT>& colors,const ARRAYS_ND_BASE<bool,TV_INT>& edge_is_blocked_x,ARRAY<bool>& color_touches_boundary)
{
    color_touches_boundary.Resize(number_of_colors);color_touches_boundary.Fill(false);
    int left_color=colors(colors.domain.min_corner.x),right_color=colors(colors.domain.max_corner.x);
    if(left_color>0)color_touches_boundary(left_color)=true;
    if(right_color>0)color_touches_boundary(right_color)=true;
}
template class FLOOD_FILL<1>;
template class FLOOD_FILL<2>;
template class FLOOD_FILL<3>;
