//#####################################################################
// Copyright 2004-2008, Ronald Fedkiw, Geoffrey Irving, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNDIRECTED_GRAPH
//#####################################################################
#include <Core/Data_Structures/DIRECTED_GRAPH.h>
#include <Core/Data_Structures/ELEMENT_ID.h>
#include <Core/Data_Structures/QUEUE.h>
#include <Core/Data_Structures/STACK.h>
#include <Core/Data_Structures/UNDIRECTED_GRAPH.h>
#include <Core/Math_Tools/exchange.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
UNDIRECTED_GRAPH_CORE::
UNDIRECTED_GRAPH_CORE(const int number_nodes,const int number_edges)
    :edges(number_edges),adjacent_edges(number_nodes)
{}
//#####################################################################
// Function Reset
//#####################################################################
void UNDIRECTED_GRAPH_CORE::
Reset()
{
    edges.Remove_All();
    for(int i=0;i<adjacent_edges.m;i++) adjacent_edges(i).Remove_All();
    adjacent_edges.Remove_All();
}
//#####################################################################
// Function Add_Edge
//#####################################################################
void UNDIRECTED_GRAPH_CORE::
Add_Edge(const int parent,const int child,const int edge)
{
    if(edges.m<=edge) edges.Resize(edge+1);
    if(adjacent_edges.m<=max(parent,child)) adjacent_edges.Resize(max(parent,child)+1);
    edges(edge).x=parent;edges(edge).y=child;
    if(parent>=0) adjacent_edges(parent).Append(edge);
    if(child>=0) adjacent_edges(child).Append(edge);
}
//#####################################################################
// Function Remove_Edge
//#####################################################################
void UNDIRECTED_GRAPH_CORE::
Remove_Edge(const int edge)
{
    ARRAY<int> &adjacent_parent=adjacent_edges(edges(edge).x),&adjacent_child=adjacent_edges(edges(edge).y);
    adjacent_parent.Remove_Index_Lazy(adjacent_parent.Find(edge));adjacent_child.Remove_Index_Lazy(adjacent_child.Find(edge));
    edges(edge)=PAIR<int,int>(-1,-1);
}
//#####################################################################
// Function Modify_Edge
//#####################################################################
void UNDIRECTED_GRAPH_CORE::
Modify_Edge(const int parent,const int child,const int edge)
{
    if(parent==edges(edge).y || child==edges(edge).x) exchange(edges(edge).x,edges(edge).y);
    if(parent!=edges(edge).x){
        ARRAY<int> &adjacent=adjacent_edges(edges(edge).x);
        adjacent.Remove_Index_Lazy(adjacent.Find(edge));
        adjacent_edges(parent).Append(edge);
        edges(edge).x=parent;}
    if(child!=edges(edge).y){
        ARRAY<int> &adjacent=adjacent_edges(edges(edge).y);
        adjacent.Remove_Index_Lazy(adjacent.Find(edge));
        adjacent_edges(child).Append(edge);
        edges(edge).y=child;}
}
//#####################################################################
// Function Breadth_First_Directed_Graph
//#####################################################################
void UNDIRECTED_GRAPH_CORE::
Breadth_First_Directed_Graph(const int root_node,DIRECTED_GRAPH_CORE& directed_graph)
{
    directed_graph.Initialize(adjacent_edges.m);ARRAY<bool,int> marked_nodes(adjacent_edges.m);ARRAY<bool,int> marked_edges(edges.m);
    QUEUE<int> queue(Value(adjacent_edges.m)+1);
    marked_nodes(root_node)=true;queue.Enqueue(root_node);
    Generate_Next_Level_Of_Breadth_First_Directed_Graph(directed_graph,marked_nodes,marked_edges,queue);
}
//#####################################################################
// Function Connected_Components
//#####################################################################
int UNDIRECTED_GRAPH_CORE::
Connected_Components(ARRAY<int,int>& component_id) // returns total number of components
{
    component_id.Resize(adjacent_edges.Size());
    int component=0;
    for(int seed=0;seed<component_id.Size();seed++) if(component_id(seed)<0){
        STACK<int> stack;
        stack.Push(seed);
        component_id(seed)=++component;
        while(!stack.Empty()){
            int current=stack.Pop();
            for(int j=0;j<adjacent_edges(current).m;j++){
                int eid=adjacent_edges(current)(j);
                int other_node=(edges(eid).x==current?edges(eid).y:edges(eid).x);
                if(component_id(other_node)<0){
                    stack.Push(other_node);
                    component_id(other_node)=component;}}}}
    return component;
}
//#####################################################################
// Function Generate_Next_Level_Of_Breadth_First_Directed_Graph
//#####################################################################
void UNDIRECTED_GRAPH_CORE::
Generate_Next_Level_Of_Breadth_First_Directed_Graph(DIRECTED_GRAPH_CORE& directed_graph,ARRAY<bool,int>& marked_nodes,ARRAY<bool,int>& marked_edges,QUEUE<int>& queue)
{
    const int node=queue.Dequeue();
    for(int i=0;i<adjacent_edges(node).m;i++){int edge=adjacent_edges(node)(i);if(!marked_edges(edge)){
        int parent=edges(edge).x,child=edges(edge).y;
        int next_node=parent==node?child:parent;
        if(next_node>=0){
            directed_graph.Add_Edge(node,next_node);marked_edges(edge)=true;
            if(!marked_nodes(next_node)){marked_nodes(next_node)=true;queue.Enqueue(next_node);}}}}
    if(!queue.Empty()) Generate_Next_Level_Of_Breadth_First_Directed_Graph(directed_graph,marked_nodes,marked_edges,queue);
}
//#####################################################################
// Function Depth_First_Directed_Graph
//#####################################################################
void UNDIRECTED_GRAPH_CORE::
Depth_First_Directed_Graph(const int root_node,DIRECTED_GRAPH_CORE& directed_graph)
{
    directed_graph.Initialize(adjacent_edges.m);ARRAY<bool,int> marked(adjacent_edges.m);
    Generate_Next_Level_Of_Depth_First_Directed_Graph(directed_graph,marked,root_node);
}
//#####################################################################
// Function Generate_Next_Level_Of_Depth_First_Directed_Graph
//#####################################################################
void UNDIRECTED_GRAPH_CORE::
Generate_Next_Level_Of_Depth_First_Directed_Graph(DIRECTED_GRAPH_CORE& directed_graph,ARRAY<bool,int>& marked,const int node)
{
    for(int i=0;i<adjacent_edges(node).m;i++){int edge=adjacent_edges(node)(i);
        int parent=edges(edge).x,child=edges(edge).y;
        int next_node=parent==node?child:parent;
        if(next_node>=0){
            directed_graph.Add_Edge(node,next_node);
            if(!marked(next_node)){marked(next_node)=true;Generate_Next_Level_Of_Depth_First_Directed_Graph(directed_graph,marked,next_node);}}}
}
//#####################################################################
// Function Read
//#####################################################################
template<class RW> void UNDIRECTED_GRAPH_CORE::
Read(std::istream& input)
{
    ARRAY<int,int> parent_node_index,child_node_index;
    Read_Binary<RW>(input,parent_node_index,child_node_index,adjacent_edges);
    if(parent_node_index.Size()!=child_node_index.Size()) PHYSBAM_FATAL_ERROR();
    edges.Resize(parent_node_index.Size(),no_init);
    for(int e=0;e<edges.m;e++){edges(e).x=parent_node_index(e);edges(e).y=child_node_index(e);}
}
//#####################################################################
// Function Write
//#####################################################################
template<class RW> void UNDIRECTED_GRAPH_CORE::
Write(std::ostream& output) const
{
    ARRAY<int,int> parent_node_index(edges.m,no_init),child_node_index(edges.m,no_init);
    for(int e=0;e<edges.m;e++){parent_node_index(e)=edges(e).x;child_node_index(e)=edges(e).y;}
    Write_Binary<RW>(output,parent_node_index,child_node_index,adjacent_edges);
}
//#####################################################################
template void UNDIRECTED_GRAPH_CORE::Read<double>(std::istream&);
template void UNDIRECTED_GRAPH_CORE::Read<float>(std::istream&);
template void UNDIRECTED_GRAPH_CORE::Write<double>(std::ostream&) const;
template void UNDIRECTED_GRAPH_CORE::Write<float>(std::ostream&) const;
}
