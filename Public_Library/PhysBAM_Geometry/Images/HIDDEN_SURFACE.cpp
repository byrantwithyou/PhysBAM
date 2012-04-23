#ifdef USE_BOOST_GEOMETRY
//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Geometry/Images/HIDDEN_SURFACE.h>
using namespace PhysBAM;
template<class TV> void HIDDEN_SURFACE<TV>::
Compute()
{
    DIRECTED_GRAPH<> graph;
    primitives.Initialize(graph);
    ARRAY<int> node_map(IDENTITY_ARRAY<>(graph.Number_Nodes()));
    Break_Graph(graph,node_map);
}
template<class TV> void HIDDEN_SURFACE<TV>::
Break_Graph(DIRECTED_GRAPH<>& graph,ARRAY<int>& node_map)
{
    ARRAY<int> components,reverse_map(node_map.m);
    int num_components=graph.Strongly_Connected_Components(components);
    ARRAY<ARRAY<int> > component_nodes(num_components);
    for(int i=0;i<components.m;i++) component_nodes(components(i)).Append(i);
    for(int i=component_nodes.m-1;i>=0;i--){
        if(component_nodes(i).m==1){
            primitives.Emit_Node(node_map(component_nodes(i)(0)));
            continue;}
        const ARRAY<int>& nodes=component_nodes(i);
        DIRECTED_GRAPH<> component(nodes.m);
        ARRAY<int> component_node_map(nodes.m);
        reverse_map.Subset(nodes)=IDENTITY_ARRAY<>(nodes.m);
        for(int j=0;j<nodes.m;j++){
            component_node_map(j)=node_map(nodes(j));
            const ARRAY<int>& children=graph.Children(nodes(j));
            for(int k=0;k<children.m;k++)
                if(components(children(k))==i)
                    component.Add_Edge(reverse_map(nodes(j)),reverse_map(children(k)));}
        Break_Component(component,component_node_map);}
}
template<class TV> void HIDDEN_SURFACE<TV>::
Break_Component(DIRECTED_GRAPH<>& graph,ARRAY<int>& node_map)
{
    DIRECTED_GRAPH<> new_graph;

    // Make this local to save memory
    {
        // Find minimum degree node
        int minimum_degree=node_map.m+1,index=-1,num_nodes=node_map.m;
        for(int i=0;i<node_map.m;i++)
            if(minimum_degree>graph.Children(i).m){
                minimum_degree=graph.Children(i).m;
                index=i;}

        // Make it easy to determine which nodes are being divided
        const ARRAY<int>& children=graph.Children(index);
        ARRAY<int> child_index(node_map.m);
        child_index.Fill(-1);
        child_index.Subset(children)=IDENTITY_ARRAY<>(children.m);

        // Divide the nodes; update node_map.
        ARRAY<int> inside_list(children.m),outside_list(children.m);
        for(int i=0;i<children.m;i++){
            int added=primitives.Divide_Primitive(node_map(children(i)),node_map(index));
            if(added<0){
                outside_list(i)=-1;
                inside_list(i)=children(i);}
            else if(added==node_map(children(i))){
                inside_list(i)=-1;
                outside_list(i)=children(i);}
            else{
                inside_list(i)=children(i);
                outside_list(i)=node_map.Append(added);}}

        new_graph.Initialize(node_map.m);

        // Node types: normal, index, inside, outside
        // add (index, inside) edges
        for(int i=0;i<inside_list.m;i++)
            if(inside_list(i)>=0)
                new_graph.Add_Edge(index,inside_list(i));

        // add (normal, normal) and (normal,index) edges
        for(int i=0;i<num_nodes;i++){
            if(i==index || child_index(i)>=0) continue;
            const ARRAY<int>& ch=graph.Children(i);
            for(int j=0;j<ch.m;j++)
                if(child_index(ch(j))<0)
                    new_graph.Add_Edge(i,ch(j));}

        for(int i=0;i<children.m;i++){
            int p=children(i);
            const ARRAY<int>& ch=graph.Children(p);
            int pin=inside_list(i),pout=outside_list(i);
            for(int j=0;j<ch.m;j++)
                if(child_index(ch(j))>=0){
                    int child=child_index(ch(j));
                    int cin=inside_list(child);
                    int cout=outside_list(child);
                    // Add (inside, inside) edges
                    if(pin>=0 && cin>=0)
                        if(primitives.Test_Edge(node_map(pin),node_map(cin)))
                            new_graph.Add_Edge(pin,cin);
                    // Add (outside, outside) edges
                    if(pout>=0 && cout>=0)
                        if(primitives.Test_Edge(node_map(pout),node_map(cout)))
                            new_graph.Add_Edge(pout,cout);}
                else{ // Add (outside, normal) edges
                    int m=ch(j);
                    if(pout>=0)
                        if(primitives.Test_Edge(node_map(pout),node_map(m)))
                            new_graph.Add_Edge(pout,m);}}

        // Add (normal, outside) edges
        for(int i=0;i<children.m;i++){
            const ARRAY<int>& pa=graph.Parents(children(i));
            int cout=outside_list(i);
            if(cout<0) continue;
            for(int j=0;j<pa.m;j++){
                int m=pa(j);
                if(child_index(m)<0 && m!=index)
                    new_graph.Add_Edge(m,cout);}}

        graph.Reset();
    }
    Break_Graph(new_graph,node_map);
}
template class HIDDEN_SURFACE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HIDDEN_SURFACE<double>;
#endif
#endif
