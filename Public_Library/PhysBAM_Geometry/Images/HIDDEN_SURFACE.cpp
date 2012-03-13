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
    ARRAY<int> components;
    int num_components=graph.Strongly_Connected_Components(components);
    ARRAY<ARRAY<int> > component_nodes(num_components);
    for(int i=0;i<components.m;i++) component_nodes(components(i)).Append(i);
    for(int i=0;i<component_nodes.m;i++){
        if(component_nodes(i).m==1){
            primitives.Emit_Node(node_map(component_nodes(i)(0)));
            continue;}
        const ARRAY<int>& nodes=component_nodes(i);
        DIRECTED_GRAPH<> component(nodes.m);
        ARRAY<int> component_node_map(nodes.m);
        for(int j=0;j<nodes.m;j++){
            component_node_map(j)=node_map(nodes(j));
            const ARRAY<int>& children=graph.Children(nodes(j));
            for(int k=0;k<children.m;k++)
                if(components(children(k))==i)
                    component.Add_Edge(nodes(j),children(k));}
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
        ARRAY<bool> is_child(node_map.m);
        is_child.Subset(children).Fill(true);

        // Divide the nodes; update node_map.
        ARRAY<ARRAY<int> > inside_list(children.m),outside_list(children.m);
        for(int i=0;i<children.m;i++){
            primitives.Divide_Primitive(children(i),index,inside_list(i),outside_list(i));
            node_map(children(i))=inside_list(i)(0);
            inside_list(i)(0)=children(i);
            for(int j=1;j<inside_list(i).m;j++){
                int k=node_map.Append(inside_list(i)(j));
                inside_list(i)(j)=k;}
            for(int j=0;j<outside_list(i).m;j++){
                int k=node_map.Append(outside_list(i)(j));
                outside_list(i)(j)=k;}}

        new_graph.Initialize(node_map.m);

        // Node types: normal, index, inside, outside
        // add (index, inside) edges
        for(int i=0;i<inside_list.m;i++)
            for(int j=0;j<inside_list(i).m;j++)
                new_graph.Add_Edge(index,inside_list(i)(j));

        // add (normal, normal) and (normal,index) edges
        for(int i=0;i<num_nodes;i++){
            if(i==index || is_child(i)) continue;
            const ARRAY<int>& ch=graph.Children(i);
            for(int j=0;j<ch.m;j++)
                if(!is_child(ch(j)))
                    new_graph.Add_Edge(i,ch(j));}

        for(int i=0;i<children.m;i++){
            int p=children(i);
            const ARRAY<int>& ch=graph.Children(p),&pin=inside_list(p),&pout=outside_list(p);
            for(int j=0;j<ch.m;j++)
                if(is_child(ch(j))){
                    // Add (inside, inside) edges
                    const ARRAY<int>& cin=inside_list(ch(j));
                    for(int k=0;k<pin.m;k++)
                        for(int m=0;m<cin.m;m++)
                            if(primitives.Test_Edge(pin(k),cin(m)))
                                new_graph.Add_Edge(pin(k),cin(m));
                    // Add (outside, outside) edges
                    const ARRAY<int>& cout=outside_list(ch(j));
                    for(int k=0;k<pout.m;k++)
                        for(int m=0;m<cout.m;m++)
                            if(primitives.Test_Edge(pout(k),cout(m)))
                                new_graph.Add_Edge(pout(k),cout(m));}
                else{ // Add (outside, normal) edges
                    int m=ch(j);
                    for(int k=0;k<pout.m;k++)
                        if(primitives.Test_Edge(pout(k),m))
                            new_graph.Add_Edge(pout(k),m);}}

        // Add (normal, outside) edges
        for(int i=0;i<children.m;i++){
            int c=children(i);
            const ARRAY<int>& pa=graph.Parents(c),&cout=outside_list(c);
            for(int j=0;j<pa.m;j++){
                int m=pa(j);
                if(!is_child(m) && m!=index)
                    for(int k=0;k<cout.m;k++)
                        if(primitives.Test_Edge(m,cout(k)))
                            new_graph.Add_Edge(m,cout(k));}}

        graph.Reset();
    }
    Break_Graph(new_graph,node_map);
}
template class HIDDEN_SURFACE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HIDDEN_SURFACE<double>;
#endif
