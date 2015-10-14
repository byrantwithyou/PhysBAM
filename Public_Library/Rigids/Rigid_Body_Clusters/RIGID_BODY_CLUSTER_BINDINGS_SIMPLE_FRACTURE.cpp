//#####################################################################
// Copyright 2008, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>::
RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,
    RIGID_BODY_CLUSTER_BINDINGS<TV>& bindings)
    :rigid_body_collection(rigid_body_collection),bindings(bindings),allowed_strain((T)FLT_MAX)
{}
//#####################################################################
// Function Initialize_Strain
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>::
Initialize_Strain(const int parent,FRACTURE_DATA& data)
{
    typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*bindings.reverse_bindings.Get(parent);
    // compute restlengths  
    for(RIGID_CLUSTER_CONSTITUENT_ID i(0);i<cluster.children.Size();i++) for(RIGID_CLUSTER_CONSTITUENT_ID j(0);j<i;j++){
            RIGID_BODY<TV> &child_1=rigid_body_collection.Rigid_Body(cluster.children(i)),&child_2=rigid_body_collection.Rigid_Body(cluster.children(j));
            data.connections.Append(VECTOR<RIGID_CLUSTER_CONSTITUENT_ID,2>(i,j));
            data.restlengths.Append((child_2.Frame().t-child_1.Frame().t).Magnitude());
            LOG::cout<<"restlength on (local coords "<<i<<","<<j<<")  "<<cluster.children(i)<<","<<cluster.children(j)<<" is "<<data.restlengths(data.restlengths.Size())<<std::endl;}
}
//#####################################################################
// Function Find_Weakest_Links
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>::
Find_Weakest_Links(int root,T min_strain,HASHTABLE<int>& visited,ARRAY<int>& edges)
{
    if(!visited.Contains(root)) visited.Insert(root);
    ARRAY<int> next_nodes;
    for(int i=0;i<graph->Adjacent_Edges(root).m;i++){
        PAIR<int,int> edge=graph->Edges(graph->Adjacent_Edges(root)(i));
        int other_node=edge.x;if(other_node==root) other_node=edge.y;
        if(visited.Contains(other_node)) continue;
        T strain=allowed_strains.Get_Default(VECTOR<int,2>(root,other_node).Sorted(),allowed_strain);
        if(strain<=min_strain){
            if(strain<min_strain){min_strain=strain;edges.Remove_All();next_nodes.Remove_All();}
            edges.Append(graph->Adjacent_Edges(root)(i));
            next_nodes.Append(other_node);}}
    for(int i=0;i<next_nodes.m;i++) Find_Weakest_Links(next_nodes(i),min_strain,visited,edges);
}
//#####################################################################
// Function Compute_New_Clusters_Based_On_Unclustered_Strain
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>::
Compute_New_Clusters_Based_On_Unclustered_Strain()
{
    parents_to_rebuild.Remove_All();
    ARRAY<int,int> components;if(graph) graph->Connected_Components(components);
    for(typename HASHTABLE<int,FRACTURE_DATA>::ITERATOR iterator(fracture_data);iterator.Valid();iterator.Next()){
        ARRAY<int> remove_connections;HASHTABLE<int> visited;bool need_rebuild=false;
        FRACTURE_DATA& data=iterator.Data();
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*bindings.reverse_bindings.Get(iterator.Key());
        for(int i=data.connections.m-1;i>=0;i--){
            const VECTOR<RIGID_CLUSTER_CONSTITUENT_ID,2>& edge=data.connections(i);
            T rl=data.restlengths(i);
            RIGID_BODY<TV> &child_1=rigid_body_collection.Rigid_Body(cluster.children(edge[0])),&child_2=rigid_body_collection.Rigid_Body(cluster.children(edge[1]));
            VECTOR<int,2> hash_index=VECTOR<int,2>(cluster.children(edge[0]),cluster.children(edge[1])).Sorted();
            T strain=abs((child_2.Frame().t-child_1.Frame().t).Magnitude()/rl-(T)1);
            T& local_decay=decay.Get_Or_Insert(hash_index,0);local_decay+=local_dt*decay_rate.Get_Or_Insert(hash_index,0);
            strain+=local_decay;
            LOG::cout<<"strain between "<<cluster.children(edge[0])<<","<<cluster.children(edge[1])<<" is "<<strain<<std::endl;
            T allowed_strain_local=allowed_strains.Get_Default(VECTOR<int,2>(cluster.children(edge[0]),cluster.children(edge[1])).Sorted(),allowed_strain);
            if(strain>allowed_strain_local){need_rebuild=true;
                if(graph){
                    ARRAY<int> break_connections;int root=cluster.children(edge[0]);
                    if(!visited.Contains(root)){
                        Find_Weakest_Links(root,allowed_strain_local,visited,break_connections);
                        for(int i=0;i<break_connections.m;i++){
                            VECTOR<int,2> edge=VECTOR<int,2>(graph->Edges(break_connections(i)).x,graph->Edges(break_connections(i)).y).Sorted();
                            graph->Remove_Edge(break_connections(i));
                            for(int j=data.connections.m-1;j>=0;j--){
                                VECTOR<int,2> tmp_edge(Value(data.connections(j)(0)),Value(data.connections(j)(1)));
                                if(edge==tmp_edge.Sorted()) remove_connections.Append(j);}}}}
                else data.connections.Remove_Index_Lazy(i);}}
        if(need_rebuild) parents_to_rebuild.Append(iterator.Key());
        remove_connections.Sort();
        for(int i=remove_connections.m-1;i>=0;i--) data.connections.Remove_Index_Lazy(remove_connections(i));}
}
//#####################################################################
// Function Create_New_Clusters
//#####################################################################
template<class TV> bool RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<TV>::
Create_New_Clusters()
{
    for(int i_dummy=0;i_dummy<parents_to_rebuild.m;i_dummy++){int parent=parents_to_rebuild(i_dummy);
        FRACTURE_DATA& data=fracture_data.Get(parent);
        typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster=*bindings.reverse_bindings.Get(parent);
        UNDIRECTED_GRAPH<RIGID_CLUSTER_CONSTITUENT_ID,int> undirected(cluster.children.Size());
        for(int i=0;i<data.connections.m;i++) undirected.Add_Edge(data.connections(i).x,data.connections(i).y,undirected.Last_Edge()+1);
        ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> components;
        int number_components=undirected.Connected_Components(components);
        if(number_components==1) continue;
        ARRAY<ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> > new_cluster_constituents(number_components);
        for(RIGID_CLUSTER_CONSTITUENT_ID i(0);i<cluster.children.Size();i++)
            new_cluster_constituents(components(i)).Append(cluster.children(i));
        bindings.Delete_Binding(parent);
        for(int i=0;i<new_cluster_constituents.m;i++) if(new_cluster_constituents(i).Size()>RIGID_CLUSTER_CONSTITUENT_ID(1)) bindings.Add_Binding(new_cluster_constituents(i));}
    return parents_to_rebuild.m>0;
}
//#####################################################################
template class RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<VECTOR<float,2> >;
template class RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE<VECTOR<double,2> >;
}
