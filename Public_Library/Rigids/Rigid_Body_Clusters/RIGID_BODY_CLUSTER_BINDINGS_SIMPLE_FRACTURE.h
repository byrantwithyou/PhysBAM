//#####################################################################
// Copyright 2008, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE__
#define __RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE__
#include <Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <Tools/Data_Structures/UNDIRECTED_GRAPH.h>
#include <Tools/Log/LOG.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>

namespace PhysBAM{

template<class TV>
class RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE:public RIGID_BODY_CLUSTER_BINDINGS<TV>::CALLBACKS
{
    typedef typename TV::SCALAR T;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& bindings;
    ARRAY<int> parents_to_rebuild;
public:
    struct FRACTURE_DATA
    {
        ARRAY<VECTOR<RIGID_CLUSTER_CONSTITUENT_ID,2> > connections;
        ARRAY<T> restlengths;
    };

    T allowed_strain;
    HASHTABLE<VECTOR<int,2>,T> allowed_strains,decay,decay_rate;
    UNDIRECTED_GRAPH<int,int>* graph;
    T local_dt;

    RIGID_BODY_CLUSTER_BINDINGS_SIMPLE_FRACTURE(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,
        RIGID_BODY_CLUSTER_BINDINGS<TV>& bindings);

    HASHTABLE<int,FRACTURE_DATA> fracture_data; // from cluster parent to FRACTURE_DATA

    void Pre_Advance_Unclustered(const T dt,const T time){}
    void Post_Advance_Unclustered(const T dt,const T time){local_dt=dt;}

    void Create_Cluster(const int parent)
    {Initialize_Strain(parent,fracture_data.Get_Or_Insert(parent));}

    void Destroy_Cluster(const int parent)
    {fracture_data.Delete(parent);}

    void Initialize_Strain(const int parent,FRACTURE_DATA& data);
    void Find_Weakest_Links(int root,T min_strain,HASHTABLE<int>& visited,ARRAY<int>& edges);
    void Compute_New_Clusters_Based_On_Unclustered_Strain();
    bool Create_New_Clusters();
//#####################################################################
};
}
#endif
