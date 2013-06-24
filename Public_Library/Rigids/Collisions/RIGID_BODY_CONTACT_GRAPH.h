//#####################################################################
// Copyright 2003-2008, Ronald Fedkiw, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_CONTACT_GRAPH
//##################################################################### 
#ifndef __RIGID_BODY_CONTACT_GRAPH__
#define __RIGID_BODY_CONTACT_GRAPH__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/DIRECTED_GRAPH.h>
#include <Tools/Data_Structures/ELEMENT_ID.h>
#include <Tools/Log/LOG.h>
#include <Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
namespace PhysBAM{

template<class TV>
class RIGID_BODY_CONTACT_GRAPH
{
public:
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles;
    DIRECTED_GRAPH<int> directed_graph; // nodes point to objects which are "above" it

    RIGID_BODY_CONTACT_GRAPH(RIGID_BODY_PARTICLES<TV>& rigid_body_particles_input) 
        :rigid_body_particles(rigid_body_particles_input),directed_graph(rigid_body_particles.Size())
    {}

    void Reset()
    {directed_graph.Reset();}

    void Initialize() // Call when number of rigid bodies has changed
    {directed_graph.Initialize(rigid_body_particles.Size());}

    void Add_Edge(const int body_below,const int body_above) 
    {directed_graph.Add_Edge(body_below,body_above);}

    int Number_Of_Levels() const
    {return directed_graph.Number_Of_Levels();}

    void Print() const
    {LOG::cout<<"CONTACT_GRAPH:"<<std::endl;
    for(int i=0;i<rigid_body_particles.Size();i++){
        LOG::cout<<"CONTACT_GRAPH\""<<rigid_body_particles.Rigid_Body(i).name<<"\" (LEVEL = "<<directed_graph.Level_Of_Node(i)<<"): ";
        if(directed_graph.Parents(i).m>0){LOG::cout<<" IS ABOVE=";RIGID_BODY<TV>::Print_Names(rigid_body_particles,directed_graph.Parents(i));}
        if(directed_graph.Children(i).m>0){LOG::cout<<" IS BELOW=";RIGID_BODY<TV>::Print_Names(rigid_body_particles,directed_graph.Children(i));}
        LOG::cout<<std::endl;}LOG::cout<<std::endl;
    for(int i=0;i<directed_graph.Number_Of_Levels();i++){
        LOG::cout<<"CONTACT_GRAPH Level "<<i<<": ";RIGID_BODY<TV>::Print_Names(rigid_body_particles,directed_graph.Nodes_In_Level(i));LOG::cout<<std::endl;}LOG::cout<<std::endl;}

    void Print_Statistics(const ARRAY<ARRAY<int> >& contact_pairs_for_level,const bool verbose=false)
    {int max_level_size=0,max_level_index=0,max_pairs_in_level=0,max_pairs_in_level_index=0,number_levels=Number_Of_Levels();int total=0;
    for(int i=0;i<number_levels;i++) {
        int level_size=directed_graph.Nodes_In_Level(i);total+=level_size;
        if(level_size > max_level_size){max_level_size=level_size;max_level_index=i;}
        if(contact_pairs_for_level(i).m > max_pairs_in_level){max_pairs_in_level=contact_pairs_for_level(i).m;max_pairs_in_level_index=i;}}
    assert(total==rigid_body_particles.Size());
    LOG::cout<<"Contact graph statistics: levels="<<number_levels<<", max_level_size="<<max_level_size<<", max_pairs_in_level="<<max_pairs_in_level<<std::endl;
    if(verbose){
        if(max_level_index > 0){
            LOG::cout<<"Bodies in biggest level (level="<<max_level_index<<"): ";
            RIGID_BODY<TV>::Print_Names(rigid_body_particles,directed_graph.Nodes_In_Level(max_level_index));LOG::cout<<std::endl;}
        if(max_pairs_in_level_index > 0){
            LOG::cout<<"Pairs in biggest level (level="<<max_pairs_in_level_index<<"): ";
            RIGID_BODY<TV>::Print_Pairs(rigid_body_particles,contact_pairs_for_level(max_pairs_in_level_index));LOG::cout<<std::endl;}}}

//#####################################################################
};   
}
#endif

