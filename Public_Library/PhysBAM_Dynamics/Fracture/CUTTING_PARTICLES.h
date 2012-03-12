//#####################################################################
// Copyright 2006, Kevin Der.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUTTING_PARTICLES
//##################################################################### 
#ifndef __CUTTING_PARTICLES__
#define __CUTTING_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
namespace PhysBAM{

class CUTTING_PARTICLES
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    enum CUTTING_PARTICLE_ID_TYPE {TET_NODE_ID,INTERSECTION_ID,TET_NODE_AND_INTERSECTION_ID};
    ARRAY<int> tet_node_indices;
    ARRAY<int> intersection_indices;
    ARRAY<CUTTING_PARTICLE_ID_TYPE> particle_ids_types;
    HASHTABLE<int,int> intersection_to_particle_id;
    HASHTABLE<int,int> tet_node_to_particle_id;

    CUTTING_PARTICLES();
    ~CUTTING_PARTICLES();

    void Preallocate(const int size);
    void Add_Tet_Node_Id(const int tet_node_index);
    void Add_Intersection_Id(const int intersection_index);
    void Add_Tet_Node_And_Intersection_Id(const int tet_node_index,const int intersection_index);
    int Particle_Id_From_Tet_Node(const int tet_node) const;
    int Particle_Id_From_Intersection(const int intersection) const;

    int Number() const
    {return particle_ids_types.m;}

    void Print() const;

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,tet_node_indices,intersection_indices,particle_ids_types,intersection_to_particle_id,tet_node_to_particle_id);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,tet_node_indices,intersection_indices,particle_ids_types,intersection_to_particle_id,tet_node_to_particle_id);}
};
}
#endif
