//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Neil Molino, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYPOTHETICAL_CUT
//##################################################################### 
#ifndef __HYPOTHETICAL_CUT__
#define __HYPOTHETICAL_CUT__

#include <Core/Utilities/NONCOPYABLE.h>
#include <Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <Deformables/Fracture/EMBEDDED_TRIANGULATED_OBJECT.h>
#include <Deformables/Fracture/EMBEDDING_POLICY.h>
#include <Solids/Fracture/HYPOTHETICAL_NODE.h>
namespace PhysBAM{

template<class TV,int d>
class HYPOTHETICAL_CUT:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename EMBEDDING_POLICY<TV,d>::EMBEDDED_OBJECT T_EMBEDDED_OBJECT;
protected:
    const T_EMBEDDED_OBJECT& embedded_object;
    T cut_quality_metric;
public:
    ARRAY<HYPOTHETICAL_NODE<TV,d> > hypothetical_nodes;
    
    HYPOTHETICAL_CUT(const T_EMBEDDED_OBJECT& embedded_object_input)
        :embedded_object(embedded_object_input),cut_quality_metric(0)
    {}

protected:
    void Add_Hypothetical_Node(const int parent1,const int parent2,const T interpolation_fraction)
    {hypothetical_nodes.Append(HYPOTHETICAL_NODE<TV,d>(embedded_object,parent1,parent2,interpolation_fraction));}
public:

//##################################################################### 
    bool Contains_Embedded_Particle(const EMBEDDED_OBJECT<TV,d>& embedded_object,const int emb_particle) const;
protected:
    TV Position(const int hypothetical_node_index) const;
public:
    void Add_Hypothetical_Nodes_To_Embedded_Object(EMBEDDED_OBJECT<TV,d>& embedded_object);
//##################################################################### 
};
}
#endif
