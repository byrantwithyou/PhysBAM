//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HIDDEN_SURFACE
//#####################################################################
#ifndef __HIDDEN_SURFACE__
#define __HIDDEN_SURFACE__

#ifdef USE_BOOST_GEOMETRY
#include <Tools/Data_Structures/DIRECTED_GRAPH.h>
#include <Geometry/Images/HIDDEN_SURFACE_PRIMITIVES.h>

namespace PhysBAM{

// Arrows point away from eye.
template<class T>
class HIDDEN_SURFACE
{
public:
    HIDDEN_SURFACE_PRIMITIVES<T>& primitives;

    HIDDEN_SURFACE(HIDDEN_SURFACE_PRIMITIVES<T>& hsp)
        :primitives(hsp)
    {}

    void Compute();
    void Break_Graph(DIRECTED_GRAPH<>& graph,ARRAY<int>& node_map);
    void Break_Component(DIRECTED_GRAPH<>& graph,ARRAY<int>& node_map);
};

}
#endif
#endif
