//#####################################################################
// Copyright 2006, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//######################################################################
// Class FREE_PARTICLES
//######################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FREE_PARTICLES<TV>::
FREE_PARTICLES()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FREE_PARTICLES<TV>::
~FREE_PARTICLES()
{}
//#####################################################################
// Function Create
//#####################################################################
template<class TV> FREE_PARTICLES<TV>* FREE_PARTICLES<TV>::
Create()
{
    return new FREE_PARTICLES;
}
//#####################################################################
// Function Mark_Nodes_Referenced
//#####################################################################
template<class TV> void FREE_PARTICLES<TV>::
Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const
{
    marks.Subset(nodes).Fill(mark);
}
//#####################################################################
template class FREE_PARTICLES<VECTOR<float,1> >;
template class FREE_PARTICLES<VECTOR<float,2> >;
template class FREE_PARTICLES<VECTOR<float,3> >;
template class FREE_PARTICLES<VECTOR<double,1> >;
template class FREE_PARTICLES<VECTOR<double,2> >;
template class FREE_PARTICLES<VECTOR<double,3> >;
