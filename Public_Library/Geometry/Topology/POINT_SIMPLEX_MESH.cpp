//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINT_SIMPLEX_MESH
//#####################################################################
#include <Geometry/Topology/POINT_SIMPLEX_MESH.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
POINT_SIMPLEX_MESH::
POINT_SIMPLEX_MESH():segment_mesh(0)
{}
//#####################################################################
// Destructor
//#####################################################################
POINT_SIMPLEX_MESH::
~POINT_SIMPLEX_MESH()
{
    delete segment_mesh;
}
//#####################################################################
// Function Initialize_Segment_Mesh
//#####################################################################
void POINT_SIMPLEX_MESH::
Initialize_Segment_Mesh()
{
    delete segment_mesh;

    ARRAY<VECTOR<int,2> > simplex_list;
    simplex_list.Resize(number_nodes-1);
    for(int i=0;i<number_nodes-1;i++) simplex_list(i)=VECTOR<int,2>(elements(i)[0],elements(i+1)[0]);
    segment_mesh=new SEGMENT_MESH(number_nodes-1,simplex_list);
}
//#####################################################################
