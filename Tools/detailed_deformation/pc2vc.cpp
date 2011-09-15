//#####################################################################
// Copyright 2006, Jerry Talton
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/SIMPLEX_INTERACTIONS.h>
#include <cstring>
#include <fstream>

using namespace PhysBAM;
using namespace std;


//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{

    if (argc!=2){cout<<"Usage: pc2vc input_colors"<<endl;exit(-1);}

    ARRAY<int> particle_colors;
    FILE_UTILITIES::Read_From_File<float>(argv[1],particle_colors);

    ARRAY<VECTOR_3D<float> > vertex_colors(particle_colors.m);

    for(int i=1;i<=vertex_colors.m;i++)vertex_colors(i)=(particle_colors(i)==1?VECTOR_3D<float>(.98,1.,.87):VECTOR_3D<float>(1,.37,.36));

    FILE_UTILITIES::Write_To_File<float>("particle_vertex_colors.col.gz",vertex_colors);
    
    return 0;
    
} 

