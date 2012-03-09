//#####################################################################
// Copyright 2003-2005, Neil Molino, Andrew Selle, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <cstring>
#include <fstream>
using namespace PhysBAM;
using namespace std;

//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    TETRAHEDRON_MESH tet_mesh;
    DEFORMABLE_PARTICLES<VECTOR<float,3> > particles;
    TETRAHEDRALIZED_VOLUME<float> tet_vol(tet_mesh,particles);
    FILE_UTILITIES::Read_From_File<float>(argv[1],tet_vol);
    std::cout<<"tet_vol.particles.number="<<tet_vol.particles.array_collection->Size()<<std::endl;
    std::cout<<"tet_vol.mesh.elements.m="<<tet_vol.mesh.elements.m<<std::endl;
    tet_vol.Update_Bounding_Box();
    std::cout<<"bounding_box size:"<<tet_vol.bounding_box->Size()<<std::endl;
    tet_vol.Initialize_Triangulated_Surface();
    tet_vol.triangulated_surface->Discard_Valence_Zero_Particles_And_Renumber();
    cout<<"triangles.m="<<tet_vol.triangulated_surface->mesh.elements.m<<endl;
    cout<<"particles.number="<<particles.array_collection->Size()<<endl;
    FILE_UTILITIES::Write_To_File<float>(argv[2],*tet_vol.triangulated_surface);
    return 0;
}
//#################################################################
