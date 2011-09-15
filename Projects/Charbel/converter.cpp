//#####################################################################
// Copyright 2009, Jon Gretarsson, Kevin Wang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

using namespace PhysBAM;
int main()
{
  STREAM_TYPE stream_type((double)0);
  int n1, n2, n3, n4, n5;
  double x1, x2, x3;
  FILE* structure = fopen("corrected_surface.top","r");
  int numNodes, numElems;
  if(fscanf(structure, "%d %d", &numNodes, &numElems)!=2) return 1;

  GEOMETRY_PARTICLES<VECTOR<double,3> >& physbam_geometry_particle
    = *new GEOMETRY_PARTICLES<VECTOR<double,3> >();
  physbam_geometry_particle.array_collection->Add_Elements(numNodes);

  for (int i=0; i<numNodes; i++) {
    if(fscanf(structure, "%d %lf %lf %lf\n", &n1, &x1, &x2, &x3)!=4) return 1;
    physbam_geometry_particle.X(i+1) = VECTOR<double,3>(x1,x2,x3);
  }

  ARRAY<VECTOR<int,3> > & physbam_triangle_list
    = *new ARRAY<VECTOR<int,3> >();
  for (int i=0; i<numElems; i++) {
    if(fscanf(structure, "%d %d %d %d %d\n", &n1, &n2, &n3, &n4, &n5)!=5) return 1;
    physbam_triangle_list.Append(VECTOR<int,3>(n3, n4, n5));
  }

  TRIANGLE_MESH& physbam_triangle_mesh
    = *new TRIANGLE_MESH(physbam_geometry_particle.array_collection->Size(), physbam_triangle_list);
  physbam_triangle_mesh.Initialize_Adjacent_Elements();

  // Construct TRIANGULATED_SURFACE.
  TRIANGULATED_SURFACE<double>& physbam_triangulated_surface
    = *new TRIANGULATED_SURFACE<double>(physbam_triangle_mesh, physbam_geometry_particle);
  physbam_triangulated_surface.Update_Triangle_List();
  FILE_UTILITIES::Write_To_File(stream_type,"corrected_surface.tri",physbam_triangulated_surface);
}
