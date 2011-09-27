//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VTK_WRITE_TRIANGULATED_SURFACE_TO_VTK_FILE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VTK_WRITE_TRIANGULATED_SURFACE_TO_VTK_FILE_HPP

#include <fstream>
#include <string>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <Jeffrey_Utilities/typeid_name.h>

namespace PhysBAM
{

template< class T >
inline void
Write_Triangulated_Surface_To_VTK_File(
    const MESH_OBJECT< VECTOR<T,3>, TRIANGLE_MESH >& mesh_object,
    const std::string& filename)
{
    std::fstream fout(filename.c_str(), std::ios_base::out);

    fout << "# vtk DataFile Version 1.0\n" // vtk file version
         << filename << '\n'               // header
         << "ASCII\n\n"                    // file format
         << "DATASET POLYDATA"             // dataset structure
         << std::endl;

    const int n_particle = mesh_object.particles.array_collection->Size();
    fout << "POINTS " << n_particle << ' ' << typeid_name<T>() << std::endl;
    for(int p = 1; p <= n_particle; ++p) {
        const VECTOR<T,3> x = mesh_object.particles.X(p);
        fout << x[1] << ' ' << x[2] << ' ' << x[3] << '\n';
    }
    fout << std::endl;

    const int n_element = mesh_object.mesh.elements.Size();
    fout << "POLYGONS " << n_element << ' ' << 4 * n_element << std::endl;
    for(int e = 1; e <= n_element; ++e) {
        const VECTOR<int,3> ee = mesh_object.mesh.elements(e);
        fout << "3 " << ee[1] - 1 << ' ' << ee[2] - 1 << ' ' << ee[3] - 1 << '\n';
    }
    fout << std::endl;
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VTK_WRITE_TRIANGULATED_SURFACE_TO_VTK_FILE_HPP
