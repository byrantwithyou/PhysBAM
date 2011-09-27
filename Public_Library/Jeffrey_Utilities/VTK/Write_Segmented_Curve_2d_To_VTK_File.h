//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VTK_WRITE_SEGMENTED_CURVE_2D_TO_VTK_FILE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VTK_WRITE_SEGMENTED_CURVE_2D_TO_VTK_FILE_HPP

#include <string>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>

namespace PhysBAM
{

template< class T >
inline void
Write_Segmented_Curve_2d_To_VTK_File(
    const MESH_OBJECT< VECTOR<T,2>, SEGMENT_MESH >& mesh_object,
    const std::string& filename)
{
    // TODO
    static_cast<void>(mesh_object);
    static_cast<void>(filename);
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VTK_WRITE_SEGMENTED_CURVE_2D_TO_VTK_FILE_HPP
