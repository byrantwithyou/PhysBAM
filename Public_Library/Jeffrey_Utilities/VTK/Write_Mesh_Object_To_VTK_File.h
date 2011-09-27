//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VTK_WRITE_MESH_OBJECT_TO_VTK_FILE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VTK_WRITE_MESH_OBJECT_TO_VTK_FILE_HPP

#include <exception>
#include <string>

#include <boost/exception/exception.hpp>
#include <boost/throw_exception.hpp>

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <Jeffrey_Utilities/VTK/Write_Segmented_Curve_2d_To_VTK_File.h>
#include <Jeffrey_Utilities/VTK/Write_Triangulated_Surface_To_VTK_File.h>

namespace PhysBAM
{

template< class T_MESH_OBJECT >
struct WRITE_MESH_OBJECT_TO_VTK_FILE_TYPE_ERROR;

template< class T, int D, class T_MESH >
inline void
Write_Mesh_Object_To_VTK_File(
    const MESH_OBJECT< VECTOR<T,D>, T_MESH >& mesh_object,
    const std::string& filename)
{
    boost::throw_exception(
        WRITE_MESH_OBJECT_TO_VTK_FILE_TYPE_ERROR<
            MESH_OBJECT< VECTOR<T,D>, T_MESH >
        >()
    );
}

template< class T >
inline void
Write_Mesh_Object_To_VTK_File(
    const MESH_OBJECT< VECTOR<T,2>, SEGMENT_MESH >& segmented_curve_2d,
    const std::string& filename)
{ Write_Segmented_Curve_2d_To_VTK_File(segmented_curve_2d, filename); }

template< class T >
inline void
Write_Mesh_Object_To_VTK_File(
    const MESH_OBJECT< VECTOR<T,3>, TRIANGLE_MESH >& triangulated_surface,
    const std::string& filename)
{ Write_Triangulated_Surface_To_VTK_File(triangulated_surface, filename); }

struct WRITE_MESH_OBJECT_TO_VTK_FILE_TYPE_ERROR_BASE
    : virtual std::exception, virtual boost::exception
{
    const char* what() const throw ( )
    { return "PhysBAM::WRITE_MESH_OBJECT_TO_VTK_FILE_TYPE_ERROR_BASE"; }
};

template< class T_MESH_OBJECT >
struct WRITE_MESH_OBJECT_TO_VTK_FILE_TYPE_ERROR
    : WRITE_MESH_OBJECT_TO_VTK_FILE_TYPE_ERROR_BASE
{
    const char* what() const throw ( )
    { return "PhysBAM::WRITE_MESH_OBJECT_TO_VTK_FILE_TYPE_ERROR< T_MESH_OBJECT >"; }
};

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_VTK_WRITE_MESH_OBJECT_TO_VTK_FILE_HPP
