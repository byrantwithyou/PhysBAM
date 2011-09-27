//#####################################################################
// Copyright 2011, Jeffrey Hellrung
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_WRITE_MESH_OBJECT_TO_FILE_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_WRITE_MESH_OBJECT_TO_FILE_HPP

#include <exception>
#include <string>

#include <boost/exception/exception.hpp>
#include <boost/spirit/home/qi/auxiliary/eoi.hpp>
#include <boost/spirit/home/qi/char/char.hpp>
#include <boost/spirit/home/qi/nonterminal/rule.hpp>
#include <boost/spirit/home/qi/operator/difference.hpp>
#include <boost/spirit/home/qi/operator/kleene.hpp>
#include <boost/spirit/home/qi/operator/optional.hpp>
#include <boost/spirit/home/qi/operator/sequence.hpp>
#include <boost/spirit/home/qi/parse.hpp>
#include <boost/spirit/home/qi/string/lit.hpp>
#include <boost/spirit/home/support/unused.hpp>
#include <boost/throw_exception.hpp>

#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/TYPED_STREAM.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_MESH_OBJECT.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <Jeffrey_Utilities/VTK/Write_Mesh_Object_To_VTK_File.h>

namespace PhysBAM
{

template< class T_MESH_OBJECT >
struct WRITE_MESH_OBJECT_TO_FILE_FILENAME_ERROR;

namespace Detail_Write_Mesh_Object_To_File
{

template< int D, class T_MESH >
inline const char*
PhysBAM_Extension();
inline bool
Filename_Has_PhysBAM_Extension(const std::string& filename, const char* const ext);
inline bool
Filename_Has_VTK_Extension(const std::string& filename);

} // namespace Detail_Write_Mesh_Object_To_File

template< class T, int D, class T_MESH >
inline void
Write_Mesh_Object_To_File(
    const MESH_OBJECT< VECTOR<T,D>, T_MESH >& mesh_object,
    const std::string& filename)
{
    const char* const physbam_ext =
        Detail_Write_Mesh_Object_To_File::PhysBAM_Extension< D, T_MESH >();
    if(Detail_Write_Mesh_Object_To_File::Filename_Has_PhysBAM_Extension(filename, physbam_ext)) {
        typedef float RW;
        static const PhysBAM::STREAM_TYPE rw((RW()));
        FILE_UTILITIES::Write_To_File(rw, filename, mesh_object);
    }
    else if(Detail_Write_Mesh_Object_To_File::Filename_Has_VTK_Extension(filename))
        Write_Mesh_Object_To_VTK_File(mesh_object, filename);
    else
        boost::throw_exception(
            WRITE_MESH_OBJECT_TO_FILE_FILENAME_ERROR<
                MESH_OBJECT< VECTOR<T,D>, T_MESH >
            >(filename)
        );
}

struct WRITE_MESH_OBJECT_TO_FILE_FILENAME_ERROR_BASE
    : virtual std::exception, virtual boost::exception
{
    std::string filename;
    explicit WRITE_MESH_OBJECT_TO_FILE_FILENAME_ERROR_BASE(const std::string& filename_)
        : filename(filename_)
    { }
    ~WRITE_MESH_OBJECT_TO_FILE_FILENAME_ERROR_BASE() throw ( )
    { }
    const char* what() const throw ( )
    { return "PhysBAM::WRITE_MESH_OBJECT_TO_FILE_FILENAME_ERROR_BASE"; }
};

template< class T_MESH_OBJECT >
struct WRITE_MESH_OBJECT_TO_FILE_FILENAME_ERROR
    : WRITE_MESH_OBJECT_TO_FILE_FILENAME_ERROR_BASE
{
    explicit WRITE_MESH_OBJECT_TO_FILE_FILENAME_ERROR(const std::string& filename)
        : WRITE_MESH_OBJECT_TO_FILE_FILENAME_ERROR_BASE(filename)
    { }
    const char* what() const throw ( )
    { return "PhysBAM::WRITE_MESH_OBJECT_TO_FILE_FILENAME_ERROR< T_MESH_OBJECT >"; }
};

namespace Detail_Write_Mesh_Object_To_File
{

#define DEFINE_PHYSBAM_EXTENSION( D, T_MESH, String ) \
template<> inline const char* PhysBAM_Extension< D, T_MESH >() { return String; }
DEFINE_PHYSBAM_EXTENSION( 2,     SEGMENT_MESH, "curve2d" )
DEFINE_PHYSBAM_EXTENSION( 2,    TRIANGLE_MESH, "tri2d"   )
DEFINE_PHYSBAM_EXTENSION( 3,     SEGMENT_MESH, "curve"   )
DEFINE_PHYSBAM_EXTENSION( 3,    TRIANGLE_MESH, "tri"     )
DEFINE_PHYSBAM_EXTENSION( 3, TETRAHEDRON_MESH, "tet"     )
#undef DEFINE_PHYSBAM_EXTENSION

namespace qi = boost::spirit::qi;

inline bool
Filename_Has_Extension(
    const std::string& filename,
    qi::rule<
        std::string::const_iterator,
        qi::unused_type ( )
    > dot_ext_eoi)
{
    std::string::const_iterator it = filename.begin();
    return qi::parse(
        it, filename.end(),
        *(qi::char_ - dot_ext_eoi) >> dot_ext_eoi
    );
}

inline bool
Filename_Has_PhysBAM_Extension(const std::string& filename, const char* const ext)
{ return Filename_Has_Extension(filename, '.' >> qi::lit(ext) >> -qi::lit(".gz") >> qi::eoi); }

inline bool
Filename_Has_VTK_Extension(const std::string& filename)
{ return Filename_Has_Extension(filename, ".vtk" >> qi::eoi); }

} // namespace Detail_Write_Mesh_Object_To_File

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_WRITE_MESH_OBJECT_TO_FILE_HPP
