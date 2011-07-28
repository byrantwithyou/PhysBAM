//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_CLOSE_MESH_OBJECT_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_CLOSE_MESH_OBJECT_HPP

#include <boost/foreach.hpp>

#include <Jeffrey_Utilities/Mesh/Remove_Degenerate_Mesh_Elements.h>
#include <Jeffrey_Utilities/VECTOR_OPS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

template< class T, int D, class T_MESH >
void
Close_Mesh_Object(MESH_OBJECT< VECTOR<T,D>, T_MESH >& mesh_object)
{
    typedef MESH_OBJECT< VECTOR<T,D>, T_MESH > MESH_OBJECT_TYPE;

    {
        struct SCOPED_MESH_OBJECT_MESH_INCIDENT_ELEMENTS
        {
            MESH_OBJECT_TYPE& mesh_object;
            const bool uninitialized;
            explicit SCOPED_MESH_OBJECT_MESH_INCIDENT_ELEMENTS(MESH_OBJECT_TYPE& mesh_object)
                : mesh_object(mesh_object),
                  uninitialized(mesh_object.mesh.incident_elements == 0)
            {
                if(uninitialized)
                    mesh_object.mesh.Initialize_Incident_Elements();
            }
            ~SCOPED_MESH_OBJECT_MESH_INCIDENT_ELEMENTS()
            {
                if(uninitialized) {
                    delete mesh_object.mesh.incident_elements;
                    mesh_object.mesh.incident_elements = 0;
                }
            }
        } _scoped_incident_elements(mesh_object);

        HASHTABLE< VECTOR<T,D>, int > index_of_x(mesh_object.particles.array_collection->Size());
        for(int p = 1; p <= mesh_object.particles.array_collection->Size(); ++p) {
            const int q = index_of_x.Get_Or_Insert(mesh_object.particles.X(p), p);
            if(q != p) {
                BOOST_FOREACH( const int e, (*mesh_object.mesh.incident_elements)(p) )
                    Replace(mesh_object.mesh.elements(e), p, q);
            }
        }
    }

    Remove_Degenerate_Mesh_Elements(mesh_object.mesh.elements);
    mesh_object.Refresh_Auxiliary_Structures();
    mesh_object.Discard_Valence_Zero_Particles_And_Renumber();
}

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_GEOMETRY_CLOSE_MESH_OBJECT_HPP
