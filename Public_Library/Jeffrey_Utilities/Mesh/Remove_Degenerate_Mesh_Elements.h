//#####################################################################
// Copyright 2011, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the
// license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################

#ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MESH_REMOVE_DEGENERATE_MESH_ELEMENTS_HPP
#define PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MESH_REMOVE_DEGENERATE_MESH_ELEMENTS_HPP

#include <Jeffrey_Utilities/Remove_If.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM
{

namespace Detail_Remove_Degenerate_Mesh_Elements
{

struct IS_DEGENERATE_MESH_ELEMENT;

} // namespace Detail_Remove_Degenerate_Mesh_Elements

template< int D >
inline void
Remove_Degenerate_Mesh_Elements(ARRAY< VECTOR<int,D> >& mesh_elements)
{
    typedef Detail_Remove_Degenerate_Mesh_Elements::IS_DEGENERATE_MESH_ELEMENT IS_DEGENERATE_MESH_ELEMENT_;
    Remove_If(mesh_elements, IS_DEGENERATE_MESH_ELEMENT_());
}

namespace Detail_Remove_Degenerate_Mesh_Elements
{

struct IS_DEGENERATE_MESH_ELEMENT
{
    typedef bool result_type;
    bool operator()(const VECTOR<int,1>& /*e*/) const
    { return false; }
    bool operator()(const VECTOR<int,2> e) const
    { return e[1] == e[2]; }
    bool operator()(const VECTOR<int,3> e) const
    { return e[1] == e[2] || e[1] == e[3] || e[2] == e[3]; }
    template< int D >
    bool operator()(const VECTOR<int,D> e) const
    {
        for(int i = 1; i < D; ++i)
            for(int j = i + 1; j <= D; ++j)
                if(e[i] == e[j])
                    return true;
    }
};

} // namespace Detail_Remove_Degenerate_Mesh_Elements

} // namespace PhysBAM

#endif // #ifndef PHYSBAM_PUBLIC_LIBRARY_JEFFREY_UTILITIES_MESH_REMOVE_DEGENERATE_MESH_ELEMENTS_HPP
