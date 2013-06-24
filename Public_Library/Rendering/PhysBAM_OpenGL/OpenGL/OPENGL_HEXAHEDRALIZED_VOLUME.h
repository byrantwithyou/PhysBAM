//#####################################################################
// Copyright 2002-2005, Eilene Hao, Geoffrey Irving, Sergey Koltakov, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_HEXAHEDRALIZED_VOLUME
//##################################################################### 
#ifndef __OPENGL_HEXAHEDRALIZED_VOLUME__
#define __OPENGL_HEXAHEDRALIZED_VOLUME__

#include <Tools/Vectors/VECTOR_3D.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Topology/HEXAHEDRON_MESH.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
namespace PhysBAM{

template<class T>
class OPENGL_HEXAHEDRALIZED_VOLUME:public OPENGL_OBJECT
{
public:
    OPENGL_MATERIAL material;
    OPENGL_MATERIAL inverted_material;
    bool use_inverted_material;
    HEXAHEDRON_MESH* hexahedron_mesh;
    const GEOMETRY_PARTICLES<VECTOR<T,3> >* particles;
    int current_hexahedron;
    bool boundary_only,draw_subsets;
    ARRAY<int> subset;
    ARRAY<int> subset_particles;
    ARRAY<VECTOR<T,3> > vectors_at_hex_centers;
    T vector_size;

    OPENGL_HEXAHEDRALIZED_VOLUME(const OPENGL_MATERIAL& material_input,const OPENGL_MATERIAL& inverted_material_input);
    OPENGL_HEXAHEDRALIZED_VOLUME(HEXAHEDRON_MESH* hexahedron_mesh_input,const GEOMETRY_PARTICLES<VECTOR<T,3> >* particles_input,const OPENGL_MATERIAL& material_input,
        const OPENGL_MATERIAL& inverted_material_input);
    virtual ~OPENGL_HEXAHEDRALIZED_VOLUME();

    void Initialize()
    {if(!hexahedron_mesh->boundary_mesh) hexahedron_mesh->Initialize_Boundary_Mesh();}

    bool Toggle_Boundary_Only()
    {boundary_only=!boundary_only;return boundary_only;}

    bool Toggle_Differentiate_Inverted()
    {use_inverted_material=!use_inverted_material;return use_inverted_material;}

//#####################################################################
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    void Draw_Subset_Particles() const;
    void Draw_Vector_At_Hex_Center() const;
    void Draw_Wireframe_Mesh(const HEXAHEDRON_MESH& hexahedron_mesh) const;
    void Draw_Boundary_Triangles(const HEXAHEDRON_MESH& hexahedron_mesh) const;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
