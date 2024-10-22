//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Geoffrey Irving, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE
//##################################################################### 
#ifndef __EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE__
#define __EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE__

#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/LOG.h>
#include <Deformables/Fracture/EMBEDDED_MATERIAL_SURFACE.h>
namespace PhysBAM{

template<class T_input>
class EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE:public EMBEDDED_MATERIAL_SURFACE<VECTOR<T_input,3>,3>
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef EMBEDDED_MATERIAL_SURFACE<TV,3> BASE;
    using BASE::particles;using BASE::material_surface_mesh;using BASE::material_surface;using BASE::embedded_object;
    using BASE::previously_perturbed;using BASE::parent_particles;using BASE::embedded_particles;

    EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE(EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_object);
    ~EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE();

//#####################################################################
    void Create_Material_Surface_From_Manifold_Embedded_Surface(bool verbose=true);
    void Conservative_Perturb_Nodes_For_Collision_Freeness(const T perturb_amount,const ARRAY<bool>& particle_on_surface);
    void Perturb_Nodes_For_Collision_Freeness(const T perturb_amount) override;
private:
    bool Center_Octahedron_In_Material(const int tetrahedron);
    void Construct_Material_Surface_Mesh() override;
    void Merge_Or_Cancel_Duplicate_Triangles();
    void Add_To_Material_Surface_Tetrahedron(const int i,const int j,const int k,const int l);
    void Add_To_Material_Surface_Tetrahedron_Face(const int i,const int j,const int l,const bool is_clockwise);
    void Add_To_Material_Surface_Triangle(int x0,int x1,int x2,const bool is_clockwise);
    void Add_To_Material_Surface_Planar_Quad(const int x0,const int x1,const int x2,const int x3,const bool is_clockwise);
    void Add_To_Material_Surface_Quad_Cut(const int il,const int jl,const int jk,const int ik,const bool is_clockwise);
    void Add_To_Material_Surface_Subtetrahedron_And_Subprism(const int tetrahedron,const int embedded_triangle1,const int i,const int j,const int k,const int l);
    void Add_To_Material_Surface_Subtetrahedron(const int i,const int j,const int k,const int l,const bool is_clockwise);
    void Add_To_Material_Surface_Subprism(const int i,const int j,const int k,const int l,const bool is_clockwise);
    void Add_To_Material_Surface_Subtetrahedrons_And_Wrick(const int tetrahedron,const int embedded_triangle1,const int embedded_triangle2,const int i,const int j,const int k,const int l);
    void Add_To_Material_Surface_Wrick(const int i,const int j,const int k,const int l,const bool is_clockwise);  
    void Add_To_Material_Surface_Wedge_On_Either_Side(const int tetrahedron,const int embedded_triangle1,const int embedded_triangle2,const int i,const int j,const int k,const int l);
    void Add_To_Material_Surface_Subwedge(const int i,const int j,const int k,const int l,const bool is_clockwise);    
    void Add_To_Material_Surface_Subtetrahedron_Tets_And_Oct_Plus_Tet(const int tetrahedron,const int i,const int j,const int k,const int l);
    void Add_To_Material_Surface_Oct_Plus_Tet(const int i,const int j,const int k,const int l,const bool is_clockwise);        // j is material
    void Add_To_Material_Surface_Half_Oct_Plus_Tet(const int i,const int j,const int k,const int l,const bool is_clockwise);   // j is material
    void Add_To_Material_Surface_Subtetrahedron_And_Wedge_And_Half_Oct_Plus_Tet(const int tetrahedron,const int i,const int j,const int k,const int l);
//#####################################################################
};
}
#endif
