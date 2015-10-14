//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Neil Molino, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_TETRAHEDRALIZED_VOLUME
//#####################################################################
#ifndef __EMBEDDED_TETRAHEDRALIZED_VOLUME__
#define __EMBEDDED_TETRAHEDRALIZED_VOLUME__

#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Deformables/Fracture/EMBEDDED_OBJECT.h>
namespace PhysBAM{

template<class T_input>
class EMBEDDED_TETRAHEDRALIZED_VOLUME:public EMBEDDED_OBJECT<VECTOR<T_input,3>,3>
{
    typedef T_input T;typedef VECTOR<T,3> TV;
public:
    typedef EMBEDDED_OBJECT<TV,3> BASE;
    using BASE::embedded_particles;using BASE::parent_particles;using BASE::Is_Parent;using BASE::Embedded_Particle_On_Segment;using BASE::simplicial_object;using BASE::embedded_mesh;
    using BASE::Embedded_Subelements_In_Element;using BASE::Number_Of_Embedded_Subelements_In_Element;using BASE::Node_In_Simplex_Is_Material;using BASE::Add_Embedded_Subelement;

    EMBEDDED_TETRAHEDRALIZED_VOLUME(TETRAHEDRALIZED_VOLUME<T>& simplicial_object_input);
    int Node_Separated_By_Embedded_Subelement(const int embedded_triangle) const;
    bool Nodes_Are_Separated_In_Simplex(const int node1,const int node2,const int tetrahedron) const;
    bool Nodes_Are_Materially_Connected_In_Simplex(const int node1,const int node2,const int simplex) const;
    int Number_Of_Embedded_Cuts(const int tetrahedron) const; // a quad only counts as one cut
    bool Cut_By_Quad(const int tetrahedron) const;
    bool Segment_Is_Broken(const int node1,const int node2) const;
    int Number_Of_Edges_With_Embedded_Particles(const int tetrahedron) const;
    int Number_Of_Tetrahedra_With_Cuts() const;
    int Number_Of_Tetrahedra_With_N_Cuts(const int n) const;

    T Fraction_Of_Tetrahedra_With_N_Cuts(const int n) const
    {return Number_Of_Tetrahedra_With_N_Cuts(n)/(T)simplicial_object.mesh.elements.m;}

    int Add_Embedded_Triangle(const int embedded_particle1,const int embedded_particle2,const int embedded_particle3)
    {return Add_Embedded_Subelement(VECTOR<int,3>(embedded_particle1,embedded_particle2,embedded_particle3));}
        
    int Add_Embedded_Triangle(const int embedded_particle1_parent1,const int embedded_particle1_parent2,const int embedded_particle2_parent1,
        const int embedded_particle2_parent2,const int embedded_particle3_parent1,const int embedded_particle3_parent2)
    {return Add_Embedded_Triangle(Embedded_Particle_On_Segment(embedded_particle1_parent1,embedded_particle1_parent2),
        Embedded_Particle_On_Segment(embedded_particle2_parent1,embedded_particle2_parent2),Embedded_Particle_On_Segment(embedded_particle3_parent1,embedded_particle3_parent2));}

//#####################################################################
};
}
#endif
