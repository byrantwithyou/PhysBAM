//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_TRIANGULATED_OBJECT
//#####################################################################
#ifndef __EMBEDDED_TRIANGULATED_OBJECT__
#define __EMBEDDED_TRIANGULATED_OBJECT__

#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Fracture/EMBEDDED_OBJECT.h>
namespace PhysBAM{

template<class TV>
class EMBEDDED_TRIANGULATED_OBJECT:public EMBEDDED_OBJECT<TV,2>
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::SEGMENTED_CURVE T_SEGMENTED_CURVE;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT T_TRIANGULATED_OBJECT;
public:
    typedef EMBEDDED_OBJECT<TV,2> BASE;
    using BASE::embedded_particles;using BASE::parent_particles;using BASE::Is_Parent;using BASE::Embedded_Particle_On_Segment;using BASE::embedded_mesh;using BASE::simplicial_object;
    using BASE::Embedded_Subelements_In_Element;using BASE::Number_Of_Embedded_Subelements_In_Element;using BASE::Node_In_Simplex_Is_Material;using BASE::Add_Embedded_Subelement;

    EMBEDDED_TRIANGULATED_OBJECT(T_TRIANGULATED_OBJECT& simplicial_object_input);

    bool Nodes_Are_Materially_Connected_In_Simplex(const int node1,const int node2,const int simplex) const
    {return Node_In_Simplex_Is_Material(node1,simplex) && Node_In_Simplex_Is_Material(node2,simplex) && !Nodes_Are_Separated_In_Simplex(node1,node2,simplex);}

    int Number_Of_Embedded_Cuts(const int triangle) // a quad only counts as one cut
    {return Number_Of_Embedded_Subelements_In_Element(triangle);}
   
    int Add_Embedded_Segment(const int embedded_particle1,const int embedded_particle2)
    {return Add_Embedded_Subelement(VECTOR<int,2>(embedded_particle1,embedded_particle2));}

    int Add_Embedded_Segment(const int embedded_particle1_parent1,const int embedded_particle1_parent2,const int embedded_particle2_parent1,
        const int embedded_particle2_parent2)
    {return Add_Embedded_Segment(Embedded_Particle_On_Segment(embedded_particle1_parent1,embedded_particle1_parent2),
        Embedded_Particle_On_Segment(embedded_particle2_parent1,embedded_particle2_parent2));}

    int Node_Separated_By_Embedded_Subelement(const int embedded_segment) const;
    bool Nodes_Are_Separated_In_Simplex(const int node1,const int node2,const int triangle) const;
    bool Segment_Is_Broken(const int node1,const int node2) const;
    int Number_Of_Edges_With_Embedded_Particles(const int triangle);
    int Embedded_Node_Common_To_Both_Segments_In_Triangle(const int triangle);
    int Isolated_Node(const int triangle);
    int Diamond_Node(const int triangle);
    T Fraction_Of_Triangles_With_N_Cuts(const int n);
//#####################################################################
};
}
#endif
