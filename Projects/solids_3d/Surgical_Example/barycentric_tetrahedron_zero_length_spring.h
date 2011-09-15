//#####################################################################
// Copyright 2006, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING  
//##################################################################### 
#ifndef __BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING__
#define __BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING__

#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/SOLIDS_FORCES.h>
namespace PhysBAM{

template<class T>
class BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING:public SOLIDS_FORCES<T,VECTOR<T,3> >
{
public:
    typedef SOLIDS_FORCES<T,VECTOR<T,3> > BASE;
    using BASE::particles;
    using BASE::CFL_elastic_time_step_of_fragment;
    using BASE::CFL_damping_time_step_of_fragment;
    using BASE::CFL_initialized_of_fragment;

public:  
    T constant_youngs_modulus; // units of force (i.e. force per unit strain)
    T constant_damping; // units of force*time (i.e. force per unit strain rate)
    ARRAY<T> youngs_modulus,damping;
    TETRAHEDRON_MESH& tetrahedron_mesh;
    ARRAY<VECTOR<int,2> >& constrained_tets;
    ARRAY<VECTOR<VECTOR<T,3>,2> >& constrained_location_weights;
    ARRAY<VECTOR<MATRIX_4X4<T>,4> > weights_outer_product;

public: 
    BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING(PARTICLES<T,VECTOR<T,3> >& particles_input,ARRAY<VECTOR<int,2> >& constrained_tets_input,ARRAY<VECTOR<VECTOR<T,3>,2> >& constrained_location_weights_input,TETRAHEDRON_MESH& mesh_input,T stiffness,T damping)
        :SOLIDS_FORCES<T,VECTOR<T,3> >(particles_input),constrained_tets(constrained_tets_input),tetrahedron_mesh(mesh_input),constrained_location_weights(constrained_location_weights_input)
    {
        Set_Stiffness(stiffness);Set_Damping(damping);
    }
    
    ~BARYCENTRIC_TETRAHEDRON_ZERO_LENGTH_SPRING()
    {}

    void Set_Stiffness(const T youngs_modulus_input)
    {constant_youngs_modulus=youngs_modulus_input;}
    
    void Set_Damping(const T rayleigh_coefficient)
    {constant_damping=rayleigh_coefficient*constant_youngs_modulus;}
    
//#####################################################################
    void Initialize_Weights_Outer_Product();
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<VECTOR<T,3> > F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const VECTOR<T,3> > V,ARRAY_VIEW<VECTOR<T,3> > F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Force_Differential(ARRAY_VIEW<const VECTOR<T,3> > dX,ARRAY_VIEW<VECTOR<T,3> > dF,const T time) const PHYSBAM_OVERRIDE;
    void Initialize_CFL() PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE;
//#####################################################################
};   
}
#endif



