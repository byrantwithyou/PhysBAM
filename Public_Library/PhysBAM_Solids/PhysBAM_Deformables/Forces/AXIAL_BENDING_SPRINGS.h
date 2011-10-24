//#####################################################################
// Copyright 2004-2008, Eran Guendelman, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AXIAL_BENDING_SPRINGS
//#####################################################################
// A spring force between the edge connecting two triangles and the cross-edge connecting their two opposite vertices.
// The spring corresponds to the shortest line connecting the two edge lines (connections can be anywhere on the infinite lines
// supporting the edges, not restricted to being on the edges themselves.
//#####################################################################
#ifndef __AXIAL_BENDING_SPRINGS__
#define __AXIAL_BENDING_SPRINGS__

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

class TRIANGLE_MESH;

template<class T_input>
class AXIAL_BENDING_SPRINGS:public DEFORMABLES_FORCES<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;using BASE::max_strain_per_time_step;using BASE::use_rest_state_for_strain_rate;using BASE::Invalidate_CFL;using BASE::cfl_number;
    typedef typename FORCE_ELEMENTS::ITERATOR SPRING_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    TRIANGLE_MESH& triangle_mesh;
    ARRAY<T> youngs_modulus; // units of force (i.e. force per unit strain)
    ARRAY<T> restlength,visual_restlength; // visual restlength corresponds to length between particles; restlength may be larger than this to avoid zero/small restlength
    ARRAY<T> damping; // units of force*time (i.e. force per unit strain rate)
    ARRAY<VECTOR<int,4> > spring_particles; // spring is shortest line between segment with particles (1,2) and segment with particles (3,4)
    ARRAY<T> attached_edge_restlength; // total rest length of edges of attached triangles which are not the center edge
    ARRAY<T> attached_edge_length; // total length of edges of attached triangles which are not the center edge
private:
    ARRAY<T> optimization_coefficient;
    ARRAY<T> optimization_current_length;
    ARRAY<TV> optimization_direction;
    ARRAY<VECTOR<T,4> > optimization_weights;
    FORCE_ELEMENTS force_springs;
public:
    bool verbose;

    AXIAL_BENDING_SPRINGS(PARTICLES<TV>& particles_input,TRIANGLE_MESH& triangle_mesh_input);

    virtual ~AXIAL_BENDING_SPRINGS();

    void Set_Stiffness(const T youngs_modulus_input)
    {ARRAYS_COMPUTATIONS::Fill(youngs_modulus,youngs_modulus_input);Invalidate_CFL();}

    void Set_Stiffness(ARRAY_VIEW<const T> youngs_modulus_input)
    {ARRAY<T>::Copy(youngs_modulus_input,youngs_modulus);Invalidate_CFL();}

    void Clamp_Restlength(const T clamped_restlength)
    {for(int i=1;i<=restlength.m;i++) restlength(i)=max(visual_restlength(i),clamped_restlength);}

    void Set_Damping(const T damping_input)
    {ARRAYS_COMPUTATIONS::Fill(damping,damping_input);Invalidate_CFL();}

    void Set_Damping(ARRAY_VIEW<const T> damping_input)
    {ARRAY<T>::Copy(damping_input,damping);Invalidate_CFL();}

//#####################################################################
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Initialize();
    void Set_Restlength_From_Particles();
    void Set_Overdamping_Fraction(const T overdamping_fraction=1); // 1 is critically damped
    void Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction); // 1 is critically damped
    void Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction=1); // 1 is critically damped
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    void Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient); // assumes mass and restlength are already defined
    void Axial_Vector(const VECTOR<int,4>& nodes,T& axial_length,TV& axial_direction,VECTOR<T,2>& weights,T& attached_edge_length) const;

    T Potential_Energy(int s,const T time) const;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    T Endpoint_Mass(int s,int b) const;
    TV Endpoint_Velocity(int s,int b) const;
    TV Endpoint_Velocity(ARRAY_VIEW<const TV> velocity,int s,int b) const;
    T Endpoint_Kinetic_Energy(int s,int b) const;
    T Endpoint_Kinetic_Energy(ARRAY_VIEW<const TV> velocity,int s,int b) const;
    T Endpoint_Kinetic_Energy(int s) const;
    T Effective_Impulse_Factor(int s) const;
//#####################################################################
};

template<class T> AXIAL_BENDING_SPRINGS<T>*
Create_Axial_Bending_Springs(PARTICLES<VECTOR<T,3> >& particles,TRIANGLE_MESH& triangle_mesh,const T clamped_restlength=(T).01,
    const T stiffness=2/(1+sqrt((T)2)),const T overdamping_fraction=2,const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,
    const bool use_rest_state_for_strain_rate=true,const bool verbose=true);

template<class T> AXIAL_BENDING_SPRINGS<T>*
Create_Axial_Bending_Springs(TRIANGULATED_SURFACE<T>& triangulated_surface,
    const T clamped_restlength=(T).01,const T stiffness=2/(1+sqrt((T)2)),const T overdamping_fraction=2,const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,
    const bool use_rest_state_for_strain_rate=true,const bool verbose=true);
}
#endif
