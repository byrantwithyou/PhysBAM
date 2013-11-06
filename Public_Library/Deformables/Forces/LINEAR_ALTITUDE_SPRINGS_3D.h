//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Igor Neverov, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_ALTITUDE_SPRINGS_3D
//#####################################################################
#ifndef __LINEAR_ALTITUDE_SPRINGS_3D__
#define __LINEAR_ALTITUDE_SPRINGS_3D__

#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Deformables/Forces/LINEAR_ALTITUDE_SPRINGS.h>
namespace PhysBAM{

template<class T_input>
class LINEAR_ALTITUDE_SPRINGS_3D:public LINEAR_ALTITUDE_SPRINGS<VECTOR<T_input,3>,3>
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;typedef typename TV::SPIN T_SPIN;
    typedef LINEAR_ALTITUDE_SPRINGS<VECTOR<T,3>,3> BASE;typedef typename BASE::SPRING_STATE SPRING_STATE;typedef typename BASE::SPRING_PARAMETER SPRING_PARAMETER;
    using BASE::use_springs_compressed_beyond_threshold;using BASE::spring_compression_fraction_threshold;using BASE::print_number_used;
    using BASE::force_elements;
public:
    using BASE::particles;using BASE::Invalidate_CFL;
    using BASE::max_strain_per_time_step;using BASE::mesh;using BASE::use_rest_state_for_strain_rate;
    using BASE::spring_states;using BASE::spring_states_all_springs;
    using BASE::use_plasticity;using BASE::use_shortest_spring_only;
    using BASE::cache_strain;using BASE::strains_of_spring;using BASE::strains_of_spring_all_springs;
    using BASE::compute_half_forces;
    typedef typename BASE::ELEMENT_ITERATOR ELEMENT_ITERATOR;
    using BASE::parameters;using BASE::Compute_Plasticity;

    LINEAR_ALTITUDE_SPRINGS_3D(DEFORMABLE_PARTICLES<TV>& particles,TETRAHEDRON_MESH& tetrahedron_mesh);
    virtual ~LINEAR_ALTITUDE_SPRINGS_3D();

//#####################################################################
    void Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient); // assumes mass and restlength are already defined
    void Set_Restlength_From_Particles();
    void Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<const TV> material_coordinates);
    void Set_Overdamping_Fraction(const T overdamping_fraction); // 1 is critically damped
    void Set_Overdamping_Fraction(const ARRAY<VECTOR<T,4> >& overdamping_fraction); // 1 is critically damped
    void Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction=1); // 1 is critically damped
    void Fill_Node_Indices(int i,int j,int k,int l,int isolated_node_number,int& node1,int& node2,int& node3,int& node4) const;
    bool Fill_Spring_State(int t,int isolated_node_number,int node1,int node2,int node3,int node4,SPRING_STATE& state);
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T scale,const T time) const PHYSBAM_OVERRIDE;
    bool Compute_Strain_Rate_And_Strain(int t,int isolated_node_number,int node1,int node2,int node3,int node4,T& strain_rate,T& strain) const;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    void Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name="") const PHYSBAM_OVERRIDE;

    T Potential_Energy(const int t,const T time) const;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};

template<class T> LINEAR_ALTITUDE_SPRINGS_3D<T>*
Create_Altitude_Springs(DEFORMABLE_PARTICLES<VECTOR<T,3> >& particles,TETRAHEDRON_MESH& mesh,
    const T stiffness=200,const T overdamping_fraction=2,const bool use_compressed_by_threshold_only=true,const T fraction_compression=.1,const bool limit_time_step_by_strain_rate=true,
    const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,const T restlength_enlargement_fraction=0,const bool verbose=true);

template<class T> LINEAR_ALTITUDE_SPRINGS_3D<T>*
Create_Altitude_Springs(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,
    const T stiffness=4/(1+sqrt((T)2)),const T overdamping_fraction=2,const bool use_compressed_by_threshold_only=true,const T fraction_compression=.1,
    const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,const T restlength_enlargement_fraction=0,
    const bool verbose=true);

}
#endif
