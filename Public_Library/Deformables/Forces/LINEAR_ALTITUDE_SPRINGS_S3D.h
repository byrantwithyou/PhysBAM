//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ronald Fedkiw, Neil Molino, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_ALTITUDE_SPRINGS_S3D
//#####################################################################
#ifndef __LINEAR_ALTITUDE_SPRINGS_S3D__
#define __LINEAR_ALTITUDE_SPRINGS_S3D__

#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Forces/LINEAR_ALTITUDE_SPRINGS.h>
namespace PhysBAM{

template<class T_input>
class LINEAR_ALTITUDE_SPRINGS_S3D:public LINEAR_ALTITUDE_SPRINGS<VECTOR<T_input,3>,2>
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
    typedef LINEAR_ALTITUDE_SPRINGS<TV,2> BASE;typedef typename BASE::SPRING_STATE SPRING_STATE;typedef typename BASE::SPRING_PARAMETER SPRING_PARAMETER;
    using BASE::use_springs_compressed_beyond_threshold;using BASE::spring_compression_fraction_threshold;using BASE::print_number_used;
protected:
    using BASE::force_elements;using BASE::cfl_number;
public:
    using BASE::particles;using BASE::Invalidate_CFL;using BASE::max_strain_per_time_step;using BASE::use_rest_state_for_strain_rate;
    using BASE::spring_states;using BASE::use_plasticity;using BASE::mesh;
    using BASE::cache_strain;using BASE::strains_of_spring;
    typedef typename BASE::ELEMENT_ITERATOR ELEMENT_ITERATOR;
public:
    using BASE::parameters;using BASE::Compute_Plasticity;

    ARRAY<bool>* triangle_inverted;

    LINEAR_ALTITUDE_SPRINGS_S3D(DEFORMABLE_PARTICLES<TV>& particles,TRIANGLE_MESH& mesh)
        :LINEAR_ALTITUDE_SPRINGS<TV,2>(particles,mesh),triangle_inverted(0)
    {}

    virtual ~LINEAR_ALTITUDE_SPRINGS_S3D()
    {}

//#####################################################################
    void Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient); // assumes mass and restlength are already defined
    void Set_Restlength_From_Particles();
    void Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<const TV> X);
    void Set_Overdamping_Fraction(const T overdamping_fraction); // 1 is critically damped
    void Set_Overdamping_Fraction(const ARRAY<VECTOR<T,3> >& overdamping_fraction); // 1 is critically damped
    void Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction=1); // 1 is critically damped
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    T CFL_Strain_Rate() const override;

    T Potential_Energy(const int t,const T time) const;
    T Potential_Energy(const T time) const override;
//#####################################################################
};

template<class T> LINEAR_ALTITUDE_SPRINGS_S3D<T>*
Create_Altitude_Springs(DEFORMABLE_PARTICLES<VECTOR<T,3> >& particles,TRIANGLE_MESH& mesh,const T stiffness=4/(1+sqrt((T)2)),
    const T overdamping_fraction=2,const bool use_compressed_by_threshold_only=true,const T fraction_compression=.1,const bool limit_time_step_by_strain_rate=true,
    const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,const T restlength_enlargement_fraction=0,const bool verbose=true)
{
    return Create_Altitude_Springs_Base(particles,mesh,stiffness,overdamping_fraction,use_compressed_by_threshold_only,fraction_compression,
        limit_time_step_by_strain_rate,max_strain_per_time_step,use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose);
}

template<class T> LINEAR_ALTITUDE_SPRINGS_S3D<T>*
Create_Altitude_Springs(TRIANGULATED_SURFACE<T>& triangulated_surface,
    const T stiffness=4/(1+sqrt((T)2)),const T overdamping_fraction=2,const bool use_compressed_by_threshold_only=true,const T fraction_compression=.1,
    const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,const T restlength_enlargement_fraction=0,
    const bool verbose=true)
{
    return Create_Altitude_Springs(dynamic_cast<DEFORMABLE_PARTICLES<VECTOR<T,3> >&>(triangulated_surface.particles),triangulated_surface.mesh,stiffness,overdamping_fraction,use_compressed_by_threshold_only,
        fraction_compression,limit_time_step_by_strain_rate,max_strain_per_time_step,use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose);
}

}
#endif
