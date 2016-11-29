//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_SPRINGS
//#####################################################################
#ifndef __LINEAR_SPRINGS__
#define __LINEAR_SPRINGS__

#include <Core/Data_Structures/FORCE_ELEMENTS.h>
#include <Core/Matrices/MATRIX_FORWARD.h>
#include <Core/Vectors/VECTOR_2D.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class LINEAR_SPRINGS:public DEFORMABLES_FORCES<TV>,public SPRINGS_TAG
{
    typedef typename TV::SCALAR T;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;using BASE::Invalidate_CFL;using BASE::cfl_number;
    using BASE::max_strain_per_time_step;using BASE::use_rest_state_for_strain_rate;using BASE::compute_half_forces;
    typedef typename FORCE_ELEMENTS::ITERATOR SEGMENT_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;
    typedef MATRIX<T,TV::m> T_MATRIX;

public:
    SEGMENT_MESH& segment_mesh;
    ARRAY<T> youngs_modulus; // units of force (i.e. force per unit strain)
    T constant_youngs_modulus; // units of force (i.e. force per unit strain)
    ARRAY<T> restlength,visual_restlength; // visual restlength corresponds to length between particles; restlength may be larger than this to avoid zero/small restlength
    ARRAY<T> damping; // units of force*time (i.e. force per unit strain rate)
    T constant_damping; // units of force*time (i.e. force per unit strain rate)
    bool use_plasticity;
    ARRAY<T> plastic_yield_strain;
    ARRAY<T> plastic_hardening;
    ARRAY<T> plastic_visual_restlength;
    T plasticity_clamp_ratio;
    bool cache_strain;
    mutable ARRAY<VECTOR<T,2> > strains_of_segment; // VECTOR<T,2>(strain_rate, strain)

    bool verbose;

protected:
    struct STATE{
        STATE()
            :coefficient(0),sqrt_coefficient(0)
        {}

        VECTOR<int,2> nodes; // copy of nodes to reduce associativity needed in cache
        T coefficient;
        T sqrt_coefficient;
        TV direction;
        TV dX;
    };
    ARRAY<STATE> states;
    ARRAY<T> current_lengths;

public:
    FORCE_ELEMENTS force_segments;
    LINEAR_SPRINGS(DEFORMABLE_PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh_input);

    virtual ~LINEAR_SPRINGS();

//#####################################################################
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    virtual void Set_Stiffness(const T youngs_modulus_input);
    virtual void Set_Stiffness(ARRAY_VIEW<const T> youngs_modulus_input);
    virtual void Set_Restlength(ARRAY_VIEW<const T> restlength_input);
    virtual void Set_Damping(const T damping_input);
    virtual void Set_Damping(ARRAY_VIEW<const T> damping_input);
    template<class T_FIELD> void Enable_Plasticity(const T_FIELD& plastic_yield_strain_input,const T_FIELD& plastic_hardening_input,const T plasticity_clamp_ratio_input=4);
    virtual void Set_Restlength_From_Particles();
    virtual void Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<TV> material_coordinates);
    virtual void Clamp_Restlength(const T clamped_restlength);
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    virtual void Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient); // assumes mass and restlength are already defined
    virtual void Set_Stiffness_Based_On_Reduced_Mass(ARRAY_VIEW<const T> scaling_coefficient); // assumes mass and restlength are already defined
    virtual void Set_Overdamping_Fraction(const T overdamping_fraction=1); // 1 is critically damped
    virtual void Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction); // 1 is critically damped
    void Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction=1); // 1 is critically damped
    void Clamp_Restlength_With_Fraction_Of_Springs(const T fraction=.01);
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    int Velocity_Dependent_Forces_Size() const override;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const override;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) override;
    T CFL_Strain_Rate() const override;
    T Average_Restlength() const;
    void Print_Restlength_Statistics() const;
    void Print_Deformation_Statistics() const;
    T Maximum_Compression_Or_Expansion_Fraction(int* index=0) const;
    T Potential_Energy(const int s,const T time) const;
    T Potential_Energy(const T time) const override;
    void Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name="") const override;

    TV Endpoint_Velocity(int s,int b) const;
    T Endpoint_Kinetic_Energy(int s,int b) const;
    T Endpoint_Kinetic_Energy(int s) const;
    typename TV::SCALAR Effective_Impulse_Factor(int s) const;
//#####################################################################
};

template<class TV> LINEAR_SPRINGS<TV>*
Create_Edge_Springs(DEFORMABLE_PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh,const typename TV::SCALAR stiffness=2e3,
    const typename TV::SCALAR overdamping_fraction=1,const bool limit_time_step_by_strain_rate=true,const typename TV::SCALAR max_strain_per_time_step=.1,
    const bool use_rest_state_for_strain_rate=true,const typename TV::SCALAR restlength_enlargement_fraction=0,const bool verbose=true);

template<class T_OBJECT> LINEAR_SPRINGS<typename T_OBJECT::VECTOR_T>*
Create_Edge_Springs(T_OBJECT& object,
    const typename T_OBJECT::SCALAR stiffness=2e3,const typename T_OBJECT::SCALAR overdamping_fraction=1,const bool limit_time_step_by_strain_rate=true,
    const typename T_OBJECT::SCALAR max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,const typename T_OBJECT::SCALAR restlength_enlargement_fraction=0,
    const bool verbose=true);

}
#endif
