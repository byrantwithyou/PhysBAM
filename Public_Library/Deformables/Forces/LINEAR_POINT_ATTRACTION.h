//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_POINT_ATTRACTION
//#####################################################################
#ifndef __LINEAR_POINT_ATTRACTION__
#define __LINEAR_POINT_ATTRACTION__

#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class LINEAR_POINT_ATTRACTION:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT T_MESH;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    T_MESH& surface;
    ARRAY<int> referenced_particles;
    T coefficient;
    TV point;
    T dt;
    bool apply_explicit_forces;
    bool apply_implicit_forces;

    LINEAR_POINT_ATTRACTION(T_MESH& mesh,const TV& pt,T coefficient_input);
    virtual ~LINEAR_POINT_ATTRACTION();

//#####################################################################
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    int Velocity_Dependent_Forces_Size() const override;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const override;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const override;
    void Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency) override;
    T CFL_Strain_Rate() const override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    T Potential_Energy(const T time) const override;
//#####################################################################
};
}
#endif
