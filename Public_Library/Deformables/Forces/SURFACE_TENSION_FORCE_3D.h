//#####################################################################
// Copyright 2015.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SURFACE_TENSION_FORCE_3D
//#####################################################################
#ifndef __SURFACE_TENSION_FORCE_3D__
#define __SURFACE_TENSION_FORCE_3D__

#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class SURFACE_TENSION_FORCE_3D:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,3> TV_INT;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    TRIANGULATED_SURFACE<T>& surface;
    T surface_tension_coefficient;
    ARRAY<T> areas;
    ARRAY<TV> normals;
    T dt;
    bool apply_explicit_forces;

    SURFACE_TENSION_FORCE_3D(TRIANGULATED_SURFACE<T>& surface_input,T surface_tension_coefficient_input);
    virtual ~SURFACE_TENSION_FORCE_3D();

//#####################################################################
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    int Velocity_Dependent_Forces_Size() const override;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const override;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const override;
    void Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency) override;
    T CFL_Strain_Rate() const override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    T Potential_Energy(const T time) const override;
//#####################################################################
};
}
#endif
