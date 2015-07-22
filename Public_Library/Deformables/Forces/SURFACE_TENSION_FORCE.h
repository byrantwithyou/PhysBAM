//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SURFACE_TENSION_FORCE
//#####################################################################
#ifndef __SURFACE_TENSION_FORCE__
#define __SURFACE_TENSION_FORCE__

#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class SURFACE_TENSION_FORCE:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    SEGMENTED_CURVE_2D<T>& surface;
    T surface_tension_coefficient;
    ARRAY<T> coefficients;
    ARRAY<TV> normal;
    ARRAY<MATRIX<T,TV::m,TV::m-1> > tangential;
    ARRAY<T> sqrt_coefficients;
    T dt;
    bool apply_explicit_forces;
    bool apply_implicit_forces;

    SURFACE_TENSION_FORCE(SEGMENTED_CURVE_2D<T>& surface_input,T surface_tension_coefficient_input);
    virtual ~SURFACE_TENSION_FORCE();

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
    void Dump_Curvatures() const;
//#####################################################################
};
}
#endif
