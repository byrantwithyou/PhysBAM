//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class B_SPLINE_CURVATURE_FORCE
//#####################################################################
#ifndef __B_SPLINE_CURVATURE_FORCE__
#define __B_SPLINE_CURVATURE_FORCE__

#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/MATRIX_FORWARD.h>
#include <Geometry/Topology_Based_Geometry/B_SPLINE.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class B_SPLINE_CURVATURE_FORCE:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;

    const B_SPLINE<TV,3>& spline;
    static const int gauss_order=7;
    ARRAY<TV> X0;
    T curvature_stiffness;
    T length_stiffness;

    struct DATA
    {
        VECTOR<T,gauss_order> kappa0,length0;
        VECTOR<TV,4> ge;
        MATRIX<MATRIX<T,TV::m>,4> he;
    };

    T pe;
    ARRAY<DATA> data;

    B_SPLINE_CURVATURE_FORCE(DEFORMABLE_PARTICLES<TV>& particles,const B_SPLINE<TV,3>& spline_input,
        T curvature_stiffness_input,T length_stiffness_input);

    virtual ~B_SPLINE_CURVATURE_FORCE();

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose=false) const override;
    void Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA> frequency) override;
    T CFL_Strain_Rate() const override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    T Potential_Energy(const T time) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
//#####################################################################
};
}
#endif
