//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class B_SPLINE_PATCH_CURVATURE_FORCE
//#####################################################################
#ifndef __B_SPLINE_PATCH_CURVATURE_FORCE__
#define __B_SPLINE_PATCH_CURVATURE_FORCE__

#include <Core/Data_Structures/FORCE_ELEMENTS.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/MATRIX_FORWARD.h>
#include <Tools/Tensors/TENSOR.h>
#include <Geometry/Topology_Based_Geometry/B_SPLINE_PATCH.h>
#include <Deformables/Constitutive_Models/MOONEY_RIVLIN_CURVATURE.h>
#include <Deformables/Forces/LAZY_HESSIAN_FORCE.h>
namespace PhysBAM{

template<class T,int gauss_order>
class B_SPLINE_PATCH_CURVATURE_FORCE:public LAZY_HESSIAN_FORCE<VECTOR<T,3> >
{
public:
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,2> IV;
    typedef DEFORMABLES_FORCES<TV> BASE;
    typedef MATRIX<T,3> TM;
    typedef VECTOR<VECTOR<T,3>,5> TM2;
    typedef MATRIX<MATRIX<T,3,3>,5,5> TT;

    using BASE::particles;

    const B_SPLINE_PATCH<TV,3>& spline;
    MOONEY_RIVLIN_CURVATURE<T> model;
    ARRAY<TV> X0;
    T density;
    bool recompute_hessian;

    struct DATA
    {
        MATRIX<VECTOR<TM,3>,gauss_order> G0_inv;
        MATRIX<VECTOR<T,3>,gauss_order> G0_det;
        MATRIX<VECTOR<VECTOR<T,16>,5>,gauss_order> A; // The 5 indexes through dw(0),dw(1),ddw(0,0),ddw(0,1),ddw(1,1).
        MATRIX<T,3,16> ge;
        MATRIX<TT,gauss_order> he;
    };

    T pe;
    ARRAY<DATA> data; 

    B_SPLINE_PATCH_CURVATURE_FORCE(DEFORMABLE_PARTICLES<TV>& particles,const B_SPLINE_PATCH<TV,3>& spline_input,const MOONEY_RIVLIN_CURVATURE<T>& model_input, T density);

    virtual ~B_SPLINE_PATCH_CURVATURE_FORCE();

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA> frequency) override;
    T CFL_Strain_Rate() const override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    T Potential_Energy(const T time) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    void Need_To_Recompute_Hessian(bool) override;
//#####################################################################
};
}
#endif
