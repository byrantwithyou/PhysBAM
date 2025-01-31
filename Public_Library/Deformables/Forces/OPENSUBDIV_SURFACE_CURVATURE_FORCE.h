//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENSUBDIV_SURFACE_CURVATURE_FORCE
//#####################################################################
#ifndef __OPENSUBDIV_SURFACE_CURVATURE_FORCE__
#define __OPENSUBDIV_SURFACE_CURVATURE_FORCE__

#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/MATRIX_FORWARD.h>
#include <Geometry/Topology_Based_Geometry/OPENSUBDIV_SURFACE.h>
#include <Deformables/Constitutive_Models/MOONEY_RIVLIN_CURVATURE.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class T,int gauss_order>
class OPENSUBDIV_SURFACE_CURVATURE_FORCE:public DEFORMABLES_FORCES<VECTOR<T,3> >
{
public:
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,2> IV;
    typedef DEFORMABLES_FORCES<TV> BASE;
    typedef MATRIX<T,3> TM;
    typedef VECTOR<VECTOR<T,3>,5> TM2;
    typedef MATRIX<MATRIX<T,3,3>,5,5> TT;

    using BASE::particles;

    const OPENSUBDIV_SURFACE<TV>& surf;
    MOONEY_RIVLIN_CURVATURE<T> model;

    struct DATA
    {
        ARRAY<VECTOR<T,3> > ge;
        MATRIX<TT,gauss_order> he;
    };

    T pe;
    ARRAY<DATA> data; 

    OPENSUBDIV_SURFACE_CURVATURE_FORCE(DEFORMABLE_PARTICLES<TV>& particles,const OPENSUBDIV_SURFACE<TV>& surf_input,const MOONEY_RIVLIN_CURVATURE<T>& model_input);

    virtual ~OPENSUBDIV_SURFACE_CURVATURE_FORCE();

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
