//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BEZIER_C2_FORCE
//#####################################################################
#ifndef __BEZIER_C2_FORCE__
#define __BEZIER_C2_FORCE__

#include <Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Geometry/Topology_Based_Geometry/BEZIER_SPLINE.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class BEZIER_C2_FORCE:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;

    const BEZIER_SPLINE<TV,3>& spline;
    ARRAY<TV> X0;
    T stiffness;

    struct DATA
    {
        VECTOR<TV,4> ge;
        MATRIX<MATRIX<T,TV::m>,4> he;
        VECTOR<int,4> pts;
    };

    VECTOR<TV,3> end_ge[2];
    MATRIX<MATRIX<T,TV::m>,3> end_he[2];
    VECTOR<int,3> end_pts[2];

    T pe;
    ARRAY<DATA> data;

    BEZIER_C2_FORCE(DEFORMABLE_PARTICLES<TV>& particles,const BEZIER_SPLINE<TV,3>& spline_input,
        T stiffness_input);

    virtual ~BEZIER_C2_FORCE();

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) PHYSBAM_OVERRIDE;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
