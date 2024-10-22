//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_FORCES
//#####################################################################
#ifndef __BW_FORCES__
#define __BW_FORCES__

#include <Core/Matrices/MATRIX_FORWARD.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

class TRIANGLE_MESH;

template<class TV,int d,int m>
class BW_FORCES:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

public:
    TRIANGLE_MESH& triangle_mesh;

protected:
    struct STATE{
        STATE()
            :stiffness_coefficient(.5),damping_coefficient(0)
        {}

        // These don't change across the entire sim
        VECTOR<int,m> nodes;
        T stiffness_coefficient;
        T damping_coefficient;

        // These are set each time UPBS is called
        VECTOR<T,d> C;
        VECTOR<T,d> C_dot;
        VECTOR<MATRIX<T,3,d>,m> dC_dx;
        MATRIX<MATRIX<T,3>,m> dC_dxi_dxj_times_C;
        MATRIX<MATRIX<T,3>,m> dC_dxi_dxj_times_C_dot;
    };
    ARRAY<STATE> states;
public:
    ARRAY<int> force_simplices;
    BW_FORCES(DEFORMABLE_PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh_input,const T stiffness_coefficient_input,const T damping_coefficient_input);

    virtual ~BW_FORCES();

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose=false) const override;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) override;
    T CFL_Strain_Rate() const override;
//#####################################################################
};
}
#endif
