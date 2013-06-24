//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ROTATION_SPRINGS
//#####################################################################
#ifndef __ROTATION_SPRINGS__
#define __ROTATION_SPRINGS__

#include <Tools/Matrices/MATRIX_POLICY.h>
#include <Tools/Vectors/VECTOR.h>
#include <Rigids/Forces_And_Torques/RIGIDS_FORCES.h>
namespace PhysBAM{

template<class TV> class MPI_SOLIDS;
template<class TV>
class ROTATION_SPRINGS:public RIGIDS_FORCES<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    enum WORKAROUND {d=TV::m};
    STATIC_ASSERT(d==3); // 2D not implemented yet
public:
    typedef RIGIDS_FORCES<TV> BASE;
    using BASE::rigid_body_collection;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;
protected:
    MPI_SOLIDS<TV>* mpi_solids;
public:
    SEGMENT_MESH& mesh;
    ARRAY<VECTOR<FRAME<TV>,2> > object_to_joint_frames;
    ARRAY<T_SPIN> angle_limits;
    ARRAY<DIAGONAL_MATRIX<T,T_SPIN::m> > stiffness; // torque / strain (past angle limits)
    ARRAY<T> damping; // torque / angular velocity
private:
    struct STATE
    {
        int edge;
        T_SPIN torque;
    };
    ARRAY<STATE> states;
public: 
    
    ROTATION_SPRINGS(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,SEGMENT_MESH& mesh);
    virtual ~ROTATION_SPRINGS();

    // unused callbacks (TODO)
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE {}
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> rigid_frequency)
        PHYSBAM_OVERRIDE {} // TODO: write correct CFL
    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const {PHYSBAM_NOT_IMPLEMENTED();}
    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE {}
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE {return FLT_MAX;}

//#####################################################################
    void Update_Position_Based_State(const T time);
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
