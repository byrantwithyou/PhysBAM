//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STANDARD_TESTS_BASE__
#define __STANDARD_TESTS_BASE__

#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <fstream>
#ifdef USE_OPENMP
#include <omp.h>
#endif
namespace PhysBAM{

template<class TV> class BACKWARD_EULER_EVOLUTION;
template<class TV> class RIGID_DEFORMABLE_PENALTY_WITH_FRICTION;
template<class TV> class RIGID_PENALTY_WITH_FRICTION;
template<class TV> class IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION;
template<class TV> class SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION;
template<class TV> class MOVE_RIGID_BODY_DIFF;

template<class TV>
class STANDARD_TESTS_BASE:public SOLIDS_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,3> TV_INT;
public:
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::data_directory;using BASE::solid_body_collection;
    using BASE::solids_parameters;using BASE::solids_evolution;
    using BASE::m;using BASE::s;using BASE::kg;

    SOLIDS_STANDARD_TESTS<TV> tests;

    bool test_forces=false;
    bool use_newmark=false,use_newmark_be=false;
    BACKWARD_EULER_EVOLUTION<TV>* backward_euler_evolution;
    bool no_line_search=false;
    bool no_descent=false;
    T unit_rho=1,unit_p=1,unit_N=1,unit_J=1;
    bool use_vanilla_newton=false;
    int threads=1;
    T rd_penalty_stiffness=0;
    T rd_penalty_friction=0.3;
    bool use_rd_penalty=false;
    RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>* rd_penalty=0;
    RIGID_PENALTY_WITH_FRICTION<TV>* rr_penalty=0;
    IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>* di_penalty=0;
    SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>* dd_penalty=0;
    ARRAY<MOVE_RIGID_BODY_DIFF<TV> > move_rb_diff;

    STANDARD_TESTS_BASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS_BASE();

    void After_Get_Initial_Data(bool automatically_add_to_collision_structures);
    void After_Initialize_Bodies();
    void Get_RD_Collision_Candidates();
    void Get_DD_Collision_Candidates();
    void Get_DI_Collision_Candidates();
    void Get_RR_Collision_Candidates();
    void Preprocess_Substep(const T dt,const T time) override;
    void Postprocess_Substep(const T dt,const T time) override;
};
}
#endif
