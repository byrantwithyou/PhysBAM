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
template<class TV> class PENALTY_FORCE_COLLECTION;
template<class TV> class MOVE_RIGID_BODY_DIFF;
template<class TV> class TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY;
template<class TV> class TRIANGLE_COLLISIONS;

template<class TV>
class STANDARD_TESTS_BASE:public SOLIDS_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::data_directory;using BASE::solid_body_collection;
    using BASE::solids_parameters;using BASE::solids_evolution;
    using BASE::m;using BASE::s;using BASE::kg;
    using BASE::viewer_dir;using BASE::stream_type;

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
    bool use_rd=false,use_rr=false,use_dd=false,use_di=false;
    bool use_rd_k=false,use_rr_k=false,use_dd_k=false,use_di_k=false;
    T rd_k=0,rr_k=0,dd_k=0,di_k=0;
    bool use_rd_mu=false,use_rr_mu=false,use_dd_mu=false,use_di_mu=false;
    T rd_mu=0,rr_mu=0,dd_mu=0,di_mu=0;
    bool use_rd_ccd=false,use_rr_ccd=false,use_di_ccd=false;
    PENALTY_FORCE_COLLECTION<TV>* pfd=0;
    ARRAY<MOVE_RIGID_BODY_DIFF<TV> > move_rb_diff;
    GRID<TV> detection_grid;
    bool use_bisection=false;
    
    STANDARD_TESTS_BASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS_BASE();

    void After_Get_Initial_Data(bool automatically_add_to_collision_structures);
    void After_Initialize_Bodies();
    void Add_Collision_Object(IMPLICIT_OBJECT<TV>* io);
    void Preprocess_Substep(const T dt,const T time) override;
    void Postprocess_Substep(const T dt,const T time) override;
    void Read_Output_Files_Solids() override;
    void Write_Output_Files() const override;
    void Init_Penalty_Collection();
};
}
#endif
