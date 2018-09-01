//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PARSE_DATA_FEM__
#define __PARSE_DATA_FEM__
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Vectors/VECTOR.h>
#include "COMMON.h"

namespace PhysBAM{
enum BC_TYPE {analytic,dirichlet_v,traction};
enum JOINT_TYPE {default_joint,corner_joint};

template<class TV> struct ANALYTIC_VECTOR;
template<class TV> struct ANALYTIC_SCALAR;

template<class TV>
struct PARSE_DATA_FEM
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    struct VERTEX_DATA
    {
        TV pt;
        BC_ID bc_id=BC_ID(0); // wall
        JOINT_TYPE joint_type=default_joint;
        ARRAY<PIPE_ID> joints;
    };

    struct BC_FUNC
    {
        BC_TYPE type=dirichlet_v;
        ANALYTIC_VECTOR<TV>* velocity=0;
        ANALYTIC_SCALAR<TV>* pressure=0;
        ANALYTIC_VECTOR<TV>* traction=0;
    };
    ANALYTIC_VECTOR<TV>* force=0;
    BC_ID analytic_bc=BC_ID(-1),wall_bc=BC_ID(0);
    ARRAY<BC_FUNC,BC_ID> bc;

    TV Velocity(const TV& X,BC_ID bc_id) const;
    TV Traction(const TV& X,const TV& N,T mu,BC_ID bc_id) const;
    TV Force(const TV& X,T mu) const;
    T Divergence(const TV& X) const;
    T Pressure(const TV& X) const;

    ARRAY<VERTEX_DATA,VERTEX_ID> pts;
    ARRAY<VECTOR<VERTEX_ID,2>,PIPE_ID> pipes;
    int half_width;
    T unit_length;

    PARSE_DATA_FEM();
    ~PARSE_DATA_FEM();
    
    void Parse_Input(const std::string& pipe_file);    
};

}
#endif
