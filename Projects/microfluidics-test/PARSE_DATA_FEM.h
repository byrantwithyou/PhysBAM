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
enum BC_TYPE {dirichlet_v,traction};
enum JOINT_TYPE {default_joint,corner_joint};

template<class TV> struct ANALYTIC_VECTOR;
template<class TV> struct ANALYTIC_SCALAR;

template<class TV,class TW>
struct PARSE_DATA_FEM
{
    typedef typename TV::SCALAR T;

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
        T flowrate;
        TW pt,dir;
        TW traction;
    };
    ANALYTIC_VECTOR<TW>* force=0;
    BC_ID wall_bc=BC_ID(0);
    ARRAY<BC_FUNC,BC_ID> bc;
    ANALYTIC_VECTOR<TW>* analytic_velocity=0;
    ANALYTIC_SCALAR<TW>* analytic_pressure=0;

    TW Velocity(const TW& X,BC_ID bc_id) const;
    TW Traction(const TW& X,const TW& N,T mu,BC_ID bc_id) const;
    TW Force(const TW& X,T mu) const;
    T Divergence(const TW& X) const;
    T Pressure(const TW& X) const;

    ARRAY<VERTEX_DATA,VERTEX_ID> pts;
    ARRAY<VECTOR<VERTEX_ID,2>,PIPE_ID> pipes;
    int half_width=1,height=0;
    T unit_length=0.1;

    PARSE_DATA_FEM();
    ~PARSE_DATA_FEM();
    
    void Parse_Input(const std::string& pipe_file);    
};

}
#endif
