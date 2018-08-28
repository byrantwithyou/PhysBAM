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
enum BC_TYPE {nobc,dirichlet_v,traction};
enum JOINT_TYPE {default_joint,corner_joint};

template<class TV>
struct PARSE_DATA_FEM
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    struct VERTEX_DATA
    {
        TV pt;
        BC_TYPE bc_type;
        TV bc;
        JOINT_TYPE joint_type;
        ARRAY<PIPE_ID> joints;
    };

    ARRAY<VERTEX_DATA,VERTEX_ID> pts;
    ARRAY<VECTOR<VERTEX_ID,2>,PIPE_ID> pipes;
    int half_width;
    T unit_length;

    void Parse_Input(const std::string& pipe_file);    
};

}
#endif
