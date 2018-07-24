//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PARSE_DATA_FEM__
#define __PARSE_DATA_FEM__
#include <Core/Vectors/VECTOR.h>
#include <functional>

namespace PhysBAM{
enum BC_TYPE {nobc,dirichlet_v,traction};

template<class TV>
struct PARSE_DATA_FEM
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    struct VERTEX_DATA
    {
        TV pt;
        BC_TYPE bc_type;
    };

    ARRAY<VERTEX_DATA> pts;
    ARRAY<VECTOR<int,2> > pipes;
    int half_width;
    T unit_length;

    void Parse_Input(const std::string& pipe_file);    
};

}
#endif
