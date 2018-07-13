//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PARSE__
#define __PARSE__
#include <Core/Math_Tools/RANGE.h>
#include <Core/Vectors/VECTOR.h>
#include <functional>

namespace PhysBAM{
enum INDEX_TYPE {fluid, wall, dirichlet, nodof};

template<class TV>
struct PARSE_DATA
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    struct VERTEX_DATA
    {
        typedef typename TV::SCALAR T;
        VECTOR<int,TV::m> pt;
        INDEX_TYPE bc_type;
        int bc_side;
        typename TV::SCALAR bc_value;
        RANGE<VECTOR<int,TV::m> > box;
        int connected_sides;
        // std::function<T(VECTOR<T,TV::m-1>)> bc_func_value;

        // T Value(const VECTOR<int,TV::m>& c) const;
        // T Value(const FACE_INDEX<TV::m>& f) const;
    };

    ARRAY<VERTEX_DATA> pts;
    ARRAY<VECTOR<int,2> > pipes;
    TV_INT box_size;
    int half_width;

    void Parse_Input(const std::string& pipe_file);    
};

}
#endif
