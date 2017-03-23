//#####################################################################
// Copyright 2017, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//##################################################################### 
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Compressible/Shallow_Water_Equations/SHALLOW_WATER_STATE.h>
#include "COMMON_DATA.h"

namespace PhysBAM
{
template<class TV> void
Setup_Common(SHALLOW_WATER_STATE<TV>& st,PARSE_ARGS& parse_args,COMMON_DATA<TV>& cd);

template<class T> void
Setup_Example(SHALLOW_WATER_STATE<VECTOR<T,2> >& st,PARSE_ARGS& parse_args,
    COMMON_DATA<VECTOR<T,2> >& cd)
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
    typedef VECTOR<T,3> T_VEC;
    Setup_Common(st,parse_args,cd);
    parse_args.Parse();

    switch(cd.test_number)
    {
        case 0:{
            st.grid.Initialize(TV_INT()+cd.resolution,RANGE<TV>::Unit_Box(),true);

            st.initialize=[&st]()
                {
                    st.U.Fill(T_VEC(1,0,0));
                };
            break;
        }
    }
}
template void Setup_Example<double>(SHALLOW_WATER_STATE<VECTOR<double,2> >&,
    PARSE_ARGS&,COMMON_DATA<VECTOR<double,2> >&);
template void Setup_Example<float>(SHALLOW_WATER_STATE<VECTOR<float,2> >&,
    PARSE_ARGS&,COMMON_DATA<VECTOR<float,2> >&);

}
