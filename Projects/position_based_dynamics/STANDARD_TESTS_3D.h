//#####################################################################
// Copyright 3015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STANDARD_TESTS_3D__
#define __STANDARD_TESTS_3D__

#include <Tools/Parsing/PARSE_ARGS.h>
#include "PBD_EXAMPLE.h"
#include "STANDARD_TESTS_BASE.h"

namespace PhysBAM{

template<class TV> class STANDARD_TESTS;

template<class T>
class STANDARD_TESTS<VECTOR<T,3> >:public STANDARD_TESTS_BASE<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef STANDARD_TESTS_BASE<TV> BASE;

public:
    using BASE::test_number;

    STANDARD_TESTS(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args);
    virtual ~STANDARD_TESTS();

    void Write_Output_Files(const int frame) PHYSBAM_OVERRIDE;
    void Read_Output_Files(const int frame) PHYSBAM_OVERRIDE;
    void Initialize() PHYSBAM_OVERRIDE;
    void Begin_Frame(const int frame) PHYSBAM_OVERRIDE;
    void End_Frame(const int frame) PHYSBAM_OVERRIDE;
    void Begin_Time_Step(const T time) PHYSBAM_OVERRIDE;
    void End_Time_Step(const T time) PHYSBAM_OVERRIDE;

//#####################################################################
};
}

#endif
