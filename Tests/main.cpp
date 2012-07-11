//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Utilities/TEST_BASE.h>
#include <cstdlib>
using namespace PhysBAM; 

int main(int argc,char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Set_Extra_Arguments(2,STRING_UTILITIES::string_sprintf("%s <suite-name> <test-number>",argv[0]));
    parse_args.Parse();
    parse_args.Print_Arguments();

    TEST_BASE::TEST_RESULT result;
    std::string test_suite=parse_args.Extra_Arg(0),test_number=parse_args.Extra_Arg(1).c_str();
    if(test_suite=="all" && test_number=="all") result=TEST_BASE::Run_All_Tests();
    else if(test_number=="all") result=TEST_BASE::Run_All_Tests(test_suite);
    else result=TEST_BASE::Run_Test(test_suite,atoi(test_number.c_str()));

    return result==TEST_BASE::success?EXIT_SUCCESS:EXIT_FAILURE;
}
//#####################################################################
