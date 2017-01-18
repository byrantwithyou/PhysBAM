//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TEST_BASE
//#####################################################################
#include <Core/Utilities/TEST_BASE.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
TEST_BASE::
TEST_BASE(const std::string& suite_name)
    :name(suite_name)
{
    Test_Registry().Set(name,this);
}
//#####################################################################
// Destructor
//#####################################################################
TEST_BASE::
~TEST_BASE()
{
}
//#####################################################################
// Function Report_Test_Result
//#####################################################################
void TEST_BASE::
Report_Test_Result(TEST_RESULT result,int n) const
{
    LOG::cout<<"Test "<<Identify_Test(n)<<(result==success?" PASSED.":result==failure?" FAILED.":" INCONCLUSIVE.")<<std::endl;
}
//#####################################################################
// Function Report_Out_Of_Range
//#####################################################################
void TEST_BASE::
Report_Out_Of_Range(int n) const
{
    LOG::cout<<"Test number "<<n<<" out of range.  Test suite "<<name<<" has "<<Number_Of_Tests()<<" tests."<<std::endl;
}
//#####################################################################
// Function Report_Registry_Failure
//#####################################################################
void TEST_BASE::
Report_Registry_Failure(std::string suite_name)
{
    LOG::cout<<"Test suite "<<suite_name<<" not registered."<<std::endl;
}
//#####################################################################
// Function Report_Result
//#####################################################################
void TEST_BASE::
Report_Result(TEST_RESULT result,const std::string& suite_name,int n,TEST_BASE* test)
{
    if(result==registry_failure) Report_Registry_Failure(suite_name);
    else if(result==out_of_range) test->Report_Out_Of_Range(n);
    else test->Report_Test_Result(result,n);
}
//#####################################################################
// Function Merge_Results
//#####################################################################
TEST_BASE::TEST_RESULT TEST_BASE::
Merge_Results(const TEST_RESULT a, const TEST_RESULT b)
{
    if(a==failure || b==failure) return failure;
    if(a==registry_failure || b==registry_failure) return registry_failure;
    if(a==out_of_range || b==out_of_range) return out_of_range;
    if(a==inconclusive || b==inconclusive) return inconclusive;
    return success;
}
//#####################################################################
// Function Identify_Test
//#####################################################################
std::string TEST_BASE::
Identify_Test(int n) const
{
    return LOG::sprintf("%s %i",name.c_str(),n);
}
//#####################################################################
// Function Lookup_Test_Suite
//#####################################################################
TEST_BASE* TEST_BASE::
Lookup_Test_Suite(std::string test_name)
{
    TEST_BASE* test;
    if(!Test_Registry().Get(test_name,test)) return 0;
    return test;
}
//#####################################################################
// Function Run_Suite_Test
//#####################################################################
TEST_BASE::TEST_RESULT TEST_BASE::
Run_Suite_Test(int n,bool verbose)
{
    if(n<0 || n>=Number_Of_Tests()){
        if(verbose) Report_Out_Of_Range(n);
        return out_of_range;}
    TEST_RESULT result=Run_Test(n);
    if(verbose) Report_Test_Result(result,n);
    return result;
}
//#####################################################################
// Function Run_All_Suite_Tests
//#####################################################################
TEST_BASE::TEST_RESULT TEST_BASE::
Run_All_Suite_Tests(bool verbose)
{
    TEST_RESULT result=success;
    for(int i=0;i<Number_Of_Tests();i++) result=Merge_Results(result,Run_Suite_Test(i,verbose));
    return result;
}
//#####################################################################
// Function Run_Test
//#####################################################################
TEST_BASE::TEST_RESULT TEST_BASE::
Run_Test(std::string test_name,int n,bool verbose)
{
    if(TEST_BASE* test=Lookup_Test_Suite(test_name)) return test->Run_Suite_Test(n,verbose);
    if(verbose) Report_Registry_Failure(test_name);
    return registry_failure;
}
//#####################################################################
// Function Run_All_Tests
//#####################################################################
TEST_BASE::TEST_RESULT TEST_BASE::
Run_All_Tests(std::string test_name,bool verbose)
{
    if(TEST_BASE* test=Lookup_Test_Suite(test_name)) return test->Run_All_Suite_Tests(verbose);
    if(verbose) Report_Registry_Failure(test_name);
    return registry_failure;
}
//#####################################################################
// Function Run_All_Tests
//#####################################################################
TEST_BASE::TEST_RESULT TEST_BASE::
Run_All_Tests(bool verbose)
{
    TEST_RESULT result=success;
    for(const auto& it:Test_Registry()) result=Merge_Results(result,it.data->Run_All_Suite_Tests(verbose));
    return result;
}
//#####################################################################
// Function Test_Registry
//#####################################################################
HASHTABLE<std::string,TEST_BASE*>& TEST_BASE::
Test_Registry()
{
    static HASHTABLE<std::string,TEST_BASE*> registry;
    return registry;
}
}
