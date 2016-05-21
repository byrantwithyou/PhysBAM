//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TEST_BASE
//#####################################################################
#ifndef __TEST_BASE__
#define __TEST_BASE__

#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Log/LOG.h>
#include <string>
namespace PhysBAM{

class TEST_BASE
{
public:
    enum TEST_RESULT{success,failure,out_of_range,registry_failure,inconclusive};

    std::string name;

    TEST_BASE(const std::string& suite_name)
        :name(suite_name)
    {
        Test_Registry().Set(name,this);
    }

    virtual ~TEST_BASE(){}

    virtual void Report_Test_Result(TEST_RESULT result,int n) const
    {LOG::cout<<"Test "<<Identify_Test(n)<<(result==success?" PASSED.":result==failure?" FAILED.":" INCONCLUSIVE.")<<std::endl;}

    virtual void Report_Out_Of_Range(int n) const
    {LOG::cout<<"Test number "<<n<<" out of range.  Test suite "<<name<<" has "<<Number_Of_Tests()<<" tests."<<std::endl;}

    static void Report_Registry_Failure(std::string suite_name)
    {LOG::cout<<"Test suite "<<suite_name<<" not registered."<<std::endl;}

    static void Report_Result(TEST_RESULT result,const std::string& suite_name,int n,TEST_BASE* test)
    {if(result==registry_failure) Report_Registry_Failure(suite_name);
    else if(result==out_of_range) test->Report_Out_Of_Range(n);
    else test->Report_Test_Result(result,n);}

    static TEST_RESULT Merge_Results(const TEST_RESULT a, const TEST_RESULT b)
    {if(a==failure || b==failure) return failure;
    if(a==registry_failure || b==registry_failure) return registry_failure;
    if(a==out_of_range || b==out_of_range) return out_of_range;
    if(a==inconclusive || b==inconclusive) return inconclusive;
    return success;}

    virtual TEST_RESULT Run_Test(int n)=0;

    virtual int Number_Of_Tests() const=0;

    virtual std::string Identify_Test(int n) const
    {return LOG::sprintf("%s %i",name.c_str(),n);}

    static TEST_BASE* Lookup_Test_Suite(std::string test_name)
    {TEST_BASE* test;if(!Test_Registry().Get(test_name,test)) return 0;return test;}

    TEST_RESULT Run_Suite_Test(int n,bool verbose=true)
    {if(n<0 || n>=Number_Of_Tests()){if(verbose) Report_Out_Of_Range(n);return out_of_range;}TEST_RESULT result=Run_Test(n);if(verbose) Report_Test_Result(result,n);return result;}

    TEST_RESULT Run_All_Suite_Tests(bool verbose=true)
    {TEST_RESULT result=success;for(int i=0;i<Number_Of_Tests();i++) result=Merge_Results(result,Run_Suite_Test(i,verbose));return result;}

    static TEST_RESULT Run_Test(std::string test_name,int n,bool verbose=true)
    {if(TEST_BASE* test=Lookup_Test_Suite(test_name)) return test->Run_Suite_Test(n,verbose);if(verbose) Report_Registry_Failure(test_name);return registry_failure;}

    static TEST_RESULT Run_All_Tests(std::string test_name,bool verbose=true)
    {if(TEST_BASE* test=Lookup_Test_Suite(test_name)) return test->Run_All_Suite_Tests(verbose);if(verbose) Report_Registry_Failure(test_name);return registry_failure;}

    static TEST_RESULT Run_All_Tests(bool verbose=true)
    {TEST_RESULT result=success;for(const auto& it:Test_Registry()) result=Merge_Results(result,it.data->Run_All_Suite_Tests(verbose));return result;}

//#####################################################################
    static HASHTABLE<std::string,TEST_BASE*>& Test_Registry();
//#####################################################################
};
}
#endif
