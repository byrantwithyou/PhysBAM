//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TEST_BASE
//#####################################################################
#ifndef __TEST_BASE__
#define __TEST_BASE__

#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/LOG.h>
#include <string>
namespace PhysBAM{

class TEST_BASE
{
public:
    enum TEST_RESULT{success,failure,out_of_range,registry_failure,inconclusive};

    std::string name;

    TEST_BASE(const std::string& suite_name);

    virtual ~TEST_BASE();

//#####################################################################
    virtual void Report_Test_Result(TEST_RESULT result,int n) const;
    virtual void Report_Out_Of_Range(int n) const;
    static void Report_Registry_Failure(std::string suite_name);
    static void Report_Result(TEST_RESULT result,const std::string& suite_name,int n,TEST_BASE* test);
    static TEST_RESULT Merge_Results(const TEST_RESULT a, const TEST_RESULT b);
    virtual TEST_RESULT Run_Test(int n)=0;
    virtual int Number_Of_Tests() const=0;
    virtual std::string Identify_Test(int n) const;
    static TEST_BASE* Lookup_Test_Suite(std::string test_name);
    TEST_RESULT Run_Suite_Test(int n,bool verbose=true);
    TEST_RESULT Run_All_Suite_Tests(bool verbose=true);
    static TEST_RESULT Run_Test(std::string test_name,int n,bool verbose=true);
    static TEST_RESULT Run_All_Tests(std::string test_name,bool verbose=true);
    static TEST_RESULT Run_All_Tests(bool verbose=true);
    static HASHTABLE<std::string,TEST_BASE*>& Test_Registry();
//#####################################################################
};
}
#endif
