//#####################################################################
// Copyright 2012, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STACK_TESTS
//#####################################################################

#include <Tools/Data_Structures/STACK.h>
#include <Tools/Utilities/TEST_BASE.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{

class STACK_TESTS:public TEST_BASE
{
    bool Test(const bool test,const char* str_test,int line,bool& ok) const
    {if(test) return true;LOG::cout<<"FAILED: "<<str_test<<"    line "<<line<<std::endl;ok=false;return false;}

public:
    STACK_TESTS()
        :TEST_BASE("stack")
    {}

    virtual ~STACK_TESTS(){}

#define TEST(x) Test((x),#x,__LINE__,ok)

    TEST_RESULT Run_Test(int n)
    {
        bool ok=true;
        STACK<int> q;

       TEST(q.Empty()==true);

        q.Push(3);
        q.Push(5);
        TEST(q.Pop()==5);
        TEST(q.Empty()==false);

        q.Remove_All();
        TEST(q.Empty()==true);
        
        STACK<int> s;
        s.Preallocate(2);
        TEST(s.Empty()==true);
        s.Push(5);
        s.Push(9);
        TEST(s.Peek()==9);

        s.Clean_Memory();
        TEST(s.Empty()==true);
        
        return ok?success:failure;
    }

    int Number_Of_Tests() const
    {return 1;}

//#####################################################################
};
static STACK_TESTS stack_tests;
}
