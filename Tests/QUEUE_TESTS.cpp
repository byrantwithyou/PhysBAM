//#####################################################################
// Copyright 2012, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUEUE_TESTS
//#####################################################################

#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Tools/Utilities/TEST_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

class QUEUE_TESTS:public TEST_BASE
{
    bool Test(const bool test,const char* str_test,int line,bool& ok) const
    {if(test) return true;LOG::cout<<"FAILED: "<<str_test<<"    line "<<line<<std::endl;ok=false;return false;}

public:
    QUEUE_TESTS()
        :TEST_BASE("queue")
    {}

    virtual ~QUEUE_TESTS(){}

#define TEST(x) Test((x),#x,__LINE__,ok)

    TEST_RESULT Run_Test(int n)
    {
        bool ok=true;
        QUEUE<int> q(3);

        TEST(q.Size()==0);
        TEST(q.Empty()==true);
        TEST(q.Full()==false);

        q.Enqueue(3);
        TEST(q.Empty()==false);
        TEST(q.Dequeue()==3);
        TEST(q.Empty()==true);
        TEST(q.Size()==0);
        q.Enqueue(3);
        q.Enqueue(7);
        q.Enqueue(11);
        TEST(q(2)==11);
        TEST(q.Size()==3);
        TEST(q.Peek()==3);
        TEST(q.Dequeue()==3);
        TEST(q.Peek()==7);
        TEST(q.Size()==2);

        q.Enqueue(15);
        TEST(q.Size()==3);
        q.Safe_Enqueue(19);
        TEST(q.Size()==4);
        TEST(q.Peek()==7);

        TEST(q(0)==7);
        TEST(q(1)==11);
        TEST(q(2)==15);
        TEST(q(3)==19);

        TEST(q.Dequeue()==7);
        TEST(q.Dequeue()==11);
        TEST(q.Dequeue()==15);
        TEST(q(0)==19);
        TEST(q.Peek()==19);
        TEST(q.Dequeue()==19);

        q.Remove_All();
        TEST(q.Size()==0);

        q.Safe_Enqueue(3);
        TEST(q.Empty()==false);
        q.Safe_Enqueue(7);
        q.Safe_Enqueue(11);
        TEST(q.Size()==3);
        q.Enqueue(44);
        
        TEST(q.Dequeue()==3);
        TEST(q.Dequeue()==7);
        q.Enqueue(15);
        TEST(q.Peek()==11);
        TEST(q(2)==15);
        TEST(q.Size()==3);

        TEST(q.Empty()==false);
        
        return ok?success:failure;
    }

    int Number_Of_Tests() const
    {return 1;}

//#####################################################################
};
static QUEUE_TESTS queue_tests;
}
