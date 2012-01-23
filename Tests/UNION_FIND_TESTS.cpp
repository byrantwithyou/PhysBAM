//#####################################################################
// Copyright 2012, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNION_FIND_TESTS
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/CONSTANT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Utilities/TEST_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

class UNION_FIND_TESTS:public TEST_BASE
{
    bool Test(const bool test,const char* str_test,int line,bool& ok) const
    {if(test) return true;LOG::cout<<"FAILED: "<<str_test<<"    line "<<line<<std::endl;ok=false;return false;}

public:
    UNION_FIND_TESTS()
        :TEST_BASE("union_find")
    {}

    virtual ~UNION_FIND_TESTS(){}

#define TEST(x) Test((x),#x,__LINE__,ok)

    TEST_RESULT Run_Test(int n)
    {
        bool ok=true;
        
        UNION_FIND<int> u0;
        u0.Initialize(2);
        TEST(u0.Size()==2);
        TEST(u0.Add_Entry()==3);

        UNION_FIND<int> u1(4);
        u1.Initialize(4);
        TEST(u1.Size()==4);
        TEST(u1.Is_Root(0)==true);
        TEST(u1.Is_Root(1)==true);
        TEST(u1.Is_Root(2)==true);
        TEST(u1.Is_Root(3)==true);
 
        TEST(u1.Union(0,3)==0);
        TEST(u1.Is_Root(0)==true);
        TEST(u1.Is_Root(3)==false);

        TEST(u1.Find(3)==0);
        TEST(u1.Find(1)==1);
        TEST(u1.Find(2)==2);
        TEST(u1.Find(0)==0);

        u1.Clear_Connectivity();
        TEST(u1.Find(3)==3);
        VECTOR<int,3> inds;
        inds(0)=1;inds(1)=2;inds(2)=3;
        int r=u1.Union(inds);
        
        TEST(r==1||r==2||r==3);
        TEST(u1.Find(3)==r);
        TEST(u1.Find(1)==r);
        TEST(u1.Find(2)==r);
        TEST(u1.Find(0)==0);

        UNION_FIND<int> u2(6);
        u2.Initialize(6);

        inds(0)=0;inds(1)=2;inds(2)=4;
        int r1=u2.Union(inds);
        TEST(r1==0||r1==2||r1==4);

        inds(0)=1;inds(1)=3;inds(2)=5;
        int r2=u2.Union(inds);
        TEST(r2==1||r2==3||r2==5);

        TEST(u2.Find(0)==r1);
        TEST(u2.Find(1)==r2);
        TEST(u2.Find(2)==r1);
        TEST(u2.Find(3)==r2);
        TEST(u2.Find(4)==r1);
        TEST(u2.Find(5)==r2);

        TEST(u1.Add_Entry()==5);
        TEST(u1.Add_Entry()==6);
        u1.Union(4,5);

        return ok?success:failure;
    }

    int Number_Of_Tests() const
    {return 1;}

//#####################################################################
};
static UNION_FIND_TESTS union_find_tests;
}
