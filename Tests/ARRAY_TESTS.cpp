//#####################################################################
// Copyright 2012, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_TESTS
//#####################################################################

#include <PhysBAM_Tools/Arrays/ARRAY_BASE.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/TEST_BASE.h>
namespace PhysBAM{

class ARRAY_TESTS:public TEST_BASE
{
    bool Test(const bool test,const char* str_test,int line,const char* str_T,bool& ok) const
    {if(test) return true;LOG::cout<<"FAILED: "<<str_test<<"    T "<<str_T<<"    line "<<line<<std::endl;ok=false;return false;}

public:
    ARRAY_TESTS()
        :TEST_BASE("array")
    {}

    virtual ~ARRAY_TESTS(){}

#define TEST(x) Test((x),#x,__LINE__,typeid(T).name(),ok)

    template<class T>
    void Test_Type(T a,T b,T c,T d,T M,T m,int argM,int argm,bool& ok)
    {
        ARRAY<T> ar1;
        ARRAY<T> ar2(4);
        
        TEST(ar1.Size()==0);
        TEST(ar2.Size()==4);
        
        ar2(0) = a;
        ar2(1) = b;
        ar2(2) = c;
        ar2(3) = d;

        ar2.Resize(5);
        TEST(ar2.Size()==5);

        ar2.Exact_Resize(4);
        TEST(ar2.Size()==4);
        
        ar1.Append(a);
        ar2.Append(a);
        TEST(ar1.Size()==1);
        TEST(ar2.Size()==5);
        TEST(ar2(4)==a);
        TEST(ar2(0)==a);
        TEST(ar2(1)==b);
        
        ar2.Append_Unique(a);
        TEST(ar2.Size()==5);

        ar2.Remove_End();
        TEST(ar2.Size()==4);
        TEST(ar2(3)==d);

        ar2.Remove_Index(0);
        ar2.Remove_Index(2);
        TEST(ar2.Size()==2);
        TEST(ar2(0)==b);
        TEST(ar2(1)==c);

        ar2.Remove_All();
        TEST(ar2.Size()==0);
        ar2.Append(a);
        ar2.Append(b);
        ar2.Append(c);
        ar2.Append(d);       

        ar2.Insert(a,0);
        TEST(ar2(0)==a);
        TEST(ar2(1)==a);
        ar2.Insert(d,5);
        TEST(ar2(5)==d);

        TEST(ar2.Pop()==d);
        TEST(ar2.Size()==5);

        ar1.Clean_Memory();
        TEST(ar1.Size()==0);
        ar1.Exchange(ar2);
        TEST(ar2.Size()==0);
        TEST(ar1.Size()==5);

        ar1.Remove_All();
        ar1.Append(a);
        ar1.Append(b);
        ar1.Append(c);
        ar1.Append(d);    
        
        TEST(ar1.Max()==M);
        TEST(ar1.Min()==m);
        TEST(ar1.Argmin()==argm);
        TEST(ar1.Argmax()==argM);

        TEST(ar1.Find(m)==argm);
        TEST(ar1.Find(M)==argM);
        TEST(ar1.Contains(a)==true);
        TEST(ar1.Contains(b)==true);
        TEST(ar1.Contains(c)==true);
        TEST(ar1.Contains(d)==true);        
        TEST(ar1.Contains(a*100)==false);       
    }

    TEST_RESULT Run_Test(int n)
    {
        bool ok=true;
        Test_Type(1,2,3,4,4,1,3,0,ok);
        Test_Type(4,3,2,1,4,1,0,3,ok);
        Test_Type(1,4,2,3,4,1,1,0,ok);
        Test_Type(2,1,4,3,4,1,2,1,ok);
        Test_Type(2,1,2,4,4,1,3,1,ok);

        return ok?success:failure;
    }

    int Number_Of_Tests() const
    {return 1;}

//#####################################################################
};
static ARRAY_TESTS array_tests;
}
