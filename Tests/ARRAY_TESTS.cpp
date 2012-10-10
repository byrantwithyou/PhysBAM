//#####################################################################
// Copyright 2012, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_TESTS
//#####################################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_BASE.h>
#include <PhysBAM_Tools/Utilities/TEST_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

class ARRAY_TESTS:public TEST_BASE
{
    bool Test(const bool test,const char* str_test,int line,bool& ok) const
    {if(test) return true;LOG::cout<<"FAILED: "<<str_test<<"    line "<<line<<std::endl;ok=false;return false;}

public:
    ARRAY_TESTS()
        :TEST_BASE("array")
    {}

    virtual ~ARRAY_TESTS(){}

#define TEST(x) Test((x),#x,__LINE__,ok)

    TEST_RESULT Run_Test(int n)
    {
        bool ok=true;
        ARRAY<int> ar1;
        ARRAY<int> ar2(4);
        
        TEST(ar1.Size()==0);
        TEST(ar2.Size()==4);
        
        ar2(0) = 7;
        ar2(1) = 5;
        ar2(2) = 11;
        ar2(3) = 6;

        ar2.Resize(5);
        TEST(ar2.Size()==5);
        TEST(ar2(4)==0);

        ar2.Exact_Resize(4);
        TEST(ar2.Size()==4);
        
        ar1.Append(7);
        ar2.Append(7);
        TEST(ar1.Size()==1);
        TEST(ar2.Size()==5);
        TEST(ar2(4)==7);
        TEST(ar2(0)==7);
        TEST(ar2(1)==5);
        
        ar2.Append_Unique(7);
        TEST(ar2.Size()==5);
        ar2.Append_Unique(8);
        TEST(ar2(5)==8);
        TEST(ar2.Size()==6);

        ar2.Remove_End();
        TEST(ar2.Size()==5);
        ar2.Remove_End();
        TEST(ar2.Size()==4);
        TEST(ar2(3)==6);

        ar2.Remove_Index(0);
        ar2.Remove_Index(2);
        TEST(ar2.Size()==2);
        TEST(ar2(0)==5);
        TEST(ar2(1)==11);

        ar2.Remove_All();
        TEST(ar2.Size()==0);
        ar2.Append(7);
        ar2.Append(5);
        ar2.Append(11);
        ar2.Append(6);       
        TEST(ar2.Size()==4);

        ar2.Insert(7,0);
        TEST(ar2(0)==7);
        TEST(ar2(1)==7);
        ar2.Insert(6,5);
        TEST(ar2(5)==6);
        TEST(ar2.Size()==6);

        TEST(ar2.Pop()==6);
        TEST(ar2.Size()==5);

        ar1.Clean_Memory();
        TEST(ar1.Size()==0);
        ar1.Exchange(ar2);
        TEST(ar2.Size()==0);
        TEST(ar1.Size()==5);

        ar1.Remove_All();
        ar1.Append(7);
        ar1.Append(5);
        ar1.Append(11);
        ar1.Append(6);    
        
        TEST(ar1.Max()==11);
        TEST(ar1.Min()==5);
        TEST(ar1(ar1.Arg_Min())==ar1.Min());
        TEST(ar1(ar1.Arg_Max())==ar1.Max());

        ar1(0)=20;
        ar1(3)=1;
        TEST(ar1.Max()==20);
        TEST(ar1.Min()==1);
        TEST(ar1(ar1.Arg_Min())==ar1.Min());
        TEST(ar1(ar1.Arg_Max())==ar1.Max());

        ar1(0)=0;
        ar1(3)=100;
        TEST(ar1.Max()==100);
        TEST(ar1.Min()==0);
        TEST(ar1(ar1.Arg_Min())==ar1.Min());
        TEST(ar1(ar1.Arg_Max())==ar1.Max());

        TEST(ar1.Find(ar1(0))==0);
        TEST(ar1.Find(ar1(3))==3);
        TEST(ar1.Find(30)==-1);
        TEST(ar1.Contains(ar1(0))==true);
        TEST(ar1.Contains(ar1(1))==true);
        TEST(ar1.Contains(ar1(3))==true);
        TEST(ar1.Find(ar1.Last())==ar1.Size()-1);

        ar2=ar1;
        TEST(ar1==ar2);
        ar2(0)++;
        TEST(ar1!=ar2);

        ar2=VECTOR<int,2>(1,3);
        TEST(ar2.Size()==2);
        TEST(ar2(0)==1);
        TEST(ar2(1)==3);

        TEST(!ar1.Valid_Index(-1));
        TEST(!ar1.Valid_Index(ar1.m));
        TEST(ar1.Valid_Index(ar1.m-1));
        TEST(ar1.Valid_Index(0));
        TEST(*ar1.Get_Array_Pointer()==ar1(0));

        ar2=ar1;
        ar1.Resize(5);
        ar2.Append(0);
        TEST(ar1==ar2);

        ar2=ar1=ARRAY<int>(VECTOR<int,3>(5,4,3));
        ar1.Append(7);
        ar1.Append(8);
        ar1.Append(9);
        ar1.Append(10);
        ar2.Append_Elements(VECTOR<int,4>(7,8,9,10));
        TEST(ar1==ar2);
        ar1.Append(11);
        ar1.Append(12);
        ar2.Append_Unique_Elements(VECTOR<int,5>(11,7,8,9,12));
        TEST(ar1==ar2);

        ar1.Remove_Index_Lazy(3);
        TEST(ar1.m==ar2.m-1);
        ar2.Remove_Index_Lazy(3);
        TEST(ar1.m==ar2.m);
        TEST(ar1==ar2);

        return ok?success:failure;
    }

    int Number_Of_Tests() const
    {return 1;}

//#####################################################################
};
static ARRAY_TESTS array_tests;
}
