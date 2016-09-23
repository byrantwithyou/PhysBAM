//#####################################################################
// Copyright 2008, Justin Solomon.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HASHTABLE_TESTS
//#####################################################################

#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Utilities/TEST_BASE.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

class HASHTABLE_TESTS:public TEST_BASE
{
    bool Test(const bool test,const char* str_test,int line,const char* str_T,const char* str_K,bool& ok) const
    {if(test) return true;LOG::cout<<"FAILED: "<<str_test<<"    T "<<str_T<<"    K "<<str_K<<"    line "<<line<<std::endl;ok=false;return false;}

public:
    HASHTABLE_TESTS()
        :TEST_BASE("hashtable")
    {}

    virtual ~HASHTABLE_TESTS(){}

#define TEST(x) Test((x),#x,__LINE__,typeid(T).name(),typeid(K).name(),ok)

    template<class T,class K>
    void Test_Type(T a,T b,T c,T d,K r,K s,K t,K u,bool& ok)
    {
        HASHTABLE<T> h;
        HASHTABLE<T,K> k;
        ARRAY<T> keys;
        ARRAY<K> data;

        TEST(h.Size()==0);
        TEST(k.Size()==0);

        h.Initialize_New_Table(50);
        TEST(h.Size()==0);
        TEST(!h.Contains(a));

        h.Insert(a);
        TEST(h.Contains(a));
        TEST(!h.Contains(b));
        TEST(h.Size()==1);

        h.Insert(b);
        TEST(h.Size()==2);
        TEST(h.Contains(a));
        TEST(h.Contains(b));
        TEST(!h.Contains(c));

        TEST(!h.Set(b));
        TEST(h.Contains(b));
        TEST(h.Size()==2);

        TEST(h.Set(c));
        TEST(h.Size()==3);
        TEST(h.Contains(c));

        TEST(!h.Contains(d));
        TEST(!h.Delete_If_Present(d));
        TEST(!h.Contains(d));
        TEST(h.Size()==3);

        TEST(h.Contains(a));
        TEST(h.Delete_If_Present(a));
        TEST(!h.Contains(a));
        TEST(h.Size()==2);

        h.Clean_Memory();
        TEST(h.Size()==0);

        VECTOR<T,3> V(a,b,c);
        h.Set_All(V);
        TEST(h.Size()==3);
        TEST(h.Contains(a));
        TEST(h.Contains(b));
        TEST(h.Contains(c));
        TEST(!h.Contains(d));

        h.Get_Keys(keys);
        TEST(keys.m==3);

        k.Insert(a,r);
        TEST(k.Contains(a));
        TEST(!k.Contains(b));
        TEST(k.Size()==1);

        k.Insert(b,s);
        TEST(k.Size()==2);
        TEST(k.Contains(a));
        TEST(k.Contains(b));
        TEST(!k.Contains(c));

        TEST(!k.Set(b,t));
        TEST(k.Contains(b));
        TEST(k.Size()==2);

        TEST(k.Set(c,t));
        TEST(k.Size()==3);
        TEST(k.Contains(c));

        TEST(!k.Contains(d));
        TEST(!k.Delete_If_Present(d));
        TEST(!k.Contains(d));
        TEST(k.Size()==3);

        TEST(k.Contains(a));
        TEST(k.Delete_If_Present(a));
        TEST(!k.Contains(a));
        TEST(k.Size()==2);

        k.Clean_Memory();
        TEST(k.Size()==0);

        TEST(k.Get_Or_Insert(a,r)==r);
        TEST(k.Contains(a));
        TEST(k.Size()==1);
        TEST(k.Get_Or_Insert(a,s)==r);
        TEST(k.Size()==1);
        TEST(k.Get_Or_Insert(a)==r);
        TEST(k.Size()==1);
        k.Get_Or_Insert(a)=s;
        TEST(k.Get_Or_Insert(a)==s);

        TEST(k.Get_Or_Insert(b)==K());
        TEST(k.Contains(b));
        TEST(k.Size()==2);

        K* x = 0;
        TEST((x = k.Get_Pointer(a)));
        TEST(*x==s);
        TEST(k.Size()==2);
        TEST(!(x = k.Get_Pointer(c)));
        TEST(k.Size()==2);

        TEST(k.Get(a)==s);
        try
        {
            TEST(k.Get(c)==s);
            TEST(!"should not get here");
        }
        catch(...)
        {
        }

        TEST(k.Get_Default(a)==s);
        TEST(k.Size()==2);
        TEST(k.Get_Default(d)==K());
        TEST(k.Size()==2);

        k.Exchange(a,b);
        TEST(k.Size()==2);
        TEST(k.Get_Or_Insert(a)==K());
        TEST(k.Get_Or_Insert(b)==s);

        keys.Remove_All();
        k.Get_Keys(keys);
        TEST(keys.m==2);

        k.Get_Data(data);
        TEST(data.m==2);
    }

    TEST_RESULT Run_Test(int n)
    {
        bool ok=true;
        HASHTABLE<int> hi;

        Test_Type(1,3,-1,10,1.0,1e7,-1.0,0.0,ok);
        Test_Type(VECTOR<int,1>(1),VECTOR<int,1>(3),VECTOR<int,1>(-1),VECTOR<int,1>(10),0,1,2,3,ok);
        Test_Type(VECTOR<int,2>(1,0),VECTOR<int,2>(0,1),VECTOR<int,2>(0,0),VECTOR<int,2>(1,1),0,1,2,3,ok);
        Test_Type(VECTOR<int,3>(1,3,0),VECTOR<int,3>(0,3,1),VECTOR<int,3>(0,0,1),VECTOR<int,3>(1,3,1),0,1,2,3,ok);
        Test_Type(VECTOR<int,4>(1,3,0,0),VECTOR<int,4>(1,0,3,1),VECTOR<int,4>(0,0,0,1),VECTOR<int,4>(0,1,3,1),0,1,2,3,ok);
        Test_Type("a","","b","aa",0,1,2,3,ok);

        return ok?success:failure;
    }

    int Number_Of_Tests() const
    {return 1;}

//#####################################################################
};
static HASHTABLE_TESTS hashtable_tests;
}
