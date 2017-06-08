//#####################################################################
// Copyright 2017, Ounan Ding, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_RANGE_ITERATOR_TESTS
//#####################################################################

#include <Core/Arrays/SORT.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE_ITERATOR.h>
#include <Core/Utilities/TEST_BASE.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/FACE_RANGE_ITERATOR.h>
namespace PhysBAM{

static bool fr_test_ok=false;
template<int d>
bool operator<(FACE_INDEX<d> a,FACE_INDEX<d> b)
{
    for(int i=0;i<d;i++)
        if(a.index(i)!=b.index(i))
            return a.index(i)<b.index(i);
    return a.axis<b.axis;
}
template<int d>
void All_Flags_Tests(RF flags)
{
    typedef VECTOR<int,d> TV_INT;
    RANGE<TV_INT> ri,ro;
    for(int i=0;i<d;i++) ri.min_corner(i)=3-i;
    for(int i=0;i<d;i++) ri.max_corner(i)=6;
    for(int i=0;i<d;i++) ro.min_corner(i)=-1-i;
    for(int i=0;i<d;i++) ro.max_corner(i)=7+2*i;

    for(int s=-1;s<2*d;s++)
    {
        ARRAY<FACE_INDEX<d> > a,b;
        for(FACE_RANGE_ITERATOR<d> it(ro,ri,flags,-1,s);it.Valid();it.Next())
            a.Append(it.face);
        for(int i=0;i<d;i++){
            bool failed_axis_test=false;
            for(FACE_RANGE_ITERATOR<d> it(ro,ri,flags,i,s);it.Valid();it.Next()){
                if(it.face.axis!=i && !failed_axis_test){
                    failed_axis_test=true;
                    fr_test_ok=false;
                    LOG::printf("axis iteration exclusion test failed; side=%i axis=%i flags=%i\n",s,i,(int)flags);}
                b.Append(it.face);}}
        if(a!=b){
            LOG::printf("axis iteration coverage test failed; side=%i flags=%i\n",s,(int)flags);
            fr_test_ok=false;}
    }
    for(int i=-1;i<d;i++)
        for(int s=-1;s<2*d;s++){
            ARRAY<FACE_INDEX<d> > a;
            for(FACE_RANGE_ITERATOR<d> it(ro,ri,flags,i,s);it.Valid();it.Next())
                a.Append(it.face);
            for(FACE_RANGE_ITERATOR<d> it(ro,ri,flags|RF::reverse,i,s);it.Prev_Valid();it.Prev()){
                if(a.Pop()!=it.face){
                    fr_test_ok=false;
                    LOG::printf("reverse test failed; side=%i axis=%i flags=%i\n",s,i,(int)flags);
                    break;}}}
}

template<int d>
void All_Side_Tests(RF flags,bool ignore_dup)
{
    typedef VECTOR<int,d> TV_INT;
    RANGE<TV_INT> ri,ro;
    for(int i=0;i<d;i++) ri.min_corner(i)=3-i;
    for(int i=0;i<d;i++) ri.max_corner(i)=6;
    for(int i=0;i<d;i++) ro.min_corner(i)=-1-i;
    for(int i=0;i<d;i++) ro.max_corner(i)=7+2*i;
    RF io[4]={RF::none,RF::skip_inner,RF::skip_outer,RF::skip_inner|RF::skip_outer};

    for(int i=0;i<4;i++){
        ARRAY<FACE_INDEX<d> > a,b,c;
        for(FACE_RANGE_ITERATOR<d> it(ro,ri,flags|io[i],-1,-1);it.Valid();it.Next())
            a.Append(it.face);
        for(int s=0;s<2*d;s++)
            for(FACE_RANGE_ITERATOR<d> it(ro,ri,flags|io[i],-1,s);it.Valid();it.Next())
                b.Append(it.face);
        if(ignore_dup){
            Get_Unique(c,a);
            a=c;
            Get_Unique(c,b);
            b=c;}
        a.Sort();
        b.Sort();
        if(a!=b){
            LOG::printf("side test failed; flags=%i\n",(int)(flags|io[i]));
            fr_test_ok=false;}}
}

template<int d>
void Delay_Tests(RF flags)
{
    typedef VECTOR<int,d> TV_INT;
    RANGE<TV_INT> ri,ro;
    for(int i=0;i<d;i++) ri.min_corner(i)=3-i;
    for(int i=0;i<d;i++) ri.max_corner(i)=6;
    for(int i=0;i<d;i++) ro.min_corner(i)=-1-i;
    for(int i=0;i<d;i++) ro.max_corner(i)=7+2*i;
    RF io[4]={RF::none,RF::skip_inner,RF::skip_outer,RF::skip_inner|RF::skip_outer};

    for(int i=0;i<4;i++){
        ARRAY<FACE_INDEX<d> > a,b;
        for(FACE_RANGE_ITERATOR<d> it(ro,ri,flags|io[i],-1,-1);it.Valid();it.Next())
            a.Append(it.face);
        for(FACE_RANGE_ITERATOR<d> it(ro,ri,flags|io[i]|RF::delay_corners,-1,-1);it.Valid();it.Next())
            b.Append(it.face);
        a.Sort();
        b.Sort();
        if(a!=b){
            LOG::printf("delay test failed; flags=%i\n",(int)(flags|io[i]));
            fr_test_ok=false;}}
}

template<int d>
void Single_Side_Tests()
{
    typedef VECTOR<int,d> TV_INT;
    RANGE<TV_INT> ri,ro;
    for(int i=0;i<d;i++) ri.min_corner(i)=3-i;
    for(int i=0;i<d;i++) ri.max_corner(i)=6;
    for(int i=0;i<d;i++) ro.min_corner(i)=-1-i;
    for(int i=0;i<d;i++) ro.max_corner(i)=7+2*i;
    RF io[4]={RF::none,RF::skip_inner,RF::skip_outer,RF::skip_inner|RF::skip_outer};

    for(int i=0;i<4;i++){
        for(int s=0;s<2*d;s++){
            ARRAY<FACE_INDEX<d> > a,b,c;
            for(FACE_RANGE_ITERATOR<d> it(ro,ri,io[i],-1,s);it.Valid();it.Next())
                a.Append(it.face);
            for(FACE_RANGE_ITERATOR<d> it(ro,ri,io[i]|RF::delay_corners,-1,s);it.Valid();it.Next())
                b.Append(it.face);
            for(FACE_RANGE_ITERATOR<d> it(ro,ri,io[i]|RF::duplicate_corners,-1,s);it.Valid();it.Next())
                c.Append(it.face);
            a.Sort();
            b.Sort();
            c.Sort();
            if(a!=b || a!=c){
                LOG::printf("single side test failed; flags=%i side=%i (sizes %i %i %i)\n",(int)(io[i]),s,a.m,b.m,c.m);
                fr_test_ok=false;}}}
}

template<int d>
void Containment_Tests(RF fi,RF fo,RF fx)
{
    typedef VECTOR<int,d> TV_INT;
    RANGE<TV_INT> ri,ro,rm;
    for(int i=0;i<d;i++) ri.min_corner(i)=3-i;
    for(int i=0;i<d;i++) ri.max_corner(i)=6;
    for(int i=0;i<d;i++) rm.min_corner(i)=-1-i;
    for(int i=0;i<d;i++) rm.max_corner(i)=7+2*i;
    for(int i=0;i<d;i++) ro.min_corner(i)=-4-i;
    for(int i=0;i<d;i++) ro.max_corner(i)=10+3*i;

    ARRAY<FACE_INDEX<d> > a,b;
    for(FACE_RANGE_ITERATOR<d> it(ro,ri,fx);it.Valid();it.Next())
        a.Append(it.face);
    for(FACE_RANGE_ITERATOR<d> it(rm,ri,fi);it.Valid();it.Next())
        b.Append(it.face);
    for(FACE_RANGE_ITERATOR<d> it(ro,rm,fo);it.Valid();it.Next())
        b.Append(it.face);
    a.Sort();
    b.Sort();
    if(a!=b){
        LOG::printf("containment test failed; flags: in=%i out=%i total=%i\n",(int)fi,(int)fo,(int)fx);
        fr_test_ok=false;}
}
                               
template<int d>
void Dim_Tests()
{
    RF io[4]={RF::none,RF::skip_inner,RF::skip_outer,RF::skip_inner|RF::skip_outer};
    RF epu[5]={RF::none,RF::duplicate_corners,RF::delay_corners,
               RF::partial_single_side,RF::partial_single_side|RF::delay_corners};
    for(int i=0;i<4;i++)
        for(int j=0;j<5;j++)
            All_Flags_Tests<d>(io[i]|epu[j]);
    All_Flags_Tests<d>(RF::interior);
    All_Flags_Tests<d>(RF::interior|RF::skip_outer);

    All_Side_Tests<d>(RF::none,true);
    All_Side_Tests<d>(RF::duplicate_corners,false);
    All_Side_Tests<d>(RF::delay_corners,true);
    All_Side_Tests<d>(RF::partial_single_side,false);
    All_Side_Tests<d>(RF::partial_single_side|RF::delay_corners,false);

    Delay_Tests<d>(RF::none);
    Delay_Tests<d>(RF::partial_single_side);

    Single_Side_Tests<d>();

    RF op_i[2]={RF::none,RF::skip_inner};
    RF op_o[2]={RF::none,RF::skip_outer};

    for(int i=0;i<2;i++)
        Containment_Tests<d>(RF::interior|RF::skip_outer,op_o[i],RF::interior|op_o[i]);
    for(int i=0;i<2;i++)
        Containment_Tests<d>(RF::interior,RF::skip_inner|op_o[i],RF::interior|op_o[i]);
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
            Containment_Tests<d>(RF::skip_outer|op_i[j],op_o[i],op_o[i]|op_i[j]);
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
            Containment_Tests<d>(op_i[j],RF::skip_inner|op_o[i],op_o[i]|op_i[j]);
}

class FACE_RANGE_ITERATOR_TESTS:public TEST_BASE
{
public:
    FACE_RANGE_ITERATOR_TESTS()
        :TEST_BASE("face_range")
    {puts("blah");}

    virtual ~FACE_RANGE_ITERATOR_TESTS(){}

    TEST_RESULT Run_Test(int n)
    {
        fr_test_ok=true;

        Dim_Tests<1>();
        Dim_Tests<2>();
        Dim_Tests<3>();
        
        return fr_test_ok?success:failure;
    }

    int Number_Of_Tests() const
    {return 1;}

//#####################################################################
};
static FACE_RANGE_ITERATOR_TESTS fr_tests;
}
