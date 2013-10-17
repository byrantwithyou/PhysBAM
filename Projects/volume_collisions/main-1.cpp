//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Auto_Diff/AUTO_HESS.h>
#include <Tools/Auto_Diff/AUTO_HESS_EXT.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Data_Structures/UNION_FIND.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/permutation.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <sstream>

using namespace PhysBAM;
using namespace PhysBAM::HETERO_DIFF;

typedef float RW;
typedef VECTOR<int,3> TV_INT;

const int max_pts=14;
typedef VECTOR<int,max_pts> PL;

const int in_flag=0x100;

void pr(const PL& p,const char* s="")
{
    for(int j=0;j<12;j++)
    {
        if(!p(j)) continue;
        for(int i=0;i<8;i++)
            if(p(j)&(1<<i))
                putchar("abcdrstu"[i]);
        if(j!=11) putchar(' ');
    }
    printf("%s",s);
}

struct HISTORY
{
    PL poly;
    int num_tests,num_ded_tests,num_safe;
    int tests[30];
    int ded_tests[14];
    int safe_pts[20];

    int p,grow_in;
    PL in,out;
    int first_poly,first_empty_in,first_empty_out;
};

ARRAY<HISTORY> final_history;

HASHTABLE<PL,int> poly_lookup;

VECTOR<int,256> maps[24*24];

void Evolve_History(HISTORY& h)
{
    PHYSBAM_ASSERT(h.num_tests<40);
    PHYSBAM_ASSERT(h.num_ded_tests<40);

    if(h.grow_in<0){
        int pt=h.poly(h.first_poly),t=(1<<h.p)|pt;
        h.poly(h.first_poly++)=0;

        HISTORY h2(h);
        h2.tests[h2.num_tests++]=in_flag|t;
        h2.in[h2.first_empty_in++]=pt;
        h2.grow_in=1;
        Evolve_History(h2);

        h.tests[h.num_tests++]=t;
        h.out[h.first_empty_out++]=pt;
        h.grow_in=0;}

    if(h.first_poly<max_pts){
        PL& pl=h.grow_in?h.in:h.out;
        int n_pl=h.grow_in?h.first_empty_in:h.first_empty_out;

        for(int i=0;i<n_pl;i++)
            for(int j=h.first_poly;j<max_pts;j++){
                int is=pl(i)&h.poly(j);
                if(is&(is-1)){
                    exchange(h.poly(j),h.poly(h.first_poly));

                    int pt=h.poly(h.first_poly),t=(1<<h.p)|pt;
                    h.poly(h.first_poly++)=0;

                    HISTORY h2(h);
                    h2.tests[h2.num_tests++]=in_flag|t;
                    h2.in[h2.first_empty_in++]=pt;
                    Evolve_History(h2);

                    h.tests[h.num_tests++]=t;
                    h.out[h.first_empty_out++]=pt;
                    Evolve_History(h);
                    return;}}}

    PL& pl=h.grow_in?h.out:h.in;
    int& n_pl=h.grow_in?h.first_empty_out:h.first_empty_in;
    int rest_flag=h.grow_in?0:in_flag;
    for(int j=h.first_poly;j<max_pts;j++){
        pl(n_pl++)=h.poly(j);
        h.ded_tests[h.num_ded_tests++]=rest_flag|(1<<h.p)|h.poly(j);}

    int num_pl_1=1;
    PL pl_0=pl,pl_1;
    pl_1(0)=pl_0(0);
    pl_0(0)=0;
    for(int i=0;i<max_pts;i++)
        for(int j=0;j<max_pts;j++){
            int is=pl_1(i)&pl_0(j);
            if(is&(is-1)){
                pl_1(num_pl_1++)=pl_0(j);
                pl_0(j)=0;}}
    if(num_pl_1!=n_pl){
        h.poly.Fill(-1);
        poly_lookup.Set(h.poly,final_history.Append(h));
        return;}

    h.poly=h.in;
    int first_empty_poly=h.first_empty_in;
    for(int i=0;i<h.first_empty_in;i++)
        for(int j=0;j<h.first_empty_out;j++){
            int is=h.in(i)&h.out(j);
            if(is&(is-1))
                h.poly[first_empty_poly++]=is|(1<<h.p);}

    for(int i=0;i<max_pts;i++)
        for(int j=i+1;j<max_pts;j++)
            for(int k=j+1;k<max_pts;k++){
                int is=h.poly(i)&h.poly(j)&h.poly(k);
                PHYSBAM_ASSERT((is&(is-1))==0);}

    h.poly.Sort();
    h.first_poly=max_pts-first_empty_poly;
    for(int i=0;i<h.first_empty_out;i++) h.safe_pts[h.num_safe++]=h.out[i];

    h.in.Fill(0);
    h.out.Fill(0);
    h.first_empty_in=0;
    h.first_empty_out=0;
    h.p++;
    h.grow_in=-1;
    if(h.p<8 && h.first_poly<max_pts)
        return Evolve_History(h);

    poly_lookup.Set(h.poly,final_history.Append(h));
    if(final_history.m%10000==0) printf("%i\n",final_history.m);
}

int main(int argc, char* argv[])
{
    for(int i=0;i<24;i++)
        for(int j=0;j<24;j++)
        {
            VECTOR<int,4> a=permute_four(VECTOR<int,4>(0,1,2,3),i);
            VECTOR<int,4> b=permute_four(VECTOR<int,4>(4,5,6,7),j);
            VECTOR<int,8> M=a.Append_Elements(b);
            for(int k=1;k<256;k++)
            {
                int r=0;
                for(int l=0;l<8;l++)
                    if(k&(1<<l))
                        r|=1<<M(l);
                maps[i*24+j](k)=r;
            }
        }

    ARRAY<HISTORY> hist_stack;
    PL init_volume;
    init_volume(max_pts-1)=14;
    init_volume(max_pts-2)=13;
    init_volume(max_pts-3)=11;
    init_volume(max_pts-4)=7;

    HISTORY init_hist = {init_volume,0,0,0,{},{},{},4,-1,PL(),PL(),max_pts-4,0,0};
    Evolve_History(init_hist);



    printf("%i\n", final_history.m);

    int num_tests=0,num_ded_tests=0,num_safe=0,num_bad=0;

    for(int i=0;i<final_history.m;i++){
        num_tests=std::max(num_tests,final_history(i).num_tests);
        num_ded_tests=std::max(num_ded_tests,final_history(i).num_ded_tests);
        num_safe=std::max(num_safe,final_history(i).num_safe);
        if(final_history(i).poly.Min()<0) num_bad++;}

    printf("%i %i %i %i\n", num_tests,num_ded_tests,num_safe,num_bad);

    return 0;
}

