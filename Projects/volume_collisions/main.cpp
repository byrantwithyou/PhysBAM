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

#ifndef COMPILE_ID_TYPES_AS_INT
PHYSBAM_DECLARE_ELEMENT_ID(POLY_LIST_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(HISTORY_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(PRIMITIVE_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
PHYSBAM_DECLARE_ELEMENT_ID(REDUCTION_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
#else
typedef int POLY_LIST_ID;
typedef int HISTORY_ID;
typedef int PRIMITIVE_ID;
typedef int REDUCTION_ID;
#endif

typedef float RW;
typedef VECTOR<int,3> TV_INT;

const int max_pts=14;
typedef VECTOR<int,max_pts> PL;

const int in_flag=0x100;

template<int n>
void pr(VECTOR<int,n> p,const char* s="")
{
    p.Sort();
    for(int j=0;j<n;j++)
    {
        if(!p(j)) continue;
        for(int i=0;i<8;i++)
            if(p(j)&(1<<i))
                putchar("abcdrstu"[i]);
        if(j!=max_pts-1) putchar(' ');
    }
    printf("%s",s);
}

struct HISTORY
{
    PL poly;
    int num_tests,num_ded_tests;
    VECTOR<int,30> tests;
    VECTOR<int,14> ded_tests;
    VECTOR<bool,256> safe;

    int p,grow_in;
    PL in,out;
    int first_poly,first_empty_in,first_empty_out;
    REDUCTION_ID reduction;
};

void pr(const HISTORY& h,const char* s="")
{
    pr(h.poly," + ");
    pr(h.in," + ");
    pr(h.out," + ");
    for(int i=0;i<h.num_tests;i++){
        int n=h.tests(i);
        putchar(n>=0?'+':'-');
        for(int j=0;j<8;j++)
            if(n&(1<<j)) putchar("abcdrstu"[j]);
        putchar(' ');}
    printf(" # ");
    for(int i=0;i<h.num_ded_tests;i++){
        int n=h.ded_tests(i);
        putchar(n>=0?'+':'-');
        for(int j=0;j<8;j++)
            if(n&(1<<j)) putchar("abcdrstu"[j]);
        putchar(' ');}
    printf("%s",s);
}

ARRAY<PL,POLY_LIST_ID> poly_list;
HASHTABLE<PL,POLY_LIST_ID> poly_lookup;
ARRAY<HISTORY,HISTORY_ID> final_history;

bool Evolve_History(HISTORY& h)
{
    PHYSBAM_ASSERT(h.num_tests<40);
    PHYSBAM_ASSERT(h.num_ded_tests<40);

    if(h.grow_in<0){
        int pt=h.poly(h.first_poly),t=(1<<h.p)|pt;
        h.poly(h.first_poly++)=0;

        HISTORY h2(h);

        for(int j=0;j<h.safe.m;j++){PHYSBAM_ASSERT((*(char*)&h.safe(j)&~1)==0);}
        for(int j=0;j<h2.safe.m;j++){PHYSBAM_ASSERT((*(char*)&h2.safe(j)&~1)==0);}

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

                    for(int j=0;j<h.safe.m;j++){PHYSBAM_ASSERT((*(char*)&h.safe(j)&~1)==0);}
                    for(int j=0;j<h2.safe.m;j++){PHYSBAM_ASSERT((*(char*)&h2.safe(j)&~1)==0);}

                    HISTORY_ID num_final=final_history.m;
                    int h2_num_tests=h2.num_tests;
                    h2.tests[h2.num_tests++]=in_flag|t;
                    h2.in[h2.first_empty_in++]=pt;
                    bool valid=Evolve_History(h2);

                    if(valid) h.tests[h.num_tests++]=t;
                    else h.ded_tests[h.num_ded_tests++]=t;
                    h.out[h.first_empty_out++]=pt;
                    if(Evolve_History(h)) return true;
                    if(!valid) return false;

                    for(HISTORY_ID i=num_final;i<final_history.m;i++){
                        HISTORY& hh=final_history(i);
                        hh.ded_tests[hh.num_ded_tests++]=in_flag|t;
                        for(int j=h2_num_tests+1;j<hh.num_tests;j++)
                            hh.tests[j]=hh.tests[j+1];
                        hh.num_tests--;}

                    return true;}}}

    PL& pl=h.grow_in?h.out:h.in;
    int& n_pl=h.grow_in?h.first_empty_out:h.first_empty_in;
    int rest_flag=h.grow_in?0:in_flag;
    for(int j=h.first_poly;j<max_pts;j++){
        pl(n_pl++)=h.poly(j);
        h.ded_tests[h.num_ded_tests++]=rest_flag|(1<<h.p)|h.poly(j);}

    if(n_pl>1){
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
        if(num_pl_1!=n_pl) return false;}

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
    for(int i=h.first_poly;i<max_pts;i++) h.safe(h.poly[i])=true;

    h.in.Fill(0);
    h.out.Fill(0);
    h.first_empty_in=0;
    h.first_empty_out=0;
    h.p++;
    h.grow_in=-1;
    if(h.p<8 && h.first_poly<max_pts)
        return Evolve_History(h);

    for(int i=0;i<8;i++)
        for(int j=i+1;j<8;j++)
            for(int k=j+1;k<8;k++){
                int ij=(1<<i)|(1<<j),ijk=ij|(1<<k);
                if(h.safe(ijk)) continue;
                int has[2]={0,0};
                for(int l=0;l<h.num_tests;l++)
                    if((h.tests[l]&ijk)==ijk)
                        has[h.tests[l]>>8]++;
                if(has[0] && has[1]) h.safe(ijk)=true;}

    final_history.Append(h);
    if(poly_lookup.Get_Or_Insert(h.poly,poly_list.m)==poly_list.m)
        poly_list.Append(h.poly);
    if(Value(final_history.m)%100000==0) LOG::printf("%P\n",final_history.m);
    return true;
}

void Deduce_Safe(HISTORY& h)
{
    // 0=unknown, 1=safe, 2=unsafe
    int state[256]={0};
    for(int i=0;i<256;i++) if(h.safe(i)) state[i]=1;
    bool progress=true;
    while(progress){
        progress=false;
        for(int a=0;a<4;a++)
            for(int u=4;u<8;u++)
                for(int b=0;b<4;b++)
                    if(b!=a){
                        int abu=(1<<a)|(1<<b)|(1<<u);
                        if(int s0=state[abu]){
                            for(int c=b+1;c<4;c++)
                                if(c!=a){
                                    int acu=(1<<a)|(1<<c)|(1<<u);
                                    if(int s1=state[acu]){
                                        int other=(15^abu^acu)|(1<<u);
                                        if(state[other]) continue;
                                        progress=true;
                                        if(s0==s1) state[other]=2;
                                        else state[other]=1;}}}}
        for(int a=4;a<8;a++)
            for(int u=0;u<4;u++)
                for(int b=4;b<8;b++)
                    if(b!=a){
                        int abu=(1<<a)|(1<<b)|(1<<u);
                        if(int s0=state[abu]){
                            for(int c=b+1;c<8;c++)
                                if(c!=a){
                                    int acu=(1<<a)|(1<<c)|(1<<u);
                                    if(int s1=state[acu]){
                                        int other=(15*16^abu^acu)|(1<<u);
                                        if(state[other]) continue;
                                        progress=true;
                                        if(s0==s1) state[other]=2;
                                        else state[other]=1;}}}}}
    for(int i=0;i<256;i++) if(state[i]==1) h.safe(i)=true;
}

// type 0: c = a, type 1: c = a + b, type 2: c = a - b, type 3: replace with rule a
struct REDUCTION_TEMPLATE
{
    POLY_LIST_ID a,b,c;
    int type,cost;
    long long deps[4];
};

struct REDUCTION
{
    REDUCTION_ID a,b;
    POLY_LIST_ID c;
    PRIMITIVE_ID prim;
    int type,cost;
    long long deps[4];
};

ARRAY<PL,PRIMITIVE_ID> primitives;

ARRAY<REDUCTION,REDUCTION_ID> reductions;

ARRAY<ARRAY<REDUCTION_ID>,POLY_LIST_ID> poly_reductions;

ARRAY<ARRAY<REDUCTION_TEMPLATE>,POLY_LIST_ID> reduction_rules;

VECTOR<int,256> maps[24*24];

void pr_red(REDUCTION_ID i)
{
    const REDUCTION& r=reductions(i);
    printf("(");
    pr(poly_list(r.c)," : ");
    for(int j=0;j<256;j++)
        if(r.deps[j/64]&(1LL<<(j%64))){
            for(int k=0;k<8;k++) if(j&(1<<k)) putchar("abcdrstu"[k]);
            putchar(' ');}

    printf(" <%i> ",r.cost);
    for(int i=0;i<4;i++) printf(" %016llx",r.deps[i]);
    if(r.type==0) printf(" PRIM)");
    else if(r.type==3){
        printf(" REPLACE ");
        pr_red(r.a);
        printf(")");}
    else{
        putchar(' ');
        pr_red(r.a);
        if(r.type==1) printf(" + ");
        else printf(" - ");
        pr_red(r.b);
        printf(")");}
}

long long num_cmp=0;
// 0=equal, 1=r better, 2=s better, 3=neither better
int cmp(const REDUCTION& r,const REDUCTION& s)
{
    num_cmp++;
    int r_worse=r.cost>s.cost;
    int s_worse=s.cost>r.cost;
    for(int j=0;j<4;j++){
        if(s.deps[j]&~r.deps[j]) s_worse=1;
        if(r.deps[j]&~s.deps[j]) r_worse=1;}
    return s_worse+2*r_worse;
}

void Insert_Reduction(const REDUCTION& nr,ARRAY<ARRAY<REDUCTION_ID> >& worklist)
{
//    if(nr.c==0) printf("INSERT\n");
    ARRAY<REDUCTION_ID>& ar=poly_reductions(nr.c);
    bool needed=true;
    bool ch=false;
    for(int i=0;i<ar.m;i++){
        REDUCTION& s=reductions(ar(i));
        if(s.cost==1000) continue;
        int cm=cmp(s,nr);
//        if(nr.c==0) printf("cmp %i\n",cm);
        if(cm<2) needed=false;
        else if(cm==2){s.cost=1000;s.a=reductions.m;s.type=3;ch=true;}}
    PHYSBAM_ASSERT(!ch || needed);
    if(needed){
        REDUCTION_ID in=reductions.Append(nr);
        if(nr.c==0) LOG::printf("add %P\n",in);
        ar.Append(in);
        if(nr.c==0) pr_red(in);
        if(nr.c==0) printf("\n");
        worklist(nr.cost).Append(in);}
    if(ar.m>=10 && 0){
        int n=0;
        for(int i=0;i<ar.m;i++){
            REDUCTION& s=reductions(ar(i));
            if(s.cost==1000) continue;
            n++;}
        if(n>=10){
            puts("DUMPING RULES");
            for(int i=0;i<ar.m;i++){
                pr_red(ar(i));
                puts("");
            }
            PHYSBAM_FATAL_ERROR("Too many reduction rules");}}
}

template<int n>
void Add_Primitive(const VECTOR<int,n>& v,ARRAY<ARRAY<REDUCTION_ID> >& worklist)
{
    PL pl(VECTOR<int,max_pts-n>().Append_Elements(v));

    for(int i=0;i<24*24;i++){
        VECTOR<int,max_pts> u(maps[i].Subset(pl));
        u.Sort();
        PRIMITIVE_ID k=primitives.Append(u);
        POLY_LIST_ID w=poly_lookup.Get(u);
        REDUCTION r={REDUCTION_ID(-1),REDUCTION_ID(-1),w,k,0,1,{}};
        for(int j=0;j<max_pts;j++)
            if(u(j))
                r.deps[u(j)/64]|=(1LL<<(u(j)&63));
        for(int j=0;j<16;j++)
            if(count_bits(j)==3){
                r.deps[0]|=1LL<<j;
                r.deps[j/4]|=(1LL<<(j*16&63));}
        Insert_Reduction(r,worklist);}
}

int main(int argc, char* argv[])
{
    for(int i=0;i<24;i++)
        for(int j=0;j<24;j++){
            VECTOR<int,4> a=permute_four(VECTOR<int,4>(0,1,2,3),i);
            VECTOR<int,4> b=permute_four(VECTOR<int,4>(4,5,6,7),j);
            VECTOR<int,8> M=a.Append_Elements(b);
            for(int k=1;k<256;k++){
                int r=0;
                for(int l=0;l<8;l++)
                    if(k&(1<<l))
                        r|=1<<M(l);
                maps[i*24+j](k)=r;}}

    PL init_volume;
    init_volume(max_pts-1)=14;
    init_volume(max_pts-2)=13;
    init_volume(max_pts-3)=11;
    init_volume(max_pts-4)=7;

    HISTORY init_hist = {init_volume,0,0,{},{},VECTOR<bool,256>(),4,-1,PL(),PL(),max_pts-4,0,0,REDUCTION_ID(-1)};

    for(int i=0;i<16;i++){
        init_hist.safe(i)=true;
        init_hist.safe(16*i)=true;}

    Evolve_History(init_hist);

    LOG::printf("%P\n", final_history.m);

    for(HISTORY_ID i;i<final_history.m;i++) Deduce_Safe(final_history(i));

    int num_tests=0,num_ded_tests=0,num_safe=0;

    for(HISTORY_ID i;i<final_history.m;i++){
        num_tests=std::max(num_tests,final_history(i).num_tests);
        num_ded_tests=std::max(num_ded_tests,final_history(i).num_ded_tests);
        num_safe=std::max(num_safe,final_history(i).safe.Number_True());}

    LOG::printf("%P %P %P\n", num_tests,num_ded_tests,num_safe);
    LOG::printf("%P %P\n",poly_list.m,poly_lookup.Size());

    poly_reductions.Resize(poly_list.m);
    reduction_rules.Resize(poly_list.m);

    ARRAY<ARRAY<REDUCTION_ID> > worklist(100);

    Add_Primitive(VECTOR<int,4>(),worklist);
    Add_Primitive(VECTOR<int,4>(7,11,13,14),worklist);
    Add_Primitive(VECTOR<int,4>(7,19,21,22),worklist);
    Add_Primitive(VECTOR<int,4>(19,35,49,50),worklist);
    Add_Primitive(VECTOR<int,4>(49,81,97,112),worklist);
    Add_Primitive(VECTOR<int,4>(112,176,208,224),worklist);
    Add_Primitive(VECTOR<int,6>(7,11,21,22,25,26),worklist);
    Add_Primitive(VECTOR<int,6>(49,81,112,161,193,224),worklist);
    Add_Primitive(VECTOR<int,6>(7,19,21,38,50,52),worklist);
    Add_Primitive(VECTOR<int,6>(19,21,35,37,50,52),worklist);
    Add_Primitive(VECTOR<int,6>(19,35,49,82,98,112),worklist);
    Add_Primitive(VECTOR<int,6>(19,35,81,82,97,98),worklist);

    HASHTABLE<PL,ARRAY<POLY_LIST_ID> > face_hash;
    for(POLY_LIST_ID i;i<poly_list.m;i++)
        for(int f=0;f<8;f++){
            int fi=0;
            PL face;
            for(int j=0;j<max_pts;j++)
                if(poly_list(i)(j)&(1<<f))
                    face(fi++)=poly_list(i)(j);
            if(fi) face_hash.Get_Or_Insert(face).Append(i);}

    printf("Do faces %i\n",face_hash.Size());
    for(HASHTABLE<PL,ARRAY<POLY_LIST_ID> >::ITERATOR it(face_hash);it.Valid();it.Next()){
        const ARRAY<POLY_LIST_ID>& ar=it.Data();
        const PL& f=it.Key();
        ARRAY<int> masks;

        int face_all=0, face_plane=255;
        for(int i=0;i<max_pts;i++)
            if(f(i)){
                face_all|=f(i);
                face_plane&=f(i);}
        PHYSBAM_ASSERT(face_plane && !(face_plane&(face_plane-1)));

        for(int i=0;i<ar.m;i++){
            int m=0;
            PL& a=poly_list(ar(i));
            for(int j=0;j<max_pts;j++) m|=a(j);
            masks.Append(m);}

        for(int i=0;i<ar.m;i++)
            for(int j=i+1;j<ar.m;j++){
                PL& a=poly_list(ar(i));
                PL& b=poly_list(ar(j));
                if((masks(i)&masks(j))!=face_all) continue;
                VECTOR<bool,256> m;
                m.Subset(a).Fill(true);
                bool bad=false;
                for(int k=0;k<max_pts;k++){
                    int l=b(k);
                    if(l && !(l&face_plane) && m(l)){
                        bad=true;
                        break;}
                    else m(l)=true;}
                if(bad) continue;
                PL nw;
                int nn=0;
                for(int k=1;k<256;k++)
                    if(!(k&face_plane) && m(k))
                        nw(nn++)=k;
                nw.Sort();

                if(POLY_LIST_ID* p=poly_lookup.Get_Pointer(nw)){
                    POLY_LIST_ID ii=ar(i), jj=ar(j), kk=*p;
                    REDUCTION_TEMPLATE ri={kk,jj,ii,2,face_plane<16};
                    REDUCTION_TEMPLATE rj={kk,ii,jj,2,face_plane<16};
                    REDUCTION_TEMPLATE rk={ii,jj,kk,1,face_plane<16};
                    reduction_rules(ii).Append(ri);
                    reduction_rules(jj).Append(rj);
                    reduction_rules(kk).Append(rk);}}}

    int max_cost=4;
    LOG::printf("Consider reductions %P\n",reductions.m);
    for(int c=0;c<worklist.m && c<max_cost;c++){
        for(int w=0;w<worklist(c).m;w++){
            REDUCTION_ID r=worklist(c)(w);
            if(reductions(r).cost>c) continue;
            POLY_LIST_ID poly=reductions(r).c;
            if(w%10000==0){
                LOG::printf("(%P %P : %P) %P : %P (%P %lli)   ",c,w,worklist(c).m,r,reduction_rules(poly).m,reductions.m,num_cmp);
                for(int i=0;i<15;i++) LOG::printf(" %P",worklist(i).m);
                LOG::printf("\n");}
            for(int s=0;s<reduction_rules(poly).m;s++){
                const REDUCTION_TEMPLATE& red=reduction_rules(poly)(s);
                for(int i=0;i<poly_reductions(red.a).m;i++){
                    REDUCTION_ID ra=poly_reductions(red.a)(i);
                    if(reductions(ra).cost>reductions(r).cost) continue;
                    if(reductions(r).cost+reductions(ra).cost>max_cost) continue;
                    REDUCTION nr={r,ra,red.b,PRIMITIVE_ID(),2,reductions(r).cost+reductions(ra).cost};
                    if(red.type==2) exchange(nr.a,nr.b);
                    for(int q=0;q<4;q++) nr.deps[q]=reductions(r).deps[q]|reductions(ra).deps[q];
                    Insert_Reduction(nr,worklist);}
                for(int i=0;i<poly_reductions(red.b).m;i++){
                    REDUCTION_ID rb=poly_reductions(red.b)(i);
                    if(reductions(rb).cost>reductions(r).cost) continue;
                    if(reductions(r).cost+reductions(rb).cost>max_cost) continue;
                    REDUCTION nr={r,rb,red.a,PRIMITIVE_ID(),3-red.type,reductions(r).cost+reductions(rb).cost};
                    for(int q=0;q<4;q++) nr.deps[q]=reductions(r).deps[q]|reductions(rb).deps[q];
                    Insert_Reduction(nr,worklist);}}}}

    // TODO: make reverse sweep through reductions in case later one is better than earlier one

    LOG::printf("Find unreduced %P\n",reductions.m);
    for(HISTORY_ID i;i<final_history.m;i++){
        HISTORY& h=final_history(i);
        if(h.poly(5)) continue;
        puts("===============");
        pr(h.poly,"\n");
        for(int k=0;k<256;k++){
            if(h.safe(k) && count_bits(k)==3){
                for(int j=0;j<8;j++)
                    if(k&(1<<j)) putchar("abcdrstu"[j]);
                putchar(' ');}}
        putchar('\n');

        pr(h);

        POLY_LIST_ID p_id=poly_lookup.Get(h.poly);
        LOG::printf("id %P : %P\n",p_id,poly_reductions(p_id).m);
        for(int j=0;j<poly_reductions(p_id).m;j++){
            REDUCTION_ID ri=poly_reductions(p_id)(j);
            const REDUCTION& r=reductions(ri);
            LOG::printf("cost (%P) %P\n",ri,r.cost);
            pr_red(ri);
            LOG::printf("\n");
            bool ok=true;
            for(int k=0;k<4;k++)
                for(int l=0;l<64;l++)
                    if(r.deps[k]&(1LL<<l))
                        if(!h.safe(64*k+l)){
                            for(int q=0;q<8;q++) if((64*k+l)&(1<<q)) putchar("abcdrstu"[q]);puts("  (UNSAFE)");
                            ok=false;}
            if(!ok) continue;
            PHYSBAM_ASSERT(h.reduction==REDUCTION_ID(-1));
            h.reduction=ri;
            break;}
        if(h.reduction==REDUCTION_ID(-1)){
            LOG::printf("Failed to reduce %P  ",h.poly);
            pr(h.poly,"\n");
            return 0;}}

    return 0;
}

