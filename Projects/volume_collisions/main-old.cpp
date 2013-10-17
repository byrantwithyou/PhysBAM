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
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;

#if 0
TV pt_from_planes(VECTOR<PLANE<T>,3> vp)
{
    TV rhs(vp(0).x1.Dot(vp(0).normal),vp(1).x1.Dot(vp(1).normal),vp(2).x1.Dot(vp(2).normal));
    MATRIX<T,TV::m> M(vp(0).normal,vp(1).normal,vp(2).normal);
    return M.Transposed().Solve_Linear_System(rhs);
}

const int max_pts=14;
HASHTABLE<VECTOR<int,max_pts> > polytopes;

VECTOR<int,256> maps[24*24];

bool Add(const VECTOR<int,max_pts>& p)
{
    if(polytopes.Contains(p)) return false;
    for(int i=0;i<24*24;i++){
        VECTOR<int,max_pts> v(maps[i].Subset(p));
        v.Sort();
        polytopes.Set(v);}

    return true;
}

int main(int argc, char* argv[])
{
    RANDOM_NUMBERS<T> rand;

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

    int num_d=0;
    for(int t=0;t<100000000;t++)
    {
        VECTOR<TV,TV::m+1> A,B;
        rand.Fill_Uniform(A,-1,1);
        rand.Fill_Uniform(B,-1,1);

        TETRAHEDRON<T> tA(A),tB(B);
        if(tA.Signed_Size()<0)
        {
            std::swap(tA.X(0),tA.X(1));
            tA.Create_Triangles();
        }
        if(tB.Signed_Size()<0)
        {
            std::swap(tB.X(0),tB.X(1));
            tB.Create_Triangles();
        }

        VECTOR<PLANE<T>,8> planes;
        for(int i=0;i<4;i++) planes(i)=tA.triangle(i).Plane();
        for(int i=0;i<4;i++) planes(i+4)=tB.triangle(i).Plane();

        TV pt=pt_from_planes(VECTOR<PLANE<T>,3>(planes.Subset(VECTOR<int,3>(2,3,6))));

        VECTOR<int,max_pts> volume;
        volume(max_pts-1)=14;
        volume(max_pts-2)=13;
        volume(max_pts-3)=11;
        volume(max_pts-4)=7;

        bool degen=false;
        for(int p=4;p<8;p++)
        {
            int keep_num=0;
            VECTOR<int,max_pts> keep;
            HASHTABLE<int,int> h;
            for(int i=0;i<volume.m;i++)
            {
                int v=volume(i),ii=0;
                if(!v) continue;
                VECTOR<PLANE<T>,3> pl;
                for(int j=0;j<8;j++)
                    if(v&(1<<j))
                        pl(ii++)=planes(j);
                PHYSBAM_ASSERT(ii==3);
                TV X=pt_from_planes(pl);
                T vol=planes(p).Signed_Distance(X);
                if(fabs(vol)<1e-4) degen=true;
                if(vol>0) continue;
                keep(keep_num++)=v;
                for(int j=0;j<8;j++)
                    if(v&(1<<j))
                        h.Get_Or_Insert(v&~(1<<j))++;
            }
            for(HASHTABLE<int,int>::ITERATOR it(h);it.Valid();it.Next())
                if(it.Data()==1)
                    keep(keep_num++)=it.Key()|(1<<p);
            volume=keep;
        }
        num_d+=degen;
        if(degen) continue;
        volume.Sort();
        Add(volume);
    }
    for(HASHTABLE<VECTOR<int,max_pts> >::ITERATOR it(polytopes);it.Valid();it.Next())
        LOG::cout<<it.Key()<<std::endl;
    
    return 0;
}
#else

typedef VECTOR<int,12> PL;
ARRAY<PL> poly;
HASHTABLE<PL,int> poly_lookup;

VECTOR<int,256> maps[24*24];

ARRAY<int> rep;
ARRAY<HASHTABLE<VECTOR<int,2> > > production_triples;

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

int ok_bits=(1<<3)|(1<<5)|(1<<9)|(1<<6)|(1<<10)|(1<<12);

HASHTABLE<int,int> safe_deduction;

int Safe_Deduction(int c,int t)
{
    if(int* r=safe_deduction.Get_Pointer((t<<16)|c)) return *r;
    int status[16]={0};
    for(int i=0;i<16;i++)
        if(c&(1<<i))
            status[i]=1;

    for(int i=0;i<4;i++)
        for(int j=i+1;j<4;j++)
        {
            int k=0;
            while(k==i || k==j) k++;
            int m=6-i-j-k;
            int a=(1<<i)|(1<<j),c=a|(1<<k),d=a|(1<<m);
            if((t&(1<<c)) && (t&(1<<d)))
                status[a]=2;
        }

    for(int r=0;r<4;r++)
        for(int i=0;i<16;i++)
            if(status[i])
                for(int j=i+1;j<16;j++)
                    if(status[j])
                    {
                        int k=i^j^15;
                        if(!(ok_bits&(1<<k))) continue;
                        PHYSBAM_ASSERT(i&j&k);
                        PHYSBAM_ASSERT((i|j|k)==15);
                        if(status[i] && status[i]==status[j]) status[k]=2;
                        if(status[i]+status[j]==3) status[k]=1;
                    }

    int d=0;
    for(int i=0;i<16;i++)
        if(status[i]==1)
            d|=1<<i;

    if((c&d)!=c) d=0;
    printf("deduce: (%i %i) %i\n",c,t,d);
    safe_deduction.Set((t<<16)|c,d);
    return d;
}

ARRAY<VECTOR<bool,256> > safe_masks;

void safe_mask(VECTOR<bool,256>& s,const PL& p)
{
    s.Fill(false);
    s.Subset(p).Fill(true);

    for(int j=4;j<8;j++)
    {
        int c=0,t=0;
        for(int i=1;i<16;i++)
        {
            if(s(i+(1<<j)))
                c|=1<<i;
        }
        for(int i=1;i<16;i++) if(s(i)) t|=1<<i;
        if(c==1576)
        {
            printf("BAD %i ",j);
            pr(p,"\b");
        }
        int d=Safe_Deduction(c,t);
        if(c && !d) s(255)=true;
        for(int i=1;i<16;i++)
            if(d&(1<<i))
                s(i+(1<<j))=true;
    }

    for(int j=0;j<4;j++)
    {
        int c=0,t=0;
        for(int i=1;i<16;i++)
        {
            if(s(16*i+(1<<j)))
                c|=1<<i;
        }
        for(int i=1;i<16;i++) if(s(i*16)) t|=1<<i;
        if(c==1576)
        {
            printf("BAD %i ",j);
            pr(p,"\n");
        }
        int d=Safe_Deduction(c,t);
        if(c && !d) s(255)=true;
        for(int i=1;i<16;i++)
            if(d&(1<<i))
                s(16*i+(1<<j))=true;
    }

    for(int i=0;i<16;i++){s(i)=true;s(i*16)=true;}
}

bool prod_safe(const PL& a,const PL& b,const PL& c)
{
    VECTOR<bool,256> s;
    safe_mask(s,c);
    if(s.Subset(a).Contains(false)) return false;
    if(s.Subset(b).Contains(false)) return false;
    return true;
}

int main(int argc, char* argv[])
{
    puts("// SETUP PERM MAP");
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

    if(0)
    {
        puts("// READ");
        std::string line;
        while(getline(std::cin,line))
        {
            std::stringstream s(line);
            PL v;
            for(int i=0;i<12;i++)
                if(!(s>>v(i)))
                {
                    v(i)=0;
                    break;
                }
            v.Sort();
            poly.Append(v);
        }

        puts("// SORT");
        poly.Sort(LEXICOGRAPHIC_COMPARE());

        FILE_UTILITIES::Write_To_File(STREAM_TYPE(0.0),"list.dat.gz",poly);
    }
    else
    {
        puts("// BINARY READ");
        FILE_UTILITIES::Read_From_File(STREAM_TYPE(0.0),"list.dat.gz",poly);
    }

    puts("// LOOKUPS");
    for(int i=0;i<poly.m;i++) poly_lookup.Set(poly(i),i);

    UNION_FIND<> union_find(poly.m);
    for(int j=0;j<poly.m;j++)
        for(int i=0;i<24*24;i++){
            PL v(maps[i].Subset(poly(j)));
            v.Sort();
            union_find.Union(j,poly_lookup.Get(v));}

    rep.Resize(poly.m,true,true,-1);
    for(int i=0;i<poly.m;i++)
    {
        int p=union_find.Find(i);
        if(rep(p)<0) rep(p)=i;
        rep(i)=rep(p);
    }

    puts("// COMPUTE SAFE MASKS");
    safe_masks.Resize(poly.m);
    for(int i=0;i<poly.m;i++)
    {
        safe_mask(safe_masks(i),poly(i));
        if(safe_masks(i)(255))
        {
            LOG::printf("INVALID CASE %P  ",poly(i));
            pr(poly(i),"\n");
            poly(i).Fill(0);
        }
    }

    puts("// COMPUTE FACES");
    HASHTABLE<PL,ARRAY<int> > face_hash;
    for(int i=0;i<poly.m;i++)
    {
        for(int f=0;f<8;f++)
        {
            int fi=0;
            PL face;
            for(int j=0;j<12;j++)
                if(poly(i)(j)&(1<<f))
                    face(fi++)=poly(i)(j);
            if(fi) face_hash.Get_Or_Insert(face).Append(i);
        }
    }

    puts("// FIND GEN TRIPLES");
    production_triples.Resize(poly.m);
    for(HASHTABLE<PL,ARRAY<int> >::ITERATOR it(face_hash);it.Valid();it.Next())
    {
        const ARRAY<int>& ar=it.Data();
        const PL& f=it.Key();
        ARRAY<int> masks;

        int face_all=0, face_plane=255;
        for(int i=0;i<12;i++)
            if(f(i))
            {
                face_all|=f(i);
                face_plane&=f(i);
            }
        if(!face_plane || (face_plane&(face_plane-1))) continue;

        for(int i=0;i<ar.m;i++)
        {
            int m=0;
            PL& a=poly(ar(i));
            for(int j=0;j<12;j++) m|=a(j);
            masks.Append(m);
        }
//        LOG::printf("%P -> %P %P %P %P\n",f,ar,masks,face_all,face_plane);

        for(int i=0;i<ar.m;i++)
            for(int j=i+1;j<ar.m;j++)
            {
                PL& a=poly(ar(i));
                PL& b=poly(ar(j));
                if((masks(i)&masks(j))!=face_all) continue;
                VECTOR<bool,256> m;
                m.Subset(a).Fill(true);
                bool bad=false;
                for(int k=0;k<12;k++)
                {
                    if(b(k) && !(b(k)&face_plane) && m(b(k)))
                    {
                        bad=true;
                        break;
                    }
                    else m(b(k))=1;
                }
                if(bad) continue;
                PL nw;
                int nn=0;
                for(int k=1;k<256;k++)
                    if(!(k&face_plane) && m(k))
                        nw(nn++)=k;
                nw.Sort();

                if(int* p=poly_lookup.Get_Pointer(nw))
                {
                    int ii=rep(ar(i)), jj=rep(ar(j)), kk=rep(*p);
                    if(prod_safe(poly(ar(i)),poly(ar(j)),nw))
                    {
                        production_triples(ii).Set(VECTOR<int,2>(jj,kk));
                        production_triples(jj).Set(VECTOR<int,2>(ii,kk));
                    }
                    if(prod_safe(nw,poly(ar(i)),poly(ar(j))))
                    {
                        production_triples(ii).Set(VECTOR<int,2>(kk,jj));
                        production_triples(kk).Set(VECTOR<int,2>(ii,jj));
                    }
                    if(prod_safe(poly(ar(j)),nw,poly(ar(i))))
                    {
                        production_triples(jj).Set(VECTOR<int,2>(kk,ii));
                        production_triples(kk).Set(VECTOR<int,2>(jj,ii));
                    }
                    // LOG::printf("// %P : ", VECTOR<int,3>(ii,jj,kk).Sorted());
                    // pr(a," + ");
                    // pr(b," -> ");
                    // pr(nw,"\n");
                }
            }
    }

    // for(int i=0;i<poly.m;i++)
    //     if(rep(i)==i)
    //         LOG::printf("%P\n",poly(i));

    puts("// GENERATE");
    ARRAY<ARRAY<int> > worklist(100);
    ARRAY<int> state(poly.m);
    ARRAY<int> cost(poly.m);
    ARRAY<VECTOR<int,2> > path(poly.m);
    PL init;
    {
        init.Fill(0);
        init(0)=7;
        init(1)=11;
        init(2)=13;
        init(3)=14;
        init.Sort();
        int i0=rep(poly_lookup.Get(init));
        worklist(1).Append(i0);
        state(i0)=1;
        cost(i0)=1;
    }

    {
        init.Fill(0);
        init(0)=7;
        init(1)=19;
        init(2)=21;
        init(3)=22;
        init.Sort();
        int i0=rep(poly_lookup.Get(init));
        worklist(1).Append(i0);
        state(i0)=1;
        cost(i0)=1;
    }

    {
        init.Fill(0);
        init(0)=19;
        init(1)=35;
        init(2)=49;
        init(3)=50;
        init.Sort();
        int i0=rep(poly_lookup.Get(init));
        worklist(1).Append(i0);
        state(i0)=1;
        cost(i0)=1;
    }

    {
        init.Fill(0);
        init(0)=49;
        init(1)=81;
        init(2)=97;
        init(3)=112;
        init.Sort();
        int i0=rep(poly_lookup.Get(init));
        worklist(1).Append(i0);
        state(i0)=1;
        cost(i0)=1;
    }

    {
        init.Fill(0);
        init(0)=112;
        init(1)=176;
        init(2)=208;
        init(3)=224;
        init.Sort();
        int i0=rep(poly_lookup.Get(init));
        worklist(1).Append(i0);
        state(i0)=1;
        cost(i0)=1;
    }

    {
        init.Fill(0);
        init(0)=7;
        init(1)=11;
        init(2)=21;
        init(3)=22;
        init(4)=25;
        init(5)=26;
        init.Sort();
        int i0=rep(poly_lookup.Get(init));
        worklist(1).Append(i0);
        state(i0)=1;
        cost(i0)=1;
    }

    {
        init.Fill(0);
        init(0)=49;
        init(1)=81;
        init(2)=112;
        init(3)=161;
        init(4)=193;
        init(5)=224;
        init.Sort();
        int i0=rep(poly_lookup.Get(init));
        worklist(1).Append(i0);
        state(i0)=1;
        cost(i0)=1;
    }

    {
        init.Fill(0);
        init(0)=7;
        init(1)=19;
        init(2)=21;
        init(3)=38;
        init(4)=50;
        init(5)=52;
        init.Sort();
        int i0=rep(poly_lookup.Get(init));
        worklist(1).Append(i0);
        state(i0)=1;
        cost(i0)=1;
    }

    {
        init.Fill(0);
        init(0)=19;
        init(1)=21;
        init(2)=35;
        init(3)=37;
        init(4)=50;
        init(5)=52;
        init.Sort();
        int i0=rep(poly_lookup.Get(init));
        worklist(1).Append(i0);
        state(i0)=1;
        cost(i0)=1;
    }

    {
        init.Fill(0);
        init(0)=19;
        init(1)=35;
        init(2)=49;
        init(3)=82;
        init(4)=98;
        init(5)=112;
        init.Sort();
        int i0=rep(poly_lookup.Get(init));
        worklist(1).Append(i0);
        state(i0)=1;
        cost(i0)=1;
    }

    {
        init.Fill(0);
        init(0)=19;
        init(1)=35;
        init(2)=81;
        init(3)=82;
        init(4)=97;
        init(5)=98;
        init.Sort();
        int i0=rep(poly_lookup.Get(init));
        worklist(1).Append(i0);
        state(i0)=1;
        cost(i0)=1;
    }

    for(int c=0;c<worklist.m;c++)
    {
        while(worklist(c).m)
        {
            int a=worklist(c).Pop();
            PHYSBAM_ASSERT(rep(a)==a);
            PHYSBAM_ASSERT(cost(a)==c);
            if(state(a)==2) continue;
            state(a)=2;

            for(HASHTABLE<VECTOR<int,2> >::ITERATOR it(production_triples(a));it.Valid();it.Next())
            {
                int r=it.Key()(0);
                int s=it.Key()(1);
                if(state(r))
                {
                    if(cost(s)==0 || cost(s)>cost(r)+cost(a))
                    {
                        cost(s)=cost(r)+cost(a);
                        path(s)=VECTOR<int,2>(a,r).Sorted();
                        PHYSBAM_ASSERT(state(s)<2);
                        state(s)=1;
                        worklist(cost(s)).Append(s);
                    }
                }
            }
        }
    }
    printf("digraph dg {\n");
    for(int i=0;i<path.m;i++)
    {
        if(path(i).x) printf("{ n%i n%i } -> n%i;\n",path(i).x,path(i).y,i);
        if(state(i)==2)
        {
            printf("n%i [label=\"", i);
            PL p=poly(i);
            for(int j=0;j<12;j++)
                if(p(j))
                {
                    for(int i=0;i<8;i++)
                        if(p(j)&(1<<i))
                            putchar("abcdrstu"[i]);
                    if(j!=11) putchar(' ');
                }
            printf("\"];\n");
        }
    }
    printf("};\n");

    for(int i=0;i<poly.m;i++)
        if(rep(i)==i && state(i)==0)
        {
            LOG::printf("// no gen %i %P  ",i,poly(i));
            pr(poly(i),"\n");
        }

    return 0;
}

#endif




