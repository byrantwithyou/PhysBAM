//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __BLOCK_MESHING_ITERATOR__
#define __BLOCK_MESHING_ITERATOR__
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <list>

namespace PhysBAM{

template<class TV> struct BLOCK_MESHING_ITERATOR;
template<class T>
struct BLOCK_MESHING_ITERATOR<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> IV;
    typedef VECTOR<int,3> IV3;

    const ARRAY<TV>& side0,&side1;
    ARRAY<T> l0,l1;
    ARRAY<TV> X0,X1;
    T h;
    int nseg,k;

    BLOCK_MESHING_ITERATOR(const ARRAY<TV>& s0,const ARRAY<TV>& s1,int n,T dx):
        side0(s0),side1(s1),l0(s0.m),l1(s1.m),X0(s0),h(dx),nseg(n),k(0)
    {
        Init();
        X1=Sample((T)1/nseg);
    }

    void Next()
    {
        k++;
        if(!Valid()) return;
        if(k==nseg-1) X0=side1;
        else X0=Sample((T)(k+1)/nseg);
        X0.Exchange(X1);
    }

    bool Valid() const {return k<nseg;}
    int First_Diagonal_Edge() const {return X0.m-1;}
    int Last_Diagonal_Edge() const {return 2*X0.m+X1.m-3;}

    void Build(ARRAY<TV>& X,ARRAY<IV3>& E,ARRAY<IV>& S,ARRAY<int>& ticks) const
    {
        X.Append_Elements(X0);
        X.Append_Elements(X1);
        auto angle=[X](int v0,int v1,int v2)
        {
            TV u=(X(v2)-X(v1)).Normalized();
            TV v=(X(v0)-X(v1)).Normalized();
            return TV::Angle_Between(u,v);
        };
        auto area=[X](int v0,int v1,int v2){return TRIANGLE_2D<T>(X(v0),X(v1),X(v2)).Signed_Area();};

        int n0=X0.m,n1=X1.m;
        E.Resize(n0+n1-2);
        S.Resize(2*(n0+n1-2)+1);
        ticks.Resize(S.m,use_init,-7);
        int i=0,j=n0,alt=0;
        S(n0-1)=IV(i,j);
        ticks(n0-1)=0;
        while(i<n0-1 || j<n0+n1-1){
            T a0=0;
            if(i+1<n0) a0=angle(j,i+1,i);
            T a1=0;
            if(j+1<n0+n1) a1=angle(j,j+1,i);
            bool sep_i=((i+1<=n0/2 && j<n0+n1/2) || (i+1>n0/2 && j>=n0+n1/2));
            bool sep_j=((i<=n0/2 && j+1<n0+n1/2) || (i>n0/2 && j+1>=n0+n1/2));
            // Choose i+1 side when
            // j+1 is not feasible, or
            // i is feasible, and edge i+1 -- j separates the borrowed dofs, and i+1 is more preferable than j+1.
            //
            // i+1 is more preferable when
            // j+1 side failed to separate dofs.
            // When both sides separate dofs, i+1 forms a less sharp triangle.
            // When there is a tie, alternate the side choice.
            if(j+1>=n0+n1 || (i+1<n0 && area(i,j+1,i+1)<1e-6) ||
                (i+1<n0 && area(j,j+1,i+1)>=1e-6 && sep_i && (!sep_j || (abs(a0-a1)<1e-6 && alt==0) || a0>a1)))
            {
                E(i+j-n0)=IV3(i+1,i,j);
                S(i)=IV(i+1,i);
                S(i+j)=IV(i+1,j);
                ticks(i)=0;
                ticks(i+j)=0;
                i++;
                alt=1;
            }
            else
            {
                E(i+j-n0)=IV3(j+1,i,j);
                S(i+j)=IV(j+1,i);
                S(n0+n1-2+j)=IV(j+1,j);
                ticks(i+j)=1;
                ticks(n0+n1-2+j)=0;
                j++;
                alt=0;
            }
        }
    }

private:
    void Init()
    {
        l0(0)=0;
        for(int j=1;j<side0.m;j++)
            l0(j)=l0(j-1)+(side0(j)-side0(j-1)).Magnitude();
        l1(0)=0;
        for(int j=1;j<side1.m;j++)
            l1(j)=l1(j-1)+(side1(j)-side1(j-1)).Magnitude();
    }

    ARRAY<TV> Sample(T s) const
    {
        auto point=[](int j,T t,const ARRAY<TV>& side)
        {
            if(j>=side.m-1) return side(j);
            TV p=side(j);
            return p+t*(side(j+1)-p);
        };
        auto loc=[point](T u,const ARRAY<TV>& side,const ARRAY<T>& l)
        {
            T dist=u*l.Last();
            auto iter=std::lower_bound(l.begin(),l.end(),dist);
            PHYSBAM_ASSERT(iter!=l.end());
            int j=iter-l.begin();
            if(*iter==dist) return point(j,0,side);
            else
            {
                T cur=*iter,prev=*(iter-1);
                return point(j-1,(dist-prev)/(cur-prev),side);
            }
        };

        std::list<PAIR<T,TV> > verts;
        TV p0=(1-s)*loc(0,side0,l0)+s*loc(0,side1,l1);
        TV p1=(1-s)*loc(1,side0,l0)+s*loc(1,side1,l1);
        verts.insert(verts.end(),{0,p0});
        verts.insert(verts.end(),{1,p1});
        T max_len=(p1-p0).Magnitude();
        while(max_len>1.5*h)
        {
            T len=0;
            for(auto k=verts.begin();k!=verts.end();k++)
            {
                if(k==verts.begin()) continue;
                auto prev=k;
                prev--;
                T u=(prev->x+k->x)*0.5;
                TV v=(1-s)*loc(u,side0,l0)+s*loc(u,side1,l1);
                verts.insert(k,{u,v});
                len=std::max(len,std::max((prev->y-v).Magnitude(),(k->y-v).Magnitude()));
            }
            max_len=len;
        }

        ARRAY<TV> ret;
        for(auto iter=verts.begin();iter!=verts.end();iter++)
            ret.Append(iter->y);
        return ret;
    }
};
}
#endif
