//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Math_Tools/pow.h>
#include <Core/Matrices/ROTATION.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Images/EPS_FILE.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Seeding/POISSON_DISK.h>
#include <string>
#include "ANALYTIC_FEM.h"
#include "CACHED_ELIMINATION_MATRIX.h"
#include "COMPONENT_LAYOUT_FEM.h"
#include "ELIMINATION_FEM.h"
#include "LAYOUT_BUILDER_FEM.h"
#include "MATRIX_CONSTRUCTION_FEM.h"
#include "EXECUTE_HELPER.h"
#include <boost/polygon/voronoi.hpp>

using namespace PhysBAM;

typedef float RW;
typedef double T;

int cache_size=-1;
std::string cache_pattern="cache-%d.txt";

template<int d>
void Solve_And_Check(COMPONENT_LAYOUT_FEM<T>& cl,const LAYOUT_BUILDER_FEM<T>& builder,
    T mu,const std::string& au,const std::string& ap)
{
    typedef VECTOR<T,d> TV;
    cl.Update_Masters();
    cl.Merge_Blocks();
    CACHED_ELIMINATION_MATRIX<T> cem;
    cem.quiet=true;
    cl.Compute_Dof_Pairs();
    cl.Fill_Reference_Ticks();
    MATRIX_CONSTRUCTION_FEM<TV> mc(cl,"canonical-%d",cache_size);
    mc.mu=mu*cl.unit_kg*pow<2-d>(cl.unit_m)/cl.unit_s;
    mc.Compute_Matrix_Blocks();
    ANALYTIC_FEM<TV> an(mc);
    an.Set_Velocity(au.c_str());
    an.Set_Pressure(ap.c_str());
    an.Compute_RHS();
    mc.Copy_To_CEM(cem);
    mc.Transform_Solution(cem,true,true);
    ELIMINATION_FEM<T> el(cl);
    el.Eliminate_Irregular_Blocks(cem);
    el.Eliminate_Non_Seperators(cem);
    cem.Full_Reordered_Elimination();
    cem.Back_Solve();
    cem.matrix_cache.Init(cache_pattern,cache_size);
    Execute_Helper(&cem,1);
    mc.Transform_Solution(cem,false,false);
    if(!an.Check_Analytic_Solution(false,false))
        printf("%s\n",builder.To_String().c_str());
}

template<int d>
void Test_Pipe(RANDOM_NUMBERS<T>& rng,T mu,T s,T m,T kg,const std::string& au,const std::string& ap)
{
    typedef VECTOR<T,2> TV;
    COMPONENT_LAYOUT_FEM<T> cl;
    cl.unit_m=m;
    cl.unit_s=s;
    cl.unit_kg=kg;
    LAYOUT_BUILDER_FEM<T> builder(cl);
    {
        T a=rng.Get_Uniform_Number(-1,0),b=rng.Get_Uniform_Number(1,2);
        T y=rng.Get_Uniform_Number(-1,1);
        builder.Set_Target_Length(0.25);
        builder.Set_Depth(1,2);
        auto cs=builder.Cross_Section(4,1.0);
        auto v0=builder.Vertex(TV(a,y));
        auto v1=builder.Vertex(TV(b,-y));
        auto bc0=builder.Set_BC(cs,v0,v1,1.0);
        auto bc1=builder.Set_BC(cs,v1,v0,TV());
        builder.Pipe(cs,bc0.x,bc1.x);
    }
    Solve_And_Check<d>(cl,builder,mu,au,ap);
}

template<int d>
void Test_Var_Size_Pipe(RANDOM_NUMBERS<T>& rng,T mu,T s,T m,T kg,const std::string& au,const std::string& ap)
{
    typedef VECTOR<T,2> TV;
    COMPONENT_LAYOUT_FEM<T> cl;
    cl.unit_m=m;
    cl.unit_s=s;
    cl.unit_kg=kg;
    LAYOUT_BUILDER_FEM<T> builder(cl);
    {
        T a=rng.Get_Uniform_Number(-2,-1),b=rng.Get_Uniform_Number(2,3);
        T y=rng.Get_Uniform_Number(-1,1);
        T offset=rng.Get_Uniform_Number(0.5,1),len=rng.Get_Uniform_Number(0.5,1.5);
        builder.Set_Target_Length(0.25);
        builder.Set_Depth(1,2);
        auto cs0=builder.Cross_Section(4,1.0);
        auto cs1=builder.Cross_Section(6,1.5);
        auto v0=builder.Vertex(TV(a,y));
        auto v1=builder.Vertex(TV(b,-y));
        auto bc0=builder.Set_BC(cs0,v0,v1,1.0);
        auto bc1=builder.Set_BC(cs1,v1,v0,TV());
        auto p=builder.Pipe(cs0,cs1,v0,v1,offset,len);
        builder.Pipe(cs0,bc0.x,p(0));
        builder.Pipe(cs1,bc1.x,p(1));
    }
    Solve_And_Check<d>(cl,builder,mu,au,ap);
}

template<int d>
void Test_Joint2(T angle,T mu,T s,T m,T kg,const std::string& au,const std::string& ap)
{
    typedef VECTOR<T,2> TV;
    COMPONENT_LAYOUT_FEM<T> cl;
    cl.unit_m=m;
    cl.unit_s=s;
    cl.unit_kg=kg;
    LAYOUT_BUILDER_FEM<T> builder(cl);
    {
        builder.Set_Target_Length(0.25);
        builder.Set_Depth(1,2);
        auto cs=builder.Cross_Section(4,1.0);
        auto v0=builder.Vertex(TV());
        TV dir(1,2);
        auto v1=builder.Vertex(dir);
        ROTATION<TV> R=ROTATION<TV>::From_Angle(angle);
        auto v2=builder.Vertex(R.Rotate(dir));
        auto bc0=builder.Set_BC(cs,v1,v0,1.0);
        auto bc1=builder.Set_BC(cs,v2,v0,TV(1,1));
        auto j=builder.Joint(cs,2,v0,{v1,v2});
        builder.Pipe(cs,bc0.x,j(0));
        builder.Pipe(cs,bc1.x,j(1));
    }
    Solve_And_Check<d>(cl,builder,mu,au,ap);
}

template<int d>
void Test_Joint3(T a0,T a1,T mu,T s,T m,T kg,const std::string& au,const std::string& ap)
{
    typedef VECTOR<T,2> TV;
    COMPONENT_LAYOUT_FEM<T> cl;
    cl.unit_m=m;
    cl.unit_s=s;
    cl.unit_kg=kg;
    LAYOUT_BUILDER_FEM<T> builder(cl);
    {
        builder.Set_Target_Length(0.25);
        builder.Set_Depth(1,2);
        auto cs=builder.Cross_Section(4,1.0);
        auto v0=builder.Vertex(TV());
        TV dir(1,2);
        auto v1=builder.Vertex(dir);
        ROTATION<TV> R0=ROTATION<TV>::From_Angle(a0);
        auto v2=builder.Vertex(R0.Rotate(dir));
        ROTATION<TV> R1=ROTATION<TV>::From_Angle(a0+a1);
        auto v3=builder.Vertex(R1.Rotate(dir));
        auto bc0=builder.Set_BC(cs,v1,v0,1.0);
        auto bc1=builder.Set_BC(cs,v2,v0,TV(1,1));
        auto bc2=builder.Set_BC(cs,v3,v0,0.5);
        auto j=builder.Joint(cs,3,v0,{v1,v2,v3});
        builder.Pipe(cs,bc0.x,j(0));
        builder.Pipe(cs,bc1.x,j(1));
        builder.Pipe(cs,bc2.x,j(2));
    }
    Solve_And_Check<d>(cl,builder,mu,au,ap);
}

template<int d>
void Test_Joint4_Right_Angle(T angle,T mu,T s,T m,T kg,const std::string& au,const std::string& ap)
{
    typedef VECTOR<T,2> TV;
    COMPONENT_LAYOUT_FEM<T> cl;
    cl.unit_m=m;
    cl.unit_s=s;
    cl.unit_kg=kg;
    LAYOUT_BUILDER_FEM<T> builder(cl);
    {
        builder.Set_Target_Length(0.25);
        builder.Set_Depth(1,2);
        auto cs=builder.Cross_Section(4,1.0);
        ROTATION<TV> R=ROTATION<TV>::From_Angle(angle);
        auto O=builder.Vertex(TV());
        TV v0=R.Rotate(TV(2,0));
        auto east=builder.Vertex(v0);
        auto north=builder.Vertex(v0.Rotate_Counterclockwise_90());
        auto west=builder.Vertex(-v0);
        auto south=builder.Vertex(v0.Rotate_Clockwise_90());
        auto bc_east=builder.Set_BC(cs,east,O,1.0);
        auto bc_north=builder.Set_BC(cs,north,O,0.8);
        auto bc_south=builder.Set_BC(cs,south,O,TV(1,1));
        auto bc_west=builder.Set_BC(cs,west,O,TV());
        auto j=builder.Joint_4_Right_Angle(cs,O,east);
        builder.Pipe(cs,bc_east.x,j(0));
        builder.Pipe(cs,bc_north.x,j(1));
        builder.Pipe(cs,bc_west.x,j(2));
        builder.Pipe(cs,bc_south.x,j(3));
    }
    Solve_And_Check<d>(cl,builder,mu,au,ap);
}

struct CROSS_DATA
{
    using VERT_ID=LAYOUT_BUILDER_FEM<T>::VERT_ID;
    using CID=LAYOUT_BUILDER_FEM<T>::CONNECTOR_ID;

    VERT_ID vid;
    VECTOR<PAIR<VERT_ID,CID>,4> c; // CID*4: +x,+y,-x,-y
    CROSS_DATA()
    {
        vid=VERT_ID(-1);
        for(int i=0;i<4;i++) c(i)={VERT_ID(-1),CID(-1)};
    }
};

// http://random.groverlab.org/candidate/59378553.html
static VECTOR<VECTOR<T,2>,2> grid0[]=
{
    {{2,9},{2,8}},{{5,9},{5,8}},
    {{0,8},{2,8}},{{2,8},{3,8}},{{3,8},{4,8}},{{5,8},{6,8}},
    {{3,8},{3,7}},{{4,8},{4,7}},{{6,8},{6,7}},{{7,8},{7,7}},
    {{0,7},{3,7}},{{4,7},{6,7}},{{6,7},{7,7}},
    {{0,7},{0,6}},{{3,7},{3,6}},{{4,7},{4,6}},{{6,7},{6,6}},{{7,7},{7,4}},
    {{0,6},{1,6}},{{1,6},{3,6}},{{3,6},{4,6}},{{5,6},{6,6}},
    {{0,6},{0,5}},{{1,6},{1,5}},
    {{0,5},{1,5}},{{1,5},{2,5}},{{2,5},{3,5}},{{3,5},{4,5}},
    {{0,5},{0,4}},{{1,5},{1,4}},{{2,5},{2,4}},{{3,5},{3,4}},{{4,5},{4,4}},{{5,5},{5,4}},{{6,5},{6,4}},
    {{0,4},{1,4}},{{2,4},{3,4}},{{4,4},{5,4}},{{5,4},{6,4}},{{6,4},{7,4}},
    {{0,4},{0,3}},{{1,4},{1,3}},{{2,4},{2,3}},{{3,4},{3,2}},{{4,4},{4,3}},{{5,4},{5,3}},{{6,4},{6,3}},{{7,4},{7,2}},
    {{0,3},{1,3}},{{1,3},{2,3}},{{5,3},{6,3}},
    {{5,3},{5,2}},{{6,3},{6,2}},
    {{0,2},{1,2}},{{2,2},{3,2}},{{3,2},{4,2}},{{5,2},{6,2}},{{6,2},{7,2}},
    {{0,2},{0,0}},{{1,2},{1,1}},{{2,2},{2,1}},{{5,2},{5,1}},{{6,2},{6,1}},{{7,2},{7,1}},
    {{1,1},{2,1}},{{2,1},{3.5,1}},{{3.5,1},{4,1}},{{6,1},{7,1}},
    {{3.5,1},{3.5,0}},{{7,1},{7,0}}
};

// http://random.groverlab.org/candidate/11460287.html
static VECTOR<VECTOR<T,2>,2> grid1[]=
{
    {{2,9},{2,8}},{{5,9},{5,8}},
    {{0,8},{1,8}},{{2,8},{3,8}},{{3,8},{5,8}},{{5,8},{6,8}},{{6,8},{7,8}},
    {{0,8},{0,7}},{{1,8},{1,7}},{{2,8},{2,7}},{{3,8},{3,7}},{{5,8},{5,7}},{{6,8},{6,7}},{{7,8},{7,7}},
    {{0,7},{1,7}},{{1,7},{2,7}},{{2,7},{3,7}},{{3,7},{4,7}},{{5,7},{6,7}},{{6,7},{7,7}},
    {{1,7},{1,6}},{{2,7},{2,6}},{{3,7},{3,5}},{{4,7},{4,5}},{{5,7},{5,6}},{{6,7},{6,6}},
    {{0,6},{1,6}},{{1,6},{2,6}},{{5,6},{6,6}},
    {{0,6},{0,4}},{{1,6},{1,5}},{{2,6},{2,5}},{{5,6},{5,5}},{{6,6},{6,5}},{{7,6},{7,4}},
    {{1,5},{2,5}},{{2,5},{3,5}},{{4,5},{5,5}},
    {{1,5},{1,4}},{{2,5},{2,4}},{{4,5},{4,4}},{{5,5},{5,4}},
    {{0,4},{1,4}},{{1,4},{2,4}},{{2,4},{4,4}},{{4,4},{5,4}},{{5,4},{6,4}},{{6,4},{7,4}},
    {{1,4},{1,3}},{{5,4},{5,3}},{{6,4},{6,3}},{{7,4},{7,3}},
    {{0,3},{1,3}},{{2,3},{3,3}},{{5,3},{6,3}},{{6,3},{7,3}},
    {{0,3},{0,2}},{{3,3},{3,2}},{{4,3},{4,2}},{{5,3},{5,1}},{{6,3},{6,2}},{{7,3},{7,2}},
    {{0,2},{1,2}},{{1,2},{3,2}},{{3,2},{4,2}},{{6,2},{7,2}},
    {{1,2},{1,1}},{{3,2},{3,1}},{{7,2},{7,1}},
    {{0,1},{1,1}},{{1,1},{3,1}},{{3,1},{3.5,1}},{{3.5,1},{4,1}},{{5,1},{7,1}},
    {{0,1},{0,0}},{{3.5,1},{3.5,0}},{{7,1},{7,0}}
};

void Generate_Random_Ortho_Grid(T scale,const VECTOR<VECTOR<T,2>,2>* edges,int ne,T mu,T s,T m,T kg,const char* filename)
{
    T w=0.004*scale;
    T dx=(0.025+w)*scale;

    typedef VECTOR<T,2> TV;
    using VERT_ID=LAYOUT_BUILDER_FEM<T>::VERT_ID;
    using CID=LAYOUT_BUILDER_FEM<T>::CONNECTOR_ID;

    COMPONENT_LAYOUT_FEM<T> cl;
    cl.unit_m=m;
    cl.unit_s=s;
    cl.unit_kg=kg;
    LAYOUT_BUILDER_FEM<T> builder(cl);
    builder.Set_Target_Length(w/2);
    builder.Set_Depth(w,1);
    auto cs=builder.Cross_Section(2,w);

    HASHTABLE<TV,CROSS_DATA> crs;
    T maxy=0,miny=1e6;
    for(int i=0;i<ne;i++)
    {
        const auto& e=edges[i];
        maxy=max(maxy,e(0).y,e(1).y);
        miny=min(miny,e(0).y,e(1).y);
        auto pr0=crs.Insert(e(0),{});
        if(pr0.y) pr0.x->vid=builder.Vertex({e(0).x*dx,e(0).y*dx});
        auto pr1=crs.Insert(e(1),{});
        if(pr1.y) pr1.x->vid=builder.Vertex({e(1).x*dx,e(1).y*dx});

        int a=-1,b=-1;
        if(e(0).x==e(1).x)
        {
            a=e(1).y>e(0).y?1:3;
            b=4-a;
        }
        else
        {
            a=e(1).x>e(0).x?0:2;
            b=2-a;
        }
        pr0.x->c(a).x=pr1.x->vid;
        pr1.x->c(b).x=pr0.x->vid;
    }

    for(auto& v:crs)
    {
        ARRAY<VERT_ID> arms;
        for(int i=0;i<4;i++) if(v.data.c(i).x>=VERT_ID())
            arms.Append(v.data.c(i).x);
        if(arms.m==0) continue;
        if(arms.m==1)
        {
            if(v.key.y==miny)
            {
                auto sink=builder.Set_BC(cs,v.data.vid,arms(0),TV());
                for(int i=0;i<4;i++) if(v.data.c(i).x>=VERT_ID())
                    v.data.c(i).y=sink.x;
            }
            else
            {
                T flowrate=0;
                if(v.key.y==maxy) flowrate=0.005;
                auto src=builder.Set_BC(cs,v.data.vid,arms(0),flowrate);
                for(int i=0;i<4;i++) if(v.data.c(i).x>=VERT_ID())
                    v.data.c(i).y=src.x;
            }
        }
        else if(arms.m<4)
        {
            auto j=builder.Joint(cs,arms.m,v.data.vid,arms);
            for(int i=0,k=0;i<4;i++) if(v.data.c(i).x>=VERT_ID())
                v.data.c(i).y=j(k++);
        }
        else if(arms.m==4)
        {
            auto j=builder.Joint_4_Right_Angle(cs,v.data.vid,arms(0));
            for(int i=0;i<4;i++)
                v.data.c(i).y=j(i);
        }
        else PHYSBAM_FATAL_ERROR();
    }

    for(int i=0;i<ne;i++)
    {
        const auto& e=edges[i];
        int a=-1,b=-1;
        if(e(0).x==e(1).x)
        {
            a=e(1).y>e(0).y?1:3;
            b=4-a;
        }
        else
        {
            a=e(1).x>e(0).x?0:2;
            b=2-a;
        }
        const auto& v0=crs.Get(e(0));
        const auto& v1=crs.Get(e(1));
        builder.Pipe(cs,v0.c(a).y,v1.c(b).y);
    }
    FILE* file=fopen(LOG::sprintf("%s.txt",filename).c_str(),"w");
    fprintf(file,"%s\n",builder.To_String().c_str());
}

void Generate_Grid(T scale,T mu,T s,T m,T kg,int n,const char* filename)
{
    typedef VECTOR<T,2> TV;
    T w=1*scale;
    T dx=w+3*w;

    COMPONENT_LAYOUT_FEM<T> cl;
    cl.unit_m=m;
    cl.unit_s=s;
    cl.unit_kg=kg;
    LAYOUT_BUILDER_FEM<T> builder(cl);
    builder.Set_Target_Length(0.5*scale);
    builder.Set_Depth(w,1);
    auto cs=builder.Cross_Section(2,w);

    using VERT_ID=LAYOUT_BUILDER_FEM<T>::VERT_ID;
    using CID=LAYOUT_BUILDER_FEM<T>::CONNECTOR_ID;

    ARRAY<PAIR<VERT_ID,VECTOR<CID,4> > > verts(n*n); // CID*4: +x,+y,-x,-y
    VERT_ID src=builder.Vertex(TV(-dx,(n-1)*dx)),sink=builder.Vertex(TV(n*dx,0));
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            verts(i*n+j)={builder.Vertex(TV(i*dx,j*dx)),{CID(-1),CID(-1),CID(-1),CID(-1)}};

    int dirs[4][2]=
    {
        {1,0},
        {0,1},
        {-1,0},
        {0,-1}
    };
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
        {
            auto& v=verts(i*n+j);
            ARRAY<VERT_ID> arms;
            for(int k=0;k<4;k++)
            {
                int nx=i+dirs[k][0],ny=j+dirs[k][1];
                if(nx>=0 && nx<n && ny>=0 && ny<n)
                    arms.Append(verts(nx*n+ny).x);
                else if(i==0 && j==n-1 && k==2)
                    arms.Append(src);
                else if(i==n-1 && j==0 && k==0)
                    arms.Append(sink);
            }
            ARRAY<CID> jt;
            if(arms.m==4)
                jt=builder.Joint_4_Right_Angle(cs,v.x,arms(0));
            else
                jt=builder.Joint(cs,arms.m,v.x,arms);
            for(int k=0,b=0;k<4;k++)
            {
                int nx=i+dirs[k][0],ny=j+dirs[k][1];
                if(nx>=0 && nx<n && ny>=0 && ny<n)
                    v.y(k)=jt(b++);
                else if((i==0 && j==n-1 && k==2) || (i==n-1 && j==0 && k==0))
                    v.y(k)=jt(b++);
            }
        }

    for(int i=1;i<n;i++)
        for(int j=0;j<n;j++)
            builder.Pipe(cs,verts((i-1)*n+j).y(0),verts(i*n+j).y(2));
    for(int i=0;i<n;i++)
        for(int j=1;j<n;j++)
            builder.Pipe(cs,verts(i*n+j-1).y(1),verts(i*n+j).y(3));

    auto bc_src=builder.Set_BC(cs,src,verts(n-1).x,0.005);
    auto bc_sink=builder.Set_BC(cs,sink,verts((n-1)*n).x,TV());
    builder.Pipe(cs,bc_src.x,verts(n-1).y(2));
    builder.Pipe(cs,verts((n-1)*n).y(0),bc_sink.x);

    FILE* file=fopen(LOG::sprintf("%s.txt",filename).c_str(),"w");
    fprintf(file,"%s\n",builder.To_String().c_str());
}

template<class VERT_ID,class CID,class T,class CS_ID>
std::tuple<VERT_ID,CID,VERT_ID,CID> Coil(LAYOUT_BUILDER_FEM<T>& builder,CS_ID cs,const VECTOR<T,2>& O,T l,T h)
{
    typedef VECTOR<T,2> TV;
    T dy=h/6,dx=l/6;
    LOG::printf("X=%P\n",O(0));
    LOG::printf("X=%P\n",O(0)+3*dx);
    LOG::printf("X=%P\n",O(0)+6*dx);

    LOG::printf("Y=%P\n",O(1)+3*dy);
    LOG::printf("Y=%P\n",O(1));
    LOG::printf("Y=%P\n",O(1)-3*dy);
    VERT_ID w[12]=
    {
        builder.Vertex(O+TV(0,2*dy)),
        builder.Vertex(O+TV(0,3*dy)),
        builder.Vertex(O+TV(dx,3*dy)),

        builder.Vertex(O+TV(2*dx,3*dy)),
        builder.Vertex(O+TV(3*dx,3*dy)),
        builder.Vertex(O+TV(3*dx,2*dy)),

        builder.Vertex(O+TV(3*dx,-2*dy)),
        builder.Vertex(O+TV(3*dx,-3*dy)),
        builder.Vertex(O+TV(4*dx,-3*dy)),

        builder.Vertex(O+TV(5*dx,-3*dy)),
        builder.Vertex(O+TV(6*dx,-3*dy)),
        builder.Vertex(O+TV(6*dx,-2*dy)),
    };
    VECTOR<CID,2> jt[4];
    for(int i=0;i<4;i++)
    {
        int j=1+3*i;
        jt[i]=builder.Joint(cs,2,w[j],{w[j-1],w[j+1]});
    }
    for(int i=1;i<4;i++)
        builder.Pipe(cs,jt[i-1](1),jt[i](0));
    return std::make_tuple(w[0],jt[0](0),w[11],jt[3](1));
};

void Generate_Simple(T scale,T mu,T s,T m,T kg,const char* filename)
{
    typedef VECTOR<T,2> TV;
    T w=1*scale;
    COMPONENT_LAYOUT_FEM<T> cl;
    cl.unit_m=m;
    cl.unit_s=s;
    cl.unit_kg=kg;
    T flowrate=0.005;
    LAYOUT_BUILDER_FEM<T> builder(cl);
    builder.Set_Target_Length(0.5*scale);
    builder.Set_Depth(w,1);
    auto cs=builder.Cross_Section(2,w);
    auto cs1=builder.Cross_Section(6,w*7);
    LOG::printf("Y=%P\n",3.5*w);;
    LOG::printf("Y=%P\n",-3.5*w);;

    using VERT_ID=LAYOUT_BUILDER_FEM<T>::VERT_ID;
    using CID=LAYOUT_BUILDER_FEM<T>::CONNECTOR_ID;

    LOG::printf("BEGIN simple\n");
    T x=0;
    LOG::printf("X=%P\n",x*scale);
    LOG::printf("Y=%P\n",2*scale);
    LOG::printf("Y=%P\n",-2*scale);
    VERT_ID v0=builder.Vertex(TV(x,2)*scale);
    VERT_ID v2=builder.Vertex(TV(x,-2)*scale);
    x+=2;
    LOG::printf("X=%P\n",x*scale);
    VERT_ID v1=builder.Vertex(TV(x,2)*scale);
    VERT_ID v3=builder.Vertex(TV(x,-2)*scale);
    x+=2;
    LOG::printf("X=%P\n",x*scale);
    LOG::printf("Y=%P\n",0);
    VERT_ID v4=builder.Vertex(TV(x,0)*scale);
    auto src0=builder.Set_BC(cs,v0,v1,flowrate);
    auto src2=builder.Set_BC(cs,v2,v3,flowrate);
    auto j1=builder.Joint(cs,2,v1,{v0,v4});
    auto j3=builder.Joint(cs,2,v3,{v2,v4});
    x+=3;
    LOG::printf("X=%P\n",x*scale);
    VERT_ID v5=builder.Vertex(TV(x,0)*scale);
    auto c0=Coil<VERT_ID,CID>(builder,cs,TV(x,0)*scale,4.5*scale,4.5*scale);
    x+=4.5;
    LOG::printf("X=%P\n",x*scale);
    auto j4=builder.Joint(cs,3,v4,{v1,v3,v5});
    auto j5=builder.Joint(cs,2,v5,{v4,std::get<0>(c0)});
    auto c1=Coil<VERT_ID,CID>(builder,cs,TV(x,0)*scale,4.5*scale,4.5*scale);
    x+=4.5;
    LOG::printf("X=%P\n",x*scale);
    VERT_ID v6=builder.Vertex(TV(x,0)*scale);
    x+=6;
    LOG::printf("X=%P\n",(x-6+2)*scale);
    LOG::printf("X=%P\n",(x-6+5)*scale);
    VERT_ID v7=builder.Vertex(TV(x,0)*scale);
    auto j6=builder.Joint(cs,2,v6,{std::get<2>(c1),v7});
    auto p67=builder.Pipe(cs,cs1,v6,v7,2*scale,3*scale);
    x+=6;
    LOG::printf("X=%P\n",(x-6+2)*scale);
    LOG::printf("X=%P\n",(x-6+5)*scale);
    VERT_ID v8=builder.Vertex(TV(x,0)*scale);
    auto p78=builder.Pipe(cs1,cs,v7,v8,2*scale,3*scale);
    x+=1;
    LOG::printf("X=%P\n",x*scale);
    LOG::printf("Y=%P\n",-6*scale);
    LOG::printf("Y=%P\n",-8*scale);
    VERT_ID v9=builder.Vertex(TV(x,0)*scale);
    VERT_ID v10=builder.Vertex(TV(x,-6)*scale);
    VERT_ID v10a=builder.Vertex(TV(x,-8)*scale);
    auto src10a=builder.Set_BC(cs,v10a,v10,flowrate);
    x+=1;
    LOG::printf("X=%P\n",x*scale);
    T y=-3;
    LOG::printf("Y=%P\n",y*scale);
    VERT_ID v11=builder.Vertex(TV(x,y)*scale);
    auto j9=builder.Joint(cs,2,v9,{v8,v11});
    auto j10=builder.Joint(cs,2,v10,{v10a,v11});
    x+=3;
    LOG::printf("X=%P\n",x*scale);
    VERT_ID v12=builder.Vertex(TV(x,y)*scale);
    auto j11=builder.Joint(cs,3,v11,{v9,v10,v12});
    auto c2=Coil<VERT_ID,CID>(builder,cs,TV(x,y)*scale,4.5*scale,4.5*scale);
    auto j12=builder.Joint(cs,2,v12,{v11,std::get<0>(c2)});
    x+=4.5;
    LOG::printf("X=%P\n",x*scale);
    auto c3=Coil<VERT_ID,CID>(builder,cs,TV(x,y)*scale,4.5*scale,4.5*scale);
    x+=4.5;
    LOG::printf("X=%P\n",x*scale);
    VERT_ID v13=builder.Vertex(TV(x,y)*scale);
    x+=2;
    LOG::printf("X=%P\n",x*scale);
    VERT_ID v14=builder.Vertex(TV(x,y)*scale);
    auto j13=builder.Joint(cs,2,v13,{std::get<2>(c3),v14});
    auto sink=builder.Set_BC(cs,v14,v13,TV());

    LOG::printf("END simple\n");

    builder.Pipe(cs,src0.x,j1(0));
    builder.Pipe(cs,src2.x,j3(0));
    builder.Pipe(cs,j1(1),j4(0));
    builder.Pipe(cs,j3(1),j4(1));
    builder.Pipe(cs,j4(2),j5(0));
    builder.Pipe(cs,j5(1),std::get<1>(c0));
    builder.Pipe(cs,std::get<3>(c0),std::get<1>(c1));
    builder.Pipe(cs,std::get<3>(c1),j6(0));
    builder.Pipe(cs,j6(1),p67(0));
    builder.Pipe(cs1,p67(1),p78(0));
    builder.Pipe(cs,p78(1),j9(0));
    builder.Pipe(cs,j9(1),j11(0));
    builder.Pipe(cs,src10a.x,j10(0));
    builder.Pipe(cs,j10(1),j11(1));
    builder.Pipe(cs,j11(2),j12(0));
    builder.Pipe(cs,j12(1),std::get<1>(c2));
    builder.Pipe(cs,std::get<3>(c2),std::get<1>(c3));
    builder.Pipe(cs,std::get<3>(c3),j13(0));
    builder.Pipe(cs,j13(1),sink.x);

    FILE* file=fopen(LOG::sprintf("%s.txt",filename).c_str(),"w");
    fprintf(file,"%s\n",builder.To_String().c_str());
}

void Generate_Fab0(T mu,T s,T m,T kg)
{
    typedef VECTOR<T,2> TV;
    T dx=1.25;

    COMPONENT_LAYOUT_FEM<T> cl;
    cl.unit_m=m;
    cl.unit_s=s;
    cl.unit_kg=kg;
    LAYOUT_BUILDER_FEM<T> builder(cl);
    builder.Set_Target_Length(dx);
    builder.Set_Depth(dx,1);
    auto cs=builder.Cross_Section(4,5);

    using VERT_ID=LAYOUT_BUILDER_FEM<T>::VERT_ID;
    using CID=LAYOUT_BUILDER_FEM<T>::CONNECTOR_ID;

    VERT_ID v0=builder.Vertex(TV(70.353,422.351));
    VERT_ID v1=builder.Vertex(TV(70.353,352.552));
    VERT_ID v2=builder.Vertex(TV(323.853,352.552));
    VERT_ID v3=builder.Vertex(TV(323.853,306.858));
    VERT_ID v4=builder.Vertex(TV(221.353,306.858));
    VERT_ID v5=builder.Vertex(TV(221.353,90.351));
    VERT_ID v6=builder.Vertex(TV(323.853,331.352));
    TV O(336.351,331.352);
    VERT_ID v7=builder.Vertex(O);

    VERT_ID v9=builder.Vertex(TV(452.354,331.352));
    VERT_ID v10=builder.Vertex(TV(452.354,292.352));
    VERT_ID v11=builder.Vertex(TV(532.353,292.352));
    VERT_ID v12=builder.Vertex(TV(532.353,241.351));

    auto src0=builder.Set_BC(cs,v0,v1,1);
    auto src5=builder.Set_BC(cs,v5,v4,1);
    auto sink12=builder.Set_BC(cs,v12,v11,TV());

    auto j1=builder.Joint(cs,2,v1,{v0,v2});
    auto j2=builder.Joint(cs,2,v2,{v1,v6});
    auto j3=builder.Joint(cs,2,v3,{v6,v4});
    auto j4=builder.Joint(cs,2,v4,{v3,v5});
    auto j6=builder.Joint(cs,3,v6,{v2,v3,v7});
    auto j10=builder.Joint(cs,2,v10,{v9,v11});
    auto j11=builder.Joint(cs,2,v11,{v10,v12});

    T l=20.390,h=45.694;
    VERT_ID q=builder.Vertex(O+TV(0,h));
    auto j7=builder.Joint(cs,2,v7,{v6,q});
    CID last=j7(1);
    VERT_ID last_v=v7;
    for(int i=0;i<5;i++)
    {
        auto c=Coil<VERT_ID,CID>(builder,cs,O,l,h);
        builder.Pipe(cs,last,std::get<1>(c));
        last=std::get<3>(c);
        last_v=std::get<2>(c);
        O+=TV(l,0);
    }
    VERT_ID w=builder.Vertex(O);
    auto jw=builder.Joint(cs,2,w,{last_v,v9});
    auto j9=builder.Joint(cs,2,v9,{w,v10});
    builder.Pipe(cs,last,jw(0));
    builder.Pipe(cs,jw(1),j9(0));

    builder.Pipe(cs,src0.x,j1(0));
    builder.Pipe(cs,j1(1),j2(0));
    builder.Pipe(cs,j2(1),j6(0));
    builder.Pipe(cs,j6(2),j7(0));
    builder.Pipe(cs,j6(1),j3(0));
    builder.Pipe(cs,j3(1),j4(0));
    builder.Pipe(cs,j4(1),src5.x);
    builder.Pipe(cs,j9(1),j10(0));
    builder.Pipe(cs,j10(1),j11(0));
    builder.Pipe(cs,j11(1),sink12.x);

    printf("%s\n",builder.To_String().c_str());
}

namespace boost {
namespace polygon {
template<>
struct geometry_concept<VECTOR<T,2> >
{
    typedef point_concept type;
};

template<>
struct point_traits<VECTOR<T,2> >
{
    typedef int coordinate_type;

    static inline coordinate_type get(const VECTOR<T,2>& point,orientation_2d orient)
    {
        return (orient==HORIZONTAL)?point.x:point.y;
    }
};
}
}

void Generate_Voronoi_Pipes(T scale,RANDOM_NUMBERS<T>& rng,T mu,T s,T m,T kg,int n,T radius,const char* filename)
{
    T w=0.2*scale;
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    using VERT_ID=LAYOUT_BUILDER_FEM<T>::VERT_ID;
    using CID=LAYOUT_BUILDER_FEM<T>::CONNECTOR_ID;

    COMPONENT_LAYOUT_FEM<T> cl;
    cl.unit_m=m;
    cl.unit_s=s;
    cl.unit_kg=kg;
    LAYOUT_BUILDER_FEM<T> builder(cl);
    builder.Set_Target_Length(w/2);
    builder.Set_Depth(w,1);
    auto cs=builder.Cross_Section(2,w);

    SPHERE<TV> sphere(TV(),radius);
    GRID<TV> grid(TV_INT()+1,RANGE<TV>(TV()-radius,TV()+radius),true);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE(0.f),grid,LOG::sprintf("out-%s",filename));
    Flush_Frame<TV>("init");

    ARRAY<TV> X;
    POISSON_DISK<TV> poisson_disk(1);
    poisson_disk.Set_Distance_By_Volume(pi*sqr(radius)/n);
    ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> > obj(sphere);
    poisson_disk.Sample(rng,obj,X);

    for(const auto& x:X)
        Add_Debug_Particle(x,VECTOR<T,3>(0.5,0.5,0.5));
    Flush_Frame<TV>("samples");

    int N=100;
    for(int i=0;i<N;i++)
    {
        ROTATION<TV> Q=ROTATION<TV>::From_Angle(2*pi/N*i);
        Add_Debug_Particle(Q.Rotate(TV(1,0)*radius),VECTOR<T,3>(1,0,0));
    }

    using boost::polygon::voronoi_diagram;
    typedef boost::polygon::voronoi_vertex<T> VORONOI_VERTEX;
    voronoi_diagram<T> vd;
    construct_voronoi(X.begin(),X.end(),&vd);
    for(int j=0;j<X.m;j++)
    {
        Add_Debug_Particle(X(j),VECTOR<T,3>(0.5,0.5,0.5));
        Add_Debug_Text(X(j),LOG::sprintf("%d",j),VECTOR<T,3>(1,1,1));
    }

    HASHTABLE<const VORONOI_VERTEX*,VERT_ID> mapping;
    HASHTABLE<VERT_ID> velocity_bc;
    ARRAY<VECTOR<VERT_ID,2> > edge_ids;
    ARRAY<VECTOR<TV,2> > edges;
    auto append_clip=[scale,radius,&edges,&edge_ids,&builder,&mapping,&velocity_bc]
        (TV a,const VORONOI_VERTEX* pa,TV b,const VORONOI_VERTEX* pb)
    {
        T ra=a.Magnitude();
        if(ra>radius)
        {
            std::swap(a,b);
            std::swap(pa,pb);
            ra=a.Magnitude();
        }
        if(ra<=radius)
        {
            TV n=b-a;
            T r=n.Normalize();
            T t=-a.Dot(n);
            T dis=sqrt(sqr(a.Dot(n))-(a.Dot(a)-sqr(radius)));
            if(t+dis>=0) t+=dis;
            else t-=dis;

            auto pr=mapping.Insert(pa,VERT_ID(-1));
            if(pr.y)
                *pr.x=builder.Vertex(a*scale);
            VERT_ID vb(-1);

            if(t<r)
            {
                TV q=a+n*t;
                vb=builder.Vertex(q*scale);
                if(q.x<=0)
                    velocity_bc.Insert(vb);
            }
            else
            {
                t=r;
                auto pr=mapping.Insert(pb,VERT_ID(-1));
                if(pr.y)
                    *pr.x=builder.Vertex(b*scale);
                vb=*pr.x;
            }
            edges.Append(VECTOR<TV,2>(a,a+n*t));
            edge_ids.Append(VECTOR<VERT_ID,2>(*pr.x,vb));
        }
    };

    for(voronoi_diagram<T>::const_edge_iterator edge=vd.edges().begin();edge!=vd.edges().end();++edge)
    {
        if(edge->is_finite())
        {
            if(edge->cell()->source_index() < edge->twin()->cell()->source_index())
            {
                TV a(edge->vertex0()->x(),edge->vertex0()->y());
                TV b(edge->vertex1()->x(),edge->vertex1()->y());
                append_clip(a,edge->vertex0(),b,edge->vertex1());
            }
        }
        else if(edge->vertex0())
        {
            TV s(edge->vertex0()->x(),edge->vertex0()->y());
            TV p0=X(edge->cell()->source_index()),p1=X(edge->twin()->cell()->source_index());
            TV d=(p1-p0).Normalized().Rotate_Counterclockwise_90();
            append_clip(s,edge->vertex0(),s+d*radius,0);
        }
    }

    T min_l=radius*2;
    for(int i=0;i<edges.m;i++)
    {
        const auto& e=edges(i);
        Add_Debug_Object(e,VECTOR<T,3>(1,1,1));
        Add_Debug_Text((e(0)+e(1))*0.5,LOG::sprintf("%d",i),VECTOR<T,3>(1,0,1));
        min_l=min(min_l,(e(0)-e(1)).Magnitude());
    }
    LOG::printf("min length: %P\n",min_l);
    Flush_Frame<TV>("voronoi");


    {
        EPS_FILE<T> eps(LOG::sprintf("%s.eps",filename));
        eps.cur_format.line_width=1.0/radius;
        eps.cur_format.point_radius=2.0/radius;
        eps.cur_format.fill_style=1;
        for(int i=0;i<edges.m;i++)
        {
            eps.cur_format.line_style=1;
            const auto& e=edges(i);
            eps.Draw_Object(e(0),e(1));

            TV v0=e(0),v1=e(1);
            T r0=e(0).Magnitude(),r1=e(1).Magnitude();
            if(r0>r1)
            {
                std::swap(v0,v1);
                std::swap(r0,r1);
            }
            if(r1>0.95*radius){
                // BC
                TV dir=(v1-v0).Normalized();
                eps.cur_format.line_style=0;
                if(v1(0)<=0)
                    eps.cur_format.fill_color=VECTOR<T,3>(0,0,1);
                else
                    eps.cur_format.fill_color=VECTOR<T,3>(1,0,0);
                eps.Draw_Object(v1+dir*0.05*radius);
            }
        }

        eps.cur_format.point_radius=1.0/radius;
        eps.cur_format.fill_color=VECTOR<T,3>(0.5,0.5,0.5);
        for(int j=0;j<X.m;j++)
        {
            eps.Draw_Object(X(j));
            LOG::printf("V%d=(%P,%P)\n",j,X(j)(0)*scale,X(j)(1)*scale);
        }
    }


    ARRAY<ARRAY<CID,VERT_ID>,VERT_ID> links(builder.verts.m);
    for(auto& l:links)
        l.Resize(builder.verts.m,use_init,CID(-7));

    for(const auto& ei:edge_ids)
        links(ei(0))(ei(1))=links(ei(1))(ei(0))=CID(-1);

    for(VERT_ID i(0);i<links.m;i++)
    {
        ARRAY<VERT_ID> arms;
        for(VERT_ID j(0);j<links.m;j++)
            if(links(i)(j)==CID(-1)) arms.Append(j);
        PHYSBAM_ASSERT(arms.m==3 || arms.m==1);
        if(arms.m==3)
        {
            auto k=builder.Joint(cs,arms.m,i,arms);
            int l=0;
            for(VERT_ID j(0);j<links.m;j++)
                if(links(i)(j)==CID(-1)) links(i)(j)=k(l++);
        }
        else
        {
            CID bc(-1);
            if(velocity_bc.Contains(i))
                bc=builder.Set_BC(cs,i,arms(0),0.005).x;
            else
                bc=builder.Set_BC(cs,i,arms(0),TV()).x;

            for(VERT_ID j(0);j<links.m;j++)
                if(links(i)(j)==CID(-1)) links(i)(j)=bc;
        }
    }
    for(const auto& ei:edge_ids)
        builder.Pipe(cs,links(ei(0))(ei(1)),links(ei(1))(ei(0)));

    FILE* file=fopen(LOG::sprintf("%s.txt",filename).c_str(),"w");
    fprintf(file,"%s\n",builder.To_String().c_str());
}

template<int d>
void Run(PARSE_ARGS& parse_args)
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    
    std::string analytic_u,analytic_p,output_dir("output");
    bool gen_tests=false;
    T mu=1;
    T s=1,m=1,kg=1;
    int seed=time(0);
    parse_args.Add("-o",&output_dir,"dir","output dir");
    parse_args.Add("-mu",&mu,"mu","viscosity");
    parse_args.Add("-m",&m,"scale","scale units of length");
    parse_args.Add("-s",&s,"scale","scale units of time");
    parse_args.Add("-kg",&kg,"scale","scale units of mass");
    parse_args.Add("-seed",&seed,"seed","random seed");
    parse_args.Add("-gen",&gen_tests,"generate tests");
    parse_args.Add("-cache_size",&cache_size,"size","number of cache entries to use");
    parse_args.Add("-cache",&cache_pattern,"pattern","printf pattern for generating cache filenames");
    parse_args.Parse(true);
    if(!gen_tests)
    {
        parse_args.Extra(&analytic_u,"program","analytic velocity");
        parse_args.Extra(&analytic_p,"program","analytic pressure");
    }
    parse_args.Parse();
    LOG::printf("seed: %P\n",seed);
    RANDOM_NUMBERS<T> rng(seed);
    if(gen_tests)
    {
        Generate_Simple(1.0/40,mu,s,m,kg,"simple");
        Generate_Random_Ortho_Grid(1/0.2,grid0,sizeof(grid0)/sizeof(grid0[0]),mu,s,m,kg,"rgrid0");
        Generate_Random_Ortho_Grid(1/0.2,grid1,sizeof(grid1)/sizeof(grid1[0]),mu,s,m,kg,"rgrid1");
        Generate_Grid(1.0/80,mu,s,m,kg,20,"grid20");
        for(int i=3;i<=10;i++)
            Generate_Grid(1,mu,s,m,kg,i,LOG::sprintf("grid%d",i).c_str());
        {
            RANDOM_NUMBERS<T> rng(4);
            Generate_Voronoi_Pipes(1.0/20,rng,mu,s,m,kg,15,10.0,"voronoi-s4");
        }
        {
            RANDOM_NUMBERS<T> rng(15);
            Generate_Voronoi_Pipes(1.0/20,rng,mu,s,m,kg,15,10.0,"voronoi-s15");
        }
        return;
    }

    for(int i=0;i<10;i++)
    {
        Test_Pipe<d>(rng,mu,m,s,kg,analytic_u,analytic_p);
    }

    for(int i=0;i<10;i++)
    {
        Test_Var_Size_Pipe<d>(rng,mu,m,s,kg,analytic_u,analytic_p);
    }

    for(int i=0;i<5;i++)
    {
        T angle=rng.Get_Uniform_Number(0.2*pi,1.8*pi);
        Test_Joint2<d>(angle,mu,m,s,kg,analytic_u,analytic_p);
        Test_Joint2<d>(-angle,mu,m,s,kg,analytic_u,analytic_p);
    }

    for(int i=0;i<10;i++)
    {
        T a0=rng.Get_Uniform_Number(0.2*pi,0.9*pi),a1=rng.Get_Uniform_Number(0.2*pi,0.9*pi);
        Test_Joint3<d>(a0,a1,mu,m,s,kg,analytic_u,analytic_p);
    }

    for(int i=0;i<10;i++)
    {
        T angle=rng.Get_Uniform_Number(-2*pi,2*pi);
        Test_Joint4_Right_Angle<d>(angle,mu,m,s,kg,analytic_u,analytic_p);
    }
}

int main(int argc, char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"use 3D");
    parse_args.Parse(true);
    if(use_3d) Run<3>(parse_args);
    else Run<2>(parse_args);
    return 0;
}

