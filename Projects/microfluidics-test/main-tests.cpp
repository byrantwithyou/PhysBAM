//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/ROTATION.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <string>
#include "ANALYTIC_FEM.h"
#include "CACHED_ELIMINATION_MATRIX.h"
#include "COMPONENT_LAYOUT_FEM.h"
#include "ELIMINATION_FEM.h"
#include "LAYOUT_BUILDER_FEM.h"
#include "MATRIX_CONSTRUCTION_FEM.h"

using namespace PhysBAM;

typedef float RW;
typedef double T;

template<int d>
void Solve_And_Check(COMPONENT_LAYOUT_FEM<T>& cl,const LAYOUT_BUILDER_FEM<T>& builder,
    T mu,const std::string& au,const std::string& ap)
{
    typedef VECTOR<T,d> TV;
    cl.Update_Masters();
    cl.Merge_Blocks();
    MATRIX_CONSTRUCTION_FEM<TV> mc(cl);
    mc.mu=mu;
    CACHED_ELIMINATION_MATRIX<T> cem;
    cem.quiet=true;
    cl.Compute_Dof_Pairs();
    cl.Fill_Reference_Ticks();
    mc.Compute_Matrix_Blocks();
    ANALYTIC_FEM<TV> an(mc);
    an.Set_Velocity(au.c_str());
    an.Set_Pressure(ap.c_str());
    an.Compute_RHS();

    mc.Dump_World_Space_System();
    mc.Copy_To_CEM(cem);

    mc.Transform_Solution(cem,true,true);
    //mc.Dump_World_Space_Vector("b");

    ELIMINATION_FEM<T> el(cl);
    el.Eliminate_Irregular_Blocks(cem);
    el.Eliminate_Non_Seperators(cem);
    cem.Full_Reordered_Elimination();
    cem.Back_Solve();
    cem.Execute_Jobs(1);
    mc.Transform_Solution(cem,false,false);
    if(!an.Check_Analytic_Solution(false))
        printf("%s\n",builder.To_String().c_str());
    //mc.Dump_World_Space_Vector("x");
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
void Run(PARSE_ARGS& parse_args)
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    
    std::string analytic_u,analytic_p,output_dir("output");
    T mu=1;
    T s=1,m=1,kg=1;
    int seed=time(0);
    parse_args.Add("-o",&output_dir,"dir","output dir");
    parse_args.Add("-mu",&mu,"mu","viscosity");
    parse_args.Add("-m",&m,"scale","scale units of length");
    parse_args.Add("-s",&s,"scale","scale units of time");
    parse_args.Add("-kg",&kg,"scale","scale units of mass");
    parse_args.Add("-seed",&seed,"seed","random seed");
    parse_args.Extra(&analytic_u,"program","analytic velocity");
    parse_args.Extra(&analytic_p,"program","analytic pressure");
    parse_args.Parse();
    LOG::printf("seed: %P\n",seed);
    RANDOM_NUMBERS<T> rng(seed);

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

