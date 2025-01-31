//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/PAIR.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Images/PNG_FILE.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include "ANALYTIC_FEM.h"
#include "CACHED_ELIMINATION_MATRIX.h"
#include "COMPONENT_LAYOUT_FEM.h"
#include "DEBUGGING_FEM.h"
#include "ELIMINATION_FEM.h"
#include "LAYOUT_BUILDER_FEM.h"
#include "MATRIX_CONSTRUCTION_FEM.h"
#include "SOLUTION_FEM.h"
#include <chrono>
#include <fstream>
#include <map>
#include <set>
#include <string>


typedef float RW;
typedef double T;

using namespace PhysBAM;

namespace PhysBAM{
extern bool use_job_timing;
}

template<int d>
void Run(PARSE_ARGS& parse_args)
{
    typedef VECTOR<T,d> TV;
    typedef MATRIX<T,d> TM;
    typedef VECTOR<T,3> TV3;
    typedef VECTOR<T,2> TV2;
    typedef VECTOR<int,2> IV2;
    typedef VECTOR<int,3> IV3;
    typedef VECTOR<T,4> RGBA;

    std::string output_dir="sol_out",sol_file,ref_sol_file;
    TV X;
    T min_color=1e-8,max_color=20;
    bool user_query=false,user_color_scale=false;
    int resolution=256;
    parse_args.Add("-o",&output_dir,"dir","output dir");
    parse_args.Add("-ref",&ref_sol_file,"file","reference solution file");
    parse_args.Add("-res",&resolution,"resolution","resolution");
    parse_args.Add("-min_color",&min_color,&user_color_scale,"value","minimum value for the coolest color");
    parse_args.Add("-max_color",&max_color,&user_color_scale,"value","maximum value for the hottest color");
    parse_args.Add("-x",&X,&user_query,"position","interpolate solution");
    parse_args.Extra(&sol_file,"file","solution file");
    parse_args.Parse();

    VIEWER_DIR viewer_dir2(output_dir);
    VIEWER_OUTPUT vo2(STREAM_TYPE(0.f),viewer_dir2);
    vo2.Use_Debug_Particles<TV2>();
    VIEWER_DIR viewer_dir3(output_dir+"/3d");
    VIEWER_OUTPUT vo3(STREAM_TYPE(0.f),viewer_dir3);
    vo3.Use_Debug_Particles<TV3>();
    
    SOLUTION_FEM<TV> sol;
    Read_From_File(sol_file,sol);
    sol.Prepare_Hierarchy();

    SOLUTION_FEM<TV> ref_sol;
    if(ref_sol_file!="")
    {
        Read_From_File(ref_sol_file,ref_sol);
        ref_sol.Prepare_Hierarchy();
    }

    auto dump_mesh=[&sol]()
    {
        for(const auto& x:sol.particles.X)
            Add_Debug_Particle(x,VECTOR<T,3>(1,1,1));
        for(const auto& t:sol.mesh.elements)
        {
            for(int i=0;i<t.m;i++)
                for(int j=i;j<t.m;j++)
                    Add_Debug_Object(VECTOR<TV,2>(sol.particles.X(t(i)),sol.particles.X(t(j))),VECTOR<T,3>(0.5,0.5,0.5));
        }
    };

    Flush_Frame("init");

    if(user_query)
    {
        dump_mesh();
        TV vel=sol.Velocity(X);
        T pres=sol.Pressure(X);
        LOG::printf("solution at %P: %P %P\n",X,vel,pres);
        Add_Debug_Particle(X,VECTOR<T,3>(1,0,0));
        Debug_Particle_Set_Attribute<TV>("V",vel);
        Debug_Particle_Set_Attribute<TV>("display_size",pres);
        Flush_Frame("solution");
    }

    INTERVAL<T> range_v=INTERVAL<T>::Empty_Box(),
        range_dv=INTERVAL<T>::Empty_Box(),
        range_p=INTERVAL<T>::Empty_Box(),
        range_dp=INTERVAL<T>::Empty_Box(),
        range_p_err=INTERVAL<T>::Empty_Box(),
        range_v_err=INTERVAL<T>::Empty_Box();
    RANGE<TV> box=RANGE<TV>::Bounding_Box(sol.particles.X);
    for(const auto& e:sol.mesh.elements)
    {
        TV x=sol.particles.X.Subset(e).Sum()/e.m;
        int elem=sol.Intersect(x);
        TM dv;
        TV dp;
        TV v=sol.Velocity(elem,x,&dv);
        T p=sol.Pressure(elem,x,&dp);
        range_v.Enlarge_To_Include_Point(v.Magnitude());
        range_dv.Enlarge_To_Include_Point(dv.Frobenius_Norm());
        range_p.Enlarge_To_Include_Point(p);
        range_dp.Enlarge_To_Include_Point(dp.Magnitude());
        if(ref_sol_file!="")
        {
            int ref_elem=ref_sol.Intersect(x);
            TV ref_v=ref_sol.Velocity(ref_elem,x);
            T ref_p=ref_sol.Pressure(ref_elem,x);
            range_v_err.Enlarge_To_Include_Point((v-ref_v).Magnitude());
            range_p_err.Enlarge_To_Include_Point(abs(p-ref_p));
        }
    }
    LOG::printf("v range: %P\n",range_v);
    LOG::printf("dv range: %P\n",range_dv);
    LOG::printf("p range: %P\n",range_p);
    LOG::printf("dp range: %P\n",range_dp);
    if(ref_sol_file!="")
    {
        LOG::printf("v error range: %P\n",range_v_err);
        LOG::printf("p error range: %P\n",range_p_err);
    }

    T dx=box.Edge_Lengths().Max()/resolution;
    int m=box.Edge_Lengths()(0)/dx,n=box.Edge_Lengths()(1)/dx;
    INTERPOLATED_COLOR_MAP<T> log10_icm;
    log10_icm.Initialize_Colors(1e-6,1,true,true,false);

    INTERPOLATED_COLOR_MAP<T> linear_icm;
    linear_icm.Initialize_Colors(0,1,false,true,false);

    // H: in degrees [0,360]
    auto HSV_To_RGB=[](T h,T s,T v)
    {
        T r=fmod(5+h/60,6),g=fmod(3+h/60,6),b=fmod(1+h/60,6);
        r=v-v*s*max(min(min(r,4-r),1.),0.);
        g=v-v*s*max(min(min(g,4-g),1.),0.);
        b=v-v*s*max(min(min(b,4-b),1.),0.);
        return VECTOR<T,3>(r,g,b);
    };

    auto normalize=[user_color_scale,min_color,max_color](T x,T mn,T mx)
    {
        if(user_color_scale)
        {
            mn=min_color;
            mx=max_color;
        }
        return (x-mn)/(mx-mn);
    };

    ARRAY<RGBA,IV2> img_v(IV2(m,n)),img_v_mag(IV2(m,n)),img_p(IV2(m,n)),img_dv(IV2(m,n)),img_dp(IV2(m,n));
    ARRAY<RGBA,IV2> img_v_err(IV2(m,n)),img_p_err(IV2(m,n));
    for(int j=0;TV::m==2 && j<n;j++) for(int i=0;i<m;i++)
    {
        TV x=box.Minimum_Corner()+TV::Axis_Vector(0)*i*dx+TV::Axis_Vector(1)*j*dx;
        int elem=sol.Intersect(x);
        if(elem<0) continue;
        TV dp;
        TM dv;
        T p=sol.Pressure(elem,x,&dp);
        TV v=sol.Velocity(elem,x,&dv);
        IV2 index(i,j);
        T angle=TV2::Oriented_Angle_Between(TV2::Axis_Vector(0),TV2(v(0),v(1)));
        if(angle<0) angle+=2*pi;
        angle=angle/pi*180;
        img_v(index)=HSV_To_RGB(angle,(v.Magnitude()-range_v.min_corner)/range_v.Size(),1).Append(1);
        img_v_mag(index)=linear_icm(normalize(v.Magnitude(),range_v.min_corner,range_v.max_corner)).Append(1);
        img_p(index)=linear_icm(normalize(p,range_p.min_corner,range_p.max_corner)).Append(1);
        img_dp(index)=log10_icm(normalize(dp.Magnitude(),range_dp.min_corner,range_dp.max_corner)).Append(1);
        img_dv(index)=log10_icm(normalize(dv.Frobenius_Norm(),range_dv.min_corner,range_dv.max_corner)).Append(1);

        if(ref_sol_file!="")
        {
            int ref_elem=ref_sol.Intersect(x);
            TV ref_v=ref_sol.Velocity(ref_elem,x);
            T ref_p=ref_sol.Pressure(ref_elem,x);
            img_v_err(index)=log10_icm(normalize((v-ref_v).Magnitude(),range_v.min_corner,range_v.max_corner)).Append(1);
            img_p_err(index)=log10_icm(normalize(abs(p-ref_p),1e-15,range_p.max_corner-range_p.min_corner)).Append(1);
        }
    }

    ARRAY<TV3,IV2> bar(IV2(1000,1));
    for(int i=0;i<1000;i++)
        bar(IV2(i,0))=linear_icm.colors.Value((i/(T)999)*(linear_icm.mx-linear_icm.mn)+linear_icm.mn);
    PNG_FILE<T>::Write(output_dir+"/bar.png",bar);
    for(int i=0;i<1000;i++)
    {
        T x=i/(T)999;
        T y=pow(10,-(1-x)*8);
        bar(IV2(i,0))=log10_icm.colors.Value(y*(log10_icm.mx-log10_icm.mn)+log10_icm.mn);
    }
    PNG_FILE<T>::Write(output_dir+"/log-bar.png",bar);

    ARRAY<RGBA,IV2> wheel(IV2(1000,1000));
    for(int i=0;i<1000;i++) for(int j=0;j<1000;j++)
    {
        TV2 dir(2*i/1000.0-1.0,2*j/1000.0-1.0);
        T s=dir.Magnitude();
        if(s>1) continue;
        T angle=TV2::Oriented_Angle_Between(TV2::Axis_Vector(0),dir);
        if(angle<0) angle+=2*pi;
        angle=angle/pi*180;
        wheel(IV2(i,j))=HSV_To_RGB(angle,s,1).Append(1);
    }
    PNG_FILE<T>::Write(output_dir+"/wheel.png",wheel);

    PNG_FILE<T>::Write(output_dir+"/v.png",img_v);
    PNG_FILE<T>::Write(output_dir+"/v-mag.png",img_v_mag);
    PNG_FILE<T>::Write(output_dir+"/p.png",img_p);
    PNG_FILE<T>::Write(output_dir+"/dp.png",img_dp);
    PNG_FILE<T>::Write(output_dir+"/dv.png",img_dv);
    if(ref_sol_file!="")
    {
        PNG_FILE<T>::Write(output_dir+"/err-p.png",img_p_err);
        PNG_FILE<T>::Write(output_dir+"/err-v.png",img_v_err);
    }
}

int main(int argc, char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"use 3D");
    parse_args.Parse(true);
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::printf("%s\n",parse_args.Print_Arguments());
    if(use_3d) Run<3>(parse_args);
    else Run<2>(parse_args);

    LOG::Finish_Logging();
    return 0;
}
