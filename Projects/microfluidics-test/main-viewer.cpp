//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/PAIR.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include "ANALYTIC_FEM.h"
#include "CACHED_ELIMINATION_MATRIX.h"
#include "COMPONENT_LAYOUT_FEM.h"
#include "DEBUGGING_FEM.h"
#include "ELIMINATION_FEM.h"
#include "LAYOUT_BUILDER_FEM.h"
#include "MATRIX_CONSTRUCTION_FEM.h"
#include "SOLUTION_FEM.h"
#include <chrono>


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

    std::string output_dir="sol_out",sol_file;
    TV X;
    T min_color=1e-8,max_color=20;
    bool dump_grad=false;
    parse_args.Add("-o",&output_dir,"dir","output dir");
    parse_args.Add("-g",&dump_grad,"dump gradients");
    parse_args.Add("-min_color",&min_color,"value","minimum value for the coolest color");
    parse_args.Add("-max_color",&max_color,"value","maximum value for the hottest color");
    parse_args.Add("-x",&X,"position","interpolate solution");
    parse_args.Extra(&sol_file,"file","solution file");
    parse_args.Parse();

    GRID<TV2> grid(IV2()+1,RANGE<TV2>::Centered_Box(),true);
    VIEWER_OUTPUT<TV2> vo2(STREAM_TYPE(0.f),grid,output_dir);
    GRID<TV3> grid3(IV3()+1,RANGE<TV3>::Centered_Box(),true);
    VIEWER_OUTPUT<TV3> vo3(STREAM_TYPE(0.f),grid3,output_dir+"/3d");
    vo2.debug_particles.debug_particles.template Add_Array<T>("display_size");
    
    SOLUTION_FEM<TV> sol;
    Read_From_File(sol_file,sol);
    sol.Prepare_Hierarchy();

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

    Flush_Frame<TV>("init");

    dump_mesh();
    TV vel=sol.Velocity(X);
    T pres=sol.Pressure(X);
    LOG::printf("solution at %P: %P %P\n",X,vel,pres);
    Add_Debug_Particle(X,VECTOR<T,3>(1,0,0));
    Debug_Particle_Set_Attribute<TV>("V",vel);
    Debug_Particle_Set_Attribute<TV>("display_size",pres);
    Flush_Frame<TV>("solution");

    if(dump_grad)
    {
        INTERPOLATED_COLOR_MAP<T> icm;
        icm.Initialize_Colors(min_color,max_color,false,true,false);
        T min_f=1e12,max_f=-1e12;
        for(const auto& e:sol.mesh.elements)
        {
            TV x=sol.particles.X.Subset(e).Sum()/e.m;
            TM g;
            sol.Velocity(x,&g);
            T n=g.Frobenius_Norm();
            min_f=min(min_f,n);
            max_f=max(max_f,n);
            Add_Debug_Particle(x,icm(n));
        }
        LOG::printf("grad V norm: %P %P\n",min_f,max_f);
        Flush_Frame<TV>("grad V");

        icm.Initialize_Colors(min_color,max_color,false,true,false);
        min_f=1e12;
        max_f=-1e12;
        for(const auto& e:sol.mesh.elements)
        {
            TV x=sol.particles.X.Subset(e).Sum()/e.m;
            TV g;
            sol.Pressure(x,&g);
            T n=g.Magnitude();
            min_f=min(min_f,n);
            max_f=max(max_f,n);
            Add_Debug_Particle(x,icm(n));
        }
        LOG::printf("grad P norm: %P %P\n",min_f,max_f);
        Flush_Frame<TV>("grad P");
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
