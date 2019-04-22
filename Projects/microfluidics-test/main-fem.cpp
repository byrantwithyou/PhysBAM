//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/PAIR.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
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
#include "MATRIX_CONSTRUCTION_FEM.h"
#include <chrono>

using namespace PhysBAM;

typedef float RW;
typedef double T;

template<int d>
void Run(PARSE_ARGS& parse_args)
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,TV::m> TV_INT;

    ARRAY<PAIR<std::chrono::steady_clock::time_point,const char*> > tm;
    tm.Preallocate(32);
    auto timer=[&tm](const char* name){tm.Append({std::chrono::steady_clock::now(),name});};
    
    T mu=1;
    timer("start");

    int threads=1;
    bool quiet=false,use_krylov=false,print_system=false;
    std::string pipe_file,output_dir="output";
    std::string analytic_u,analytic_p;
    T s=1,m=1,kg=1;
    parse_args.Add("-o",&output_dir,"dir","output dir");
    parse_args.Add("-mu",&mu,"mu","viscosity");
    parse_args.Add("-q",&quiet,"disable diagnostics; useful for timing");
    parse_args.Add("-k",&use_krylov,"solve with Krylov method");
    parse_args.Add("-d",&print_system,"dump the system to be solved");
    parse_args.Add("-m",&m,"scale","scale units of length");
    parse_args.Add("-s",&s,"scale","scale units of time");
    parse_args.Add("-kg",&kg,"scale","scale units of mass");
    parse_args.Add("-threads",&threads,"num","number of threads to use");
    parse_args.Add("-u",&analytic_u,"program","analytic velocity");
    parse_args.Add("-p",&analytic_p,"program","analytic pressure");
    parse_args.Extra(&pipe_file,"file","file describing pipes");
    parse_args.Parse();

    timer("args");

    COMPONENT_LAYOUT_FEM<T> cl;
    // A
    cl.unit_m=m;
    cl.unit_s=s;
    cl.unit_kg=kg;
//    cl.mu=mu*kg/s;
    cl.Parse_Input(pipe_file);

    timer("parse input");

    GRID<TV> grid(TV_INT()+1,cl.Compute_Bounding_Box(),true);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE(0.f),grid,output_dir);
    vo.debug_particles.debug_particles.template Add_Array<T>("display_size");
    DEBUGGING_FEM<T> debug(cl);
    
    if(!quiet)
    {
        Flush_Frame<TV>("init");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            debug.Visualize_Block_State(b);
            debug.Visualize_Ticks(b,false);
        }
        Flush_Frame<TV>("blocks");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            debug.Visualize_Block_State(b);
            Flush_Frame<TV>(LOG::sprintf("block %P (%P)",b,cl.blocks(b).block).c_str());
        }
    }

    timer("setup viewing");

    cl.Update_Masters();

    if(!quiet)
    {
        Flush_Frame<TV>("after master");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            debug.Visualize_Block_State(b);
        Flush_Frame<TV>("blocks");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            debug.Visualize_Block_State(b);
            Flush_Frame<TV>(LOG::sprintf("block %P (%P)",b,cl.blocks(b).block).c_str());
        }
    }

    timer("masters");
    
    cl.Merge_Blocks();

    timer("merge blocks");

    if(!quiet)
    {
        Flush_Frame<TV>("after merge");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            debug.Visualize_Block_State(b);
        Flush_Frame<TV>("blocks");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            debug.Visualize_Block_State(b);
            Flush_Frame<TV>(LOG::sprintf("block %P (%P)",b,cl.blocks(b).block).c_str());
        }
    }

    MATRIX_CONSTRUCTION_FEM<TV> mc(cl);
    CACHED_ELIMINATION_MATRIX<T> cem;
    cem.quiet=quiet;
    cl.Compute_Dof_Pairs();
    cl.Fill_Reference_Ticks();
    if(!quiet)
    {
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            debug.Visualize_Block_State(b);
            debug.Visualize_Ticks(b,true);
        }
        Flush_Frame<TV>("ref ticks");
    }
    mc.Compute_Matrix_Blocks();
    ANALYTIC_FEM<TV>* an=0;
    if(analytic_p.size() && analytic_u.size())
    {
        an=new ANALYTIC_FEM<TV>(mc);
        an->Set_Velocity(analytic_u.c_str());
        an->Set_Pressure(analytic_p.c_str());
        an->Compute_RHS();
    }
    else
    {
        mc.Compute_RHS();
    }
    mc.Dump_World_Space_System();
    mc.Copy_To_CEM(cem);

    mc.Transform_Solution(cem,true,true);
    mc.Dump_World_Space_Vector("b");
    if(!quiet)
    {
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            debug.Visualize_Solution(mc.rhs_block_list(b),b,true);
        Flush_Frame<TV>("rhs blocks");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            debug.Visualize_Block_Dofs(b);
    }

    timer("compute matrix");

    ELIMINATION_FEM<T> el(cl);
    
    printf("==========================\n");
    cem.Print_Current();
    printf("==========================\n");
    el.Eliminate_Irregular_Blocks(cem);

    timer("elim irreg");

    cem.Print_Current();
    printf("==========================\n");
    el.Eliminate_Non_Seperators(cem);

    timer("elim non sep");

    cem.Print_Current();
    printf("==========================\n");
    cem.Full_Reordered_Elimination();
    cem.Print_Current();
    printf("==========================\n");

    timer("elim 3");

    cem.Back_Solve();

    timer("back solve");

    cem.Execute_Jobs(threads);

    timer("exec jobs");

    mc.Transform_Solution(cem,false,false);
    if(an) an->Check_Analytic_Solution();
    if(!quiet)
    {
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            debug.Visualize_Solution(mc.rhs_block_list(b),b,true);
        Flush_Frame<TV>("solution");
    }
    mc.Dump_World_Space_Vector("x");
    debug.Visualize_Flat_Dofs();

    for(int i=1;i<tm.m;i++)
        printf("%20s %5.0f ms\n",tm(i).y,
        std::chrono::duration_cast<std::chrono::duration<double> >(tm(i).x-tm(i-1).x).count()*1000);
}

int main(int argc, char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    bool use_3d=false;
    bool use_fem=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"use 3D");
    parse_args.Add("-fem",&use_fem,"use FEM");
    parse_args.Parse(true);
    LOG::Initialize_Logging(false,false,1<<30,true);
    /*if(use_3d) Run<3>(parse_args);
      else*/ Run<2>(parse_args);

    LOG::Finish_Logging();
    return 0;
}

