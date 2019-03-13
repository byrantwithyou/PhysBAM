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
#include "CACHED_ELIMINATION_MATRIX.h"
#include "COMPONENT_LAYOUT_FEM.h"
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
    parse_args.Add("-o",&output_dir,"dir","output dir");
    parse_args.Add("-mu",&mu,"mu","viscosity");
    parse_args.Add("-q",&quiet,"disable diagnostics; useful for timing");
    parse_args.Add("-k",&use_krylov,"solve with Krylov method");
    parse_args.Add("-p",&print_system,"dump the system to be solved");
    parse_args.Add("-threads",&threads,"num","number of threads to use");
    parse_args.Extra(&pipe_file,"file","file describing pipes");
    parse_args.Parse();

    timer("args");

    COMPONENT_LAYOUT_FEM<TV> cl;
    cl.mu=mu;
    cl.Parse_Input(pipe_file);

    timer("parse input");

    GRID<TV> grid(TV_INT()+1,cl.Compute_Bounding_Box(),true);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE(0.f),grid,output_dir);
    vo.debug_particles.debug_particles.template Add_Array<T>("display_size");

    if(!quiet)
    {
        Flush_Frame<TV>("init");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            cl.Visualize_Block_State(b);
        Flush_Frame<TV>("blocks");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            cl.Visualize_Block_State(b);
            Flush_Frame<TV>(LOG::sprintf("block %P (%P)",b,cl.blocks(b).block).c_str());
        }
    }

    timer("setup viewing");

    cl.Update_Masters();

    if(!quiet)
    {
        Flush_Frame<TV>("init");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            cl.Visualize_Block_State(b);
        Flush_Frame<TV>("blocks");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            cl.Visualize_Block_State(b);
            Flush_Frame<TV>(LOG::sprintf("block %P (%P)",b,cl.blocks(b).block).c_str());
        }
    }

    timer("masters");
    
    cl.Compute();

    timer("compute");

    CACHED_ELIMINATION_MATRIX<T> cem;
    cem.quiet=quiet;
    cl.Compute_Matrix_Blocks(cem);

    timer("compute matrix");

    printf("==========================\n");
    cem.Print_Current();
    printf("==========================\n");
    cl.Eliminate_Irregular_Blocks(cem);

    timer("elim irreg");

    cem.Print_Current();
    printf("==========================\n");
    cl.Eliminate_Non_Seperators(cem);

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

    cl.Transform_Solution(cem);
    if(!quiet)
    {
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            cl.Visualize_Solution(b);
        Flush_Frame<TV>("blocks");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            cl.Visualize_Solution(b);
            Flush_Frame<TV>(LOG::sprintf("block %P (%P)",b,cl.blocks(b).block).c_str());
        }
    }

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

