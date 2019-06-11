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
    typedef VECTOR<T,3> TV3;
    typedef VECTOR<T,2> TV2;
    typedef VECTOR<int,2> IV2;
    typedef VECTOR<int,3> IV3;

    ARRAY<PAIR<std::chrono::steady_clock::time_point,const char*> > tm;
    tm.Preallocate(32);
    auto timer=[&tm](const char* name){tm.Append({std::chrono::steady_clock::now(),name});};
    
    T mu=1;
    timer("start");

    int threads=1;
    bool quiet=false,use_krylov=false,print_system=false,pinv=false,stats_only=false,force_blk_ref=false,illus_domain=false,
        illus_meshing=false,dump_solution=false,rounded_corner=false;
    TV2 min_corner,max_corner;
    std::string pipe_file,output_dir="output";
    std::string analytic_u,analytic_p;
    std::string sol_file;
    T s=1,m=1,kg=1;
    int hl_dof=-1;
    int refine=1;
    parse_args.Add("-o",&output_dir,"dir","output dir");
    parse_args.Add("-check_sol",&sol_file,"file","compare with saved solution");
    parse_args.Add("-mu",&mu,"mu","viscosity");
    parse_args.Add("-q",&quiet,"disable diagnostics; useful for timing");
    parse_args.Add("-dump_sol",&dump_solution,"dump solution");
    parse_args.Add("-pinv",&pinv,"perform pseudo inverse");
    parse_args.Add("-force_blk_ref",&force_blk_ref,"force every block a reference block");
    parse_args.Add("-stats",&stats_only,"show statistics");
    parse_args.Add("-illus",&illus_domain,"dump domain");
    parse_args.Add("-illus_meshing",&illus_meshing,"dump meshing");
    parse_args.Add("-min_corner",&min_corner,"point","min corner");
    parse_args.Add("-max_corner",&max_corner,"point","max corner");
    parse_args.Add("-rounded_corner",&rounded_corner,"use rounded corner");
    parse_args.Add("-k",&use_krylov,"solve with Krylov method");
    parse_args.Add("-d",&print_system,"dump the system to be solved");
    parse_args.Add("-m",&m,"scale","scale units of length");
    parse_args.Add("-s",&s,"scale","scale units of time");
    parse_args.Add("-kg",&kg,"scale","scale units of mass");
    parse_args.Add("-threads",&threads,"num","number of threads to use");
    parse_args.Add("-hl",&hl_dof,"dof","highlight global dof");
    parse_args.Add("-refine",&refine,"factor","refine factor");
    parse_args.Add("-u",&analytic_u,"program","analytic velocity");
    parse_args.Add("-p",&analytic_p,"program","analytic pressure");
    parse_args.Add("-timing",&use_job_timing,"use job timing");
    parse_args.Extra(&pipe_file,"file","file describing pipes");
    parse_args.Parse();

    timer("args");

    COMPONENT_LAYOUT_FEM<T> cl;
    cl.unit_m=m;
    cl.unit_s=s;
    cl.unit_kg=kg;
    cl.force_blk_ref=force_blk_ref;
    LAYOUT_BUILDER_FEM<T> builder(cl);
    builder.refine=refine;
    builder.comp_joint.rounded_corner=rounded_corner;
    builder.From_File(pipe_file);

    timer("parse input");

    GRID<TV2> grid(IV2()+1,cl.Compute_Bounding_Box(),true);
    VIEWER_OUTPUT<TV2> vo2(STREAM_TYPE(0.f),grid,output_dir);
    GRID<TV3> grid3(IV3()+1,{grid.domain.min_corner.Append(0),grid.domain.max_corner.Append(cl.depth)},true);
    VIEWER_OUTPUT<TV3> vo3(STREAM_TYPE(0.f),grid3,output_dir+"/3d");
    vo2.debug_particles.debug_particles.template Add_Array<T>("display_size");
    // (local,global)
    vo3.debug_particles.debug_particles.template Add_Array<VECTOR<int,2>>("e");
    vo3.debug_particles.debug_particles.template Add_Array<VECTOR<int,2>>("v");
    vo3.debug_particles.debug_particles.template Add_Array<VECTOR<int,2>>("p");
    DEBUGGING_FEM<T> debug(cl);
    
    if(!quiet)
    {
        Flush_Frame<TV2>("init");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            debug.Visualize_Block_State(b);
            debug.Visualize_Ticks(b,false);
        }
        Flush_Frame<TV2>("blocks");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            debug.Visualize_Block_State(b);
            Flush_Frame<TV2>(LOG::sprintf("block %P (%P)",b,cl.blocks(b).block).c_str());
        }
    }

    timer("setup viewing");

    cl.Update_Masters();

    if(!quiet)
    {
        Flush_Frame<TV2>("after master");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            debug.Visualize_Block_State(b);
        Flush_Frame<TV2>("blocks");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            debug.Visualize_Block_State(b);
            Flush_Frame<TV2>(LOG::sprintf("block %P (%P)",b,cl.blocks(b).block).c_str());
        }
    }

    timer("masters");
    
    cl.Merge_Blocks();

    timer("merge blocks");

    if(!quiet)
    {
        Flush_Frame<TV2>("after merge");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            debug.Visualize_Block_State(b);
        Flush_Frame<TV2>("blocks");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            debug.Visualize_Block_State(b);
            Flush_Frame<TV2>(LOG::sprintf("block %P (%P)",b,cl.blocks(b).block).c_str());
        }
    }

    MATRIX_CONSTRUCTION_FEM<TV> mc(cl);
    mc.mu=mu*cl.unit_kg*pow<2-TV::m>(cl.unit_m)/cl.unit_s;
    CACHED_ELIMINATION_MATRIX<T> cem;
    cem.quiet=quiet;
    cem.pinv_last_blk=pinv;
    cl.Compute_Dof_Pairs();
    cl.Fill_Reference_Ticks();
    if(!quiet)
    {
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            debug.Visualize_Block_State(b);
            debug.Visualize_Ticks(b,true);
            if(d==3)
            {
                debug.Visualize_Tetrahedron(b);
                Flush_Frame<TV>(LOG::sprintf("block %P",b).c_str());
            }
        }
        Flush_Frame<TV2>("ref ticks");

        if(d==3)
            debug.Visualize_Tetrahedron_Dofs(mc);
    }
    if(hl_dof>=0)
    {
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            BLOCK_ID bb;
            int vep,r,dim;
            std::tie(bb,vep,r,dim)=mc.Inverse_DOF_Lookup(hl_dof);
            if(bb!=b) continue;

            debug.Visualize_Block_State(b);
            if(d==3)
                debug.Visualize_Tetrahedron(b);
            debug.Highlight_Dof<d>(b,vep,r,dim);
            Flush_Frame<TV>("highlight dof");
        }
    }
    mc.Print_Statistics();
    LOG::printf("canonical-j2: %d\ncanonical-j3-avg: %d\ncanonical-j3-small: %d\ncanonical-j4: %d\n",
        builder.comp_joint.num_j2,builder.comp_joint.num_j3_avg,builder.comp_joint.num_j3_small,builder.comp_joint.num_j4);
    if(illus_domain)
        debug.Visualize_Domain("domain.eps",RANGE<TV2>::Empty_Box());
    if(illus_meshing)
    {
        RANGE<TV2> range(min_corner,max_corner);
        debug.Visualize_Domain("domain-anno.eps",range);
        debug.Visualize_Meshing("meshing.eps",range);
    }

    if(stats_only || illus_domain || illus_meshing) return;

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
    if(!quiet) mc.Dump_World_Space_System();
    mc.Copy_To_CEM(cem);

    if(!quiet)
    {
        mc.Transform_Solution(cem,true,true);
        mc.Dump_World_Space_Vector("b");
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
    if(an) an->Check_Analytic_Solution(!quiet);
    if(!sol_file.empty())
    {
        SOLUTION_FEM<TV> sol;
        Read_From_File(sol_file,sol);
        sol.Prepare_Hierarchy();
        auto fv=[&sol](const TV& X){return sol.Velocity(X);};
        auto fp=[&sol](const TV& X){return sol.Pressure(X);};
        Check_Solution<T,TV::m>(mc,fv,fp,!quiet);
    }

    if(!quiet)
    {
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            debug.Visualize_Solution(mc.rhs_block_list(b),b,true);
        Flush_Frame<TV>("solution");
        mc.Dump_World_Space_Vector("x");
        debug.Visualize_Flat_Dofs();
    }

    for(int i=1;i<tm.m;i++)
        LOG::printf("%20s %5.0f ms\n",tm(i).y,
        std::chrono::duration_cast<std::chrono::duration<double> >(tm(i).x-tm(i-1).x).count()*1000);

    if(quiet)
    {
        Create_Directory(output_dir);
        Create_Directory(output_dir+"/common");
        LOG::Instance()->Copy_Log_To_File(output_dir+"/common/log.txt",false);
    }

    if(dump_solution)
    {
        SOLUTION_FEM<TV> solution;
        solution.Build(mc);
        Write_To_File(STREAM_TYPE((double)0),output_dir+"/sol",solution);
    }
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
    LOG::printf("%s\n",parse_args.Print_Arguments());
    if(use_3d) Run<3>(parse_args);
    else Run<2>(parse_args);

    LOG::Finish_Logging();
    return 0;
}
