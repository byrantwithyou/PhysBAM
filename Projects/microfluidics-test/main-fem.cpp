//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/PAIR.h>
#include <Core/Math_Tools/pow.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <Tools/Krylov_Solvers/MINRES.h>
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
#include "KRYLOV_SOLVER_FEM.h"
#include "LAYOUT_BUILDER_FEM.h"
#include "MATRIX_CONSTRUCTION_FEM.h"
#include "SOLUTION_FEM.h"
#include <chrono>
#if USE_MKL
#include <mkl.h>
#endif
#ifdef USE_OPENMP
#include <omp.h>
#endif

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
    int kn=100000;
    bool quiet=false,use_krylov=false,print_system=false,pinv=false,stats_only=false,force_blk_ref=false,illus_domain=false,
        illus_zoom=false,dump_solution=false,rounded_corner=false;
    bool visualize_blocks=true;
    TV2 min_corner,max_corner;
    IV2 illus_size(128,128);
    bool illus_fill=false,use_mkl_sparse=false,dump_sysbin=false;
    std::string pipe_file,output_dir="output";
    std::string analytic_u,analytic_p;
    std::string sol_file;
    T s=1,m=1,kg=1;
    int hl_dof=-1;
    int refine=1;
    int cache_size=-1;
    std::string cache_pattern="cache-%d.txt";
    parse_args.Add("-o",&output_dir,"dir","output dir");
    parse_args.Add("-check_sol",&sol_file,"file","compare with saved solution");
    parse_args.Add("-mu",&mu,"mu","viscosity");
    parse_args.Add("-q",&quiet,"disable diagnostics; useful for timing");
    parse_args.Add("-dump_sol",&dump_solution,"dump solution");
    parse_args.Add("-dump_sysbin",&dump_sysbin,"dump system in binary format");
    parse_args.Add("-pinv",&pinv,"perform pseudo inverse");
    parse_args.Add("-force_blk_ref",&force_blk_ref,"force every block a reference block");
    parse_args.Add("-stats",&stats_only,"show statistics");
    parse_args.Add("-illus",&illus_domain,"dump domain");
    parse_args.Add("-illus_size",&illus_size,"dimension","illustration size");
    parse_args.Add("-illus_fill",&illus_fill,"use filled illustration");
    parse_args.Add("-illus_zoom",&illus_zoom,"magnify part");
    parse_args.Add("-min_corner",&min_corner,"point","min corner");
    parse_args.Add("-max_corner",&max_corner,"point","max corner");
    parse_args.Add("-rounded_corner",&rounded_corner,"use rounded corner");
    parse_args.Add("-k",&use_krylov,"solve with Krylov method");
    parse_args.Add("-kn",&kn,"number","max krylov iterations");
    parse_args.Add("-mklsp",&use_mkl_sparse,"use MKL sparse matrix-vector multiply");
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
    parse_args.Add("-cache_size",&cache_size,"size","number of cache entries to use");
    parse_args.Add("-cache",&cache_pattern,"pattern","printf pattern for generating cache filenames");
    parse_args.Add_Not("-no_blocks",&visualize_blocks,"visualize per-black data");
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

    RANGE<TV2> box=cl.Compute_Bounding_Box();
    LOG::printf("bounding box: %P\n",box);
    GRID<TV2> grid(IV2()+1,box,true);
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
        if(visualize_blocks)
        {
            for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            {
                debug.Visualize_Block_State(b);
                Flush_Frame<TV2>(LOG::sprintf("block %P (%P)",b,cl.blocks(b).block).c_str());
            }
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
        if(visualize_blocks)
        {
            for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            {
                debug.Visualize_Block_State(b);
                Flush_Frame<TV2>(LOG::sprintf("block %P (%P)",b,cl.blocks(b).block).c_str());
            }
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
        if(visualize_blocks)
        {
            for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            {
                debug.Visualize_Block_State(b);
                Flush_Frame<TV2>(LOG::sprintf("block %P (%P)",b,cl.blocks(b).block).c_str());
            }
        }
    }

    CACHED_ELIMINATION_MATRIX<T> cem;
    cem.quiet=quiet;
    cem.pinv_last_blk=pinv;
    cl.Compute_Dof_Pairs();
    cl.Fill_Reference_Ticks();
    MATRIX_CONSTRUCTION_FEM<TV> mc(cl,cache_pattern+"-canonical",cache_size);
    mc.mu=mu*cl.unit_kg*pow<2-TV::m>(cl.unit_m)/cl.unit_s;
    if(!quiet)
    {
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
        {
            debug.Visualize_Block_State(b);
            debug.Visualize_Ticks(b,true);
            if(d==3)
            {
                debug.Visualize_Tetrahedron(b);
                if(visualize_blocks)
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

    if(illus_domain || illus_zoom)
        if(!Directory_Exists(output_dir))
            Create_Directory(output_dir);
    if(illus_domain)
        debug.Visualize_Domain(output_dir+"/domain.eps",illus_fill,illus_size,RANGE<TV2>::Empty_Box(),RANGE<TV2>::Empty_Box());
    if(illus_zoom)
    {
        RANGE<TV2> range(min_corner,max_corner);
        debug.Visualize_Domain(output_dir+"/domain-anno.eps",illus_fill,illus_size,range,RANGE<TV2>::Empty_Box());
        debug.Visualize_Domain(output_dir+"/domain-part.eps",false,illus_size,RANGE<TV2>::Empty_Box(),range);
        debug.Visualize_Meshing(output_dir+"/meshing.eps",illus_size,range);
    }

    if(stats_only || illus_domain || illus_zoom) return;

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

    timer("compute matrix");

    if(!quiet || use_krylov || dump_sysbin)
        mc.Transform_Solution(cem,true,true);
    if(!quiet)
    {
        mc.Dump_World_Space_Vector("b");
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            debug.Visualize_Solution(mc.rhs_block_list(b),b,true);
        Flush_Frame<TV>("rhs blocks");
        if(visualize_blocks)
            for(BLOCK_ID b(0);b<cl.blocks.m;b++)
                debug.Visualize_Block_Dofs(b);
    }

    timer("transform matrix");

    if(quiet)
    {
        Create_Directory(output_dir);
        Create_Directory(output_dir+"/common");
        LOG::Instance()->Copy_Log_To_File(output_dir+"/common/log.txt",false);
    }

    ARRAY<T> krylov_sol;
    if(use_krylov || dump_sysbin)
    {
        ARRAY<int,BLOCK_ID> first[3];
        int size=mc.Compute_Global_Dof_Mapping(first);
        ARRAY<T> b;
        mc.Dump_World_Space_Vector(first,size,b);
        SPARSE_MATRIX_FLAT_MXN<T> SM;
        mc.Dump_World_Space_System(first,size,SM);

        if(use_krylov)
        {
            typedef KRYLOV_VECTOR_FEM<T> KRY_VEC;
            typedef KRYLOV_SYSTEM_FEM<T> KRY_MAT;
            KRY_MAT sys(SM,threads,use_mkl_sparse);
            KRY_VEC rhs(threads),sol(threads);
            rhs.v=b;
            sol.v.Resize(size);
            timer("prep krylov");

            ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
            //sys.Test_System(sol);
            MINRES<T> mr;
            timer("test krylov");
            bool converged=mr.Solve(sys,sol,rhs,av,1e-12,0,kn);
            timer("krylov solve");
            if(!converged) LOG::printf("KRYLOV SOLVER DID NOT CONVERGE.\n");
            krylov_sol=sol.v;
        }
        else if(dump_sysbin)
        {
            ARRAY<int> col(SM.A.m);
            ARRAY<T> entries(SM.A.m);
            for(int i=0;i<SM.A.m;i++)
            {
                col(i)=SM.A(i).j;
                entries(i)=SM.A(i).a;
            }
            auto write_bin=[&output_dir](const auto& arr,const char* name)
            {
                // size values
                std::FILE* file=std::fopen((output_dir+name).c_str(),"wb");
                std::fwrite(&arr.m,sizeof(arr.m),1,file);
                std::fwrite(arr.base_pointer,sizeof(arr(0)),arr.m,file);
                std::fclose(file);
            };
            write_bin(col,"/col.bin");
            write_bin(SM.offsets,"/offsets.bin");
            write_bin(b,"/rhs.bin");
            // dim nnz values
            std::FILE* entries_file=std::fopen((output_dir+"/entries.bin").c_str(),"wb");
            std::fwrite(&SM.m,sizeof(SM.m),1,entries_file);
            std::fwrite(&entries.m,sizeof(entries.m),1,entries_file);
            std::fwrite(entries.base_pointer,sizeof(T),entries.m,entries_file);
            std::fclose(entries_file);
        }
    }

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

    cem.matrix_cache.Init(cache_pattern,cache_size);
    // Let job scheduler spawn threads.
#ifdef USE_OPENMP
        omp_set_num_threads(1);
#endif
#ifdef USE_MKL
        mkl_set_num_threads(1);
#endif
    cem.Execute_Jobs(threads);

    timer("exec jobs");

    mc.Transform_Solution(cem,false,false);
    if(an) an->Check_Analytic_Solution(!quiet,pinv);
    if(!sol_file.empty())
    {
        SOLUTION_FEM<TV> sol;
        Read_From_File(sol_file,sol);
        sol.Prepare_Hierarchy();
        auto fv=[&sol](const TV& X){return sol.Velocity(X);};
        auto fp=[&sol](const TV& X){return sol.Pressure(X);};
        Check_Solution<T,TV::m>(mc,fv,fp,!quiet,pinv);
    }

    if(!quiet)
    {
        for(BLOCK_ID b(0);b<cl.blocks.m;b++)
            debug.Visualize_Solution(mc.rhs_block_list(b),b,true);
        Flush_Frame<TV>("solution");
        mc.Dump_World_Space_Vector("x");
        debug.Visualize_Flat_Dofs();
    }

    if(use_krylov || dump_sysbin)
    {
        ARRAY<int,BLOCK_ID> first[3];
        int size=mc.Compute_Global_Dof_Mapping(first);
        ARRAY<T> world_sol;
        mc.Dump_World_Space_Vector(first,size,world_sol);
        if(use_krylov)
            LOG::printf("Krylov diff %P\n",(world_sol-krylov_sol).Max_Abs());
        else if(dump_sysbin)
        {
            std::FILE* file=std::fopen((output_dir+"/x.bin").c_str(),"wb");
            std::fwrite(&world_sol.m,sizeof(world_sol.m),1,file);
            std::fwrite(world_sol.base_pointer,sizeof(T),world_sol.m,file);
            std::fclose(file);
        }
    }

    for(int i=1;i<tm.m;i++)
        LOG::printf("%20s %5.0f ms\n",tm(i).y,
        std::chrono::duration_cast<std::chrono::duration<double> >(tm(i).x-tm(i-1).x).count()*1000);

    if(dump_solution)
    {
        SOLUTION_FEM<TV> solution;
        solution.Build(mc,an);
        LOG::printf("sol maxv %P maxp %P\n",solution.Max_Velocity_Magnitude(),solution.Max_Pressure());
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
