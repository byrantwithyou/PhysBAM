//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <Tools/Krylov_Solvers/MINRES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_RANGE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Analytic_Tests/ANALYTIC_SCALAR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_VECTOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "CACHED_ELIMINATION_MATRIX.h"
#include "FEM_MESHING_TESTS.h"
#include "FEM_TABLE.h"
#include "FLAT_SYSTEM.h"
#include "FLAT_SYSTEM_FEM.h"
#include "FLAT_SYSTEM_FEM_EXTRUDED.h"
#include "FLUID_LAYOUT.h"
#include "FLUID_LAYOUT_FEM.h"
#include "FLUID_LAYOUT_FEM_EXTRUDED.h"
#include "EXECUTE_HELPER.h"
#include <chrono>
#include <fstream>
#include <map>
#include <set>
#include <string>

using namespace PhysBAM;

typedef float RW;
typedef double T;

void Run_Meshing_Tests(int seed)
{
    typedef VECTOR<T,2> TV;
    LOG::printf("seed: %d\n",seed);
    Test_Degree3_Joint<TV>(default_joint,0.2,20,seed);
    Test_Degree2_Joint<TV>(default_joint,-pi+pi/20,pi-pi/20,pi/40);
    Test_Degree2_Joint<TV>(corner_joint,-pi+pi/20,pi-pi/20,pi/40);
    Test_Degree2_Circle<TV>(default_joint,0.2,0.9,0.1);
    Test_Degree2_Circle<TV>(corner_joint,0.2,0.9,0.1);
}

template<int d>
void Run_FEM(PARSE_ARGS& parse_args)
{
    typedef VECTOR<T,2> TV2;
    typedef VECTOR<T,d> TV;
    T mu=1;
    std::string pipe_file;
    bool run_tests=false,use_krylov=false;
    int seed=time(0),threads=1;
    VIEWER_DIR viewer_dir("output");
    parse_args.Add("-mu",&mu,"mu","viscosity");
    parse_args.Add("-tests",&run_tests,"run FEM tests");
    parse_args.Add("-k",&use_krylov,"solve with Krylov method");
    parse_args.Add("-o",&viewer_dir.output_directory,"dir","output dir");
    parse_args.Add("-threads",&threads,"num","number of threads to use");
    parse_args.Parse(true);
    if(!run_tests)
        parse_args.Extra(&pipe_file,"file","file describing pipes");
    else parse_args.Add("-seed",&seed,"seed","random seed");
    parse_args.Parse();

    VIEWER_OUTPUT vo(STREAM_TYPE(0.f),viewer_dir);
    Use_Debug_Particles<TV>();

    if(d==2 && run_tests){
        Run_Meshing_Tests(seed);
        return;}

    PARSE_DATA_FEM<TV2,TV> pd;
    pd.Parse_Input(pipe_file);

    FLUID_LAYOUT_FEM<TV> fl;
    fl.Dump_Input(pd);
    Flush_Frame("init");
    fl.Compute(pd);
    fl.Print_Statistics();
    fl.Dump_Mesh();
    Flush_Frame("meshing");
    fl.Dump_Layout();
    Flush_Frame("blocks");
    fl.Dump_Node_Blocks();
    Flush_Frame("nodes");
    fl.Dump_Edge_Blocks();
    Flush_Frame("edges");
    fl.Dump_Dofs();
    Flush_Frame("dofs");
    
    if(d==3) return;
    ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> > coded_entries;
    ARRAY<T,CODE_ID> code_values;
    ARRAY<T,DOF_ID> rhs_vector(fl.num_dofs);
    Generate_Discretization(coded_entries,code_values,fl,pd,mu,rhs_vector);

    ARRAY<T,DOF_ID> sol;
    if(pd.analytic_velocity && pd.analytic_pressure)
    {
        ARRAY<T,DOF_ID> res(rhs_vector);
        sol.Resize(rhs_vector.m);
        for(PARTICLE_ID p(0);p<fl.Number_Particles();p++)
        {
            TV X=fl.X(p),U=pd.Velocity(X,BC_ID());
            sol(fl.pressure_dofs(p))=pd.Pressure(X);
            DOF_ID dof=fl.vel_node_dofs(p);
            if(dof<DOF_ID()) continue;
            for(int i=0;i<2;i++)
                sol(dof+i)=U(i);
        }
        for(EDGE_ID e(0);e<fl.Number_Edges();e++)
        {
            DOF_ID dof=fl.vel_edge_dofs(e);
            if(dof<DOF_ID()) continue;
            auto ed=fl.Edge(e);
            TV X=(T).5*(fl.X(ed.x)+fl.X(ed.y));
            TV U=pd.Velocity(X,BC_ID());
            for(int i=0;i<2;i++)
                sol(dof+i)=U(i);
        }
        //LOG::printf("sol: %P\n",sol);
        //LOG::printf("rhs: %P\n",rhs_vector);
        for(auto t:coded_entries)
            res(t.x)-=code_values(t.z)*sol(t.y);
        //LOG::printf("res: %P\n",res);
        LOG::printf("rhs error: %P\n",res.Max_Abs());
    }

    SYSTEM_MATRIX_HELPER<T> MH;
    for(auto& e:coded_entries) MH.data.Append({Value(e.x),Value(e.y),code_values(e.z)});
    ARRAY<T,DOF_ID> sol_vector;
    if(use_krylov)
        Solve_And_Display_Solution(fl,pd,MH,rhs_vector,&sol_vector);

    CACHED_ELIMINATION_MATRIX<T> elim_mat;
    elim_mat.orig_sizes.Resize(Value(fl.blocks.m));
    for(BLOCK_ID i(0);i<fl.blocks.m;i++) elim_mat.orig_sizes(Value(i))=fl.blocks(i).num_dofs;
    elim_mat.Fill_Blocks(fl.dof_map,coded_entries,code_values,rhs_vector);
    elim_mat.Unpack_Vector(fl.dof_map,elim_mat.test_sol,sol_vector);
    elim_mat.Fill_Orig_Rows();
    elim_mat.Reduce_Rows_By_Frequency(0,Value(fl.blocks.m)-Value(fl.pipes.m),Value(fl.blocks.m));
    elim_mat.Reduce_Rows_By_Frequency(Value(fl.blocks.m)-Value(fl.pipes.m),Value(fl.blocks.m),3);
    elim_mat.Full_Reordered_Elimination();
    elim_mat.Back_Solve();
    Execute_Helper(&elim_mat,threads);
    ARRAY<T,DOF_ID> elim_sol;
    elim_mat.Pack_Vector(fl.dof_map,elim_sol,elim_mat.rhs);
    if(use_krylov) LOG::printf("ANS DIFF: %g\n",(elim_sol-sol_vector).Max_Abs());

    if(pd.analytic_velocity && pd.analytic_pressure){
        ARRAY<T,DOF_ID> ksol_error(sol_vector),elim_sol_error(elim_sol);
        ksol_error-=sol;
        elim_sol_error-=sol;
        //LOG::printf("sol: %P\n",sol_vector);
        LOG::printf("elim sol error: %P\n",elim_sol_error.Max_Abs());
        if(use_krylov) LOG::printf("krylov sol error: %P\n",ksol_error.Max_Abs());}
    LOG::Instance()->Copy_Log_To_File(viewer_dir.output_directory+"/common/log.txt",false);
}

template<int d>
void Run(PARSE_ARGS& parse_args)
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,TV::m> TV_INT;

    T mu=1;
    T dx=.1;
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();

    int threads=1;
    bool quiet=false,use_krylov=false,print_system=false;
    std::string pipe_file;
    parse_args.Add("-mu",&mu,"mu","viscosity");
    parse_args.Add("-q",&quiet,"disable diagnostics; useful for timing");
    parse_args.Add("-k",&use_krylov,"solve with Krylov method");
    parse_args.Add("-p",&print_system,"dump the system to be solved");
    parse_args.Add("-threads",&threads,"num","number of threads to use");
    parse_args.Extra(&pipe_file,"file","file describing pipes");
    parse_args.Parse();

    PARSE_DATA<TV> pd;
    pd.Parse_Input(pipe_file);
    
    GRID<TV> grid(pd.box_size,RANGE<TV>(TV(),TV(pd.box_size)*dx),true);
    VIEWER_DIR viewer_dir("output");
    VIEWER_OUTPUT vo(STREAM_TYPE(0.f),viewer_dir);
    Use_Debug_Particles<TV>();
    if(!quiet) Flush_Frame("init");

    FLUID_LAYOUT<TV> fl(grid);
    fl.quiet=quiet;
    
    fl.Compute(pd);
    if(!quiet){
        fl.Dump_Layout();
        Flush_Frame("grid setup");
        fl.Dump_Dofs();
        Flush_Frame("grid dofs");
        fl.Dump_Blocks();
        Flush_Frame("grid blocks");}

    SYSTEM_MATRIX_HELPER<T> MH;
    ARRAY<TRIPLE<DOF_ID,DOF_ID,CODE_ID> > coded_entries;
    ARRAY<T,CODE_ID> code_values;
    ARRAY<T,DOF_ID> rhs_vector,sol_vector;
    Compute_Full_Matrix(grid,coded_entries,code_values,rhs_vector,fl,mu);
    for(auto e:coded_entries) MH.data.Append({Value(e.x),Value(e.y),code_values(e.z)});
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    if(use_krylov || print_system)
        Solve_And_Display_Solution(grid,fl,MH,rhs_vector,&sol_vector,
            use_krylov,print_system,quiet);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    LOG::printf("total dofs: %i\n",rhs_vector.m);
    LOG::printf("blocks: %i\n",fl.blocks.m);
    
    CACHED_ELIMINATION_MATRIX<T> elim_mat;
    elim_mat.quiet=quiet;
    elim_mat.orig_sizes.Resize(fl.blocks.m);
    for(int i=0;i<fl.blocks.m;i++) elim_mat.orig_sizes(i)=fl.blocks(i).num_dofs;

    elim_mat.Fill_Blocks(fl.dof_map,coded_entries,code_values,rhs_vector);

    elim_mat.Unpack_Vector(fl.dof_map,elim_mat.test_sol,sol_vector);

    elim_mat.Fill_Orig_Rows();

    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
    elim_mat.Reduce_Rows_By_Frequency(0,fl.num_vertex_blocks,fl.blocks.m);
    elim_mat.Reduce_Rows_By_Frequency(fl.num_vertex_blocks,fl.blocks.m,3);
    elim_mat.Full_Reordered_Elimination();
    elim_mat.Back_Solve();
    std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
    Execute_Helper(&elim_mat,threads);
    std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();

    ARRAY<T,DOF_ID> elim_sol;
    elim_mat.Pack_Vector(fl.dof_map,elim_sol,elim_mat.rhs);
    if(!quiet && use_krylov) LOG::printf("ANS DIFF: %g\n",(elim_sol-sol_vector).Max_Abs());
    
    if(!quiet){
        ARRAY<T,FACE_INDEX<TV::m> > face_velocity(grid,1);
        for(FACE_ITERATOR<TV> it(grid,1);it.Valid();it.Next()){
            auto& uf=fl.used_faces(it.Full_Index());
            if(uf.type!=fluid) continue;
            face_velocity(it.Full_Index())=elim_mat.vector_list(elim_mat.rhs(uf.block_id))(uf.block_dof);}
        Flush_Frame("elim solve");}
    std::chrono::steady_clock::time_point t6 = std::chrono::steady_clock::now();

    printf("grid setup    %5.0f ms\n",
        std::chrono::duration_cast<std::chrono::duration<double> >(t1-t0).count()*1000);
    printf("minres solve  %5.0f ms\n",
        std::chrono::duration_cast<std::chrono::duration<double> >(t2-t1).count()*1000);
    printf("setup blocks  %5.0f ms\n",
        std::chrono::duration_cast<std::chrono::duration<double> >(t3-t2).count()*1000);
    printf("elimination   %5.0f ms\n",
        std::chrono::duration_cast<std::chrono::duration<double> >(t4-t3).count()*1000);
    printf("execute jobs  %5.0f ms\n",
        std::chrono::duration_cast<std::chrono::duration<double> >(t5-t4).count()*1000);
    printf("finish        %5.0f ms\n",
        std::chrono::duration_cast<std::chrono::duration<double> >(t6-t5).count()*1000);
}

int main(int argc, char* argv[])
{
    bool use_3d=false;
    bool use_fem=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"use 3D");
    parse_args.Add("-fem",&use_fem,"use FEM");
    parse_args.Parse(true);
    LOG::Initialize_Logging(false,false,1<<30,true);
    if(use_fem){
        if(use_3d) Run_FEM<3>(parse_args);
        else Run_FEM<2>(parse_args);}
    else{
        if(use_3d) Run<3>(parse_args);
        else Run<2>(parse_args);}

    LOG::Finish_Logging();
    return 0;
}

