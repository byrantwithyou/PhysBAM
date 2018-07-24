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
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include "CACHED_ELIMINATION_MATRIX.h"
#include "FLAT_SYSTEM.h"
#include "FLUID_LAYOUT.h"
#include "FLUID_LAYOUT_FEM.h"

using namespace PhysBAM;

typedef float RW;
typedef double T;

void Run_FEM(PARSE_ARGS& parse_args)
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    T mu=1;

    std::string pipe_file;
    parse_args.Add("-mu",&mu,"mu","viscosity");
    parse_args.Extra(&pipe_file,"file","file describing pipes");
    parse_args.Parse();
    GRID<TV> grid;
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE(0.f),grid,"output");

    PARSE_DATA_FEM<TV> pd;
    pd.Parse_Input(pipe_file);

    FLUID_LAYOUT_FEM<TV> fl;
    fl.Compute(pd);
    Flush_Frame<TV>("init");
    fl.Dump_Mesh();
    Flush_Frame<TV>("meshing");
    fl.Dump_Layout();
    Flush_Frame<TV>("blocks");
}

template<int d>
void Run(PARSE_ARGS& parse_args)
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,TV::m> TV_INT;

    T mu=1;
    T dx=.1;

    int threads=1;
    bool quiet=false;
    std::string pipe_file;
    parse_args.Add("-mu",&mu,"mu","viscosity");
    parse_args.Add("-q",&quiet,"disable diagnostics; useful for timing");
    parse_args.Add("-threads",&threads,"num","number of threads to use");
    parse_args.Extra(&pipe_file,"file","file describing pipes");
    parse_args.Parse();

    PARSE_DATA<TV> pd;
    pd.Parse_Input(pipe_file);
    
    GRID<TV> grid(pd.box_size,RANGE<TV>(TV(),TV(pd.box_size)*dx),true);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE(0.f),grid,"output");
    if(!quiet) Flush_Frame<TV>("init");

    FLUID_LAYOUT<TV> fl(grid);
    
    fl.Compute(pd);
    if(!quiet){
        fl.Dump_Layout();
        Flush_Frame<TV>("grid setup");
        fl.Dump_Dofs();
        Flush_Frame<TV>("grid dofs");
        fl.Dump_Blocks();
        Flush_Frame<TV>("grid blocks");}

    SYSTEM_MATRIX_HELPER<T> MH;
    ARRAY<VECTOR<int,3> > coded_entries;
    ARRAY<T> code_values;
    ARRAY<T> rhs_vector,sol_vector;
    Compute_Full_Matrix(grid,coded_entries,code_values,rhs_vector,fl,mu);
    for(auto e:coded_entries) MH.data.Append({e.x,e.y,code_values(e.z)});
    if(!quiet) Solve_And_Display_Solution(grid,fl,MH,rhs_vector,&sol_vector);
    LOG::printf("total dofs: %i\n",rhs_vector.m);
    LOG::printf("blocks: %i\n",fl.blocks.m);
    
    CACHED_ELIMINATION_MATRIX<T> elim_mat;
    elim_mat.quiet=quiet;
    elim_mat.orig_sizes.Resize(fl.blocks.m);
    for(int i=0;i<fl.blocks.m;i++) elim_mat.orig_sizes(i)=fl.blocks(i).num_dofs;

    elim_mat.Fill_Blocks(fl.dof_map,coded_entries,code_values,rhs_vector);

    elim_mat.Unpack_Vector(fl.dof_map,elim_mat.test_sol,sol_vector);

    elim_mat.Fill_Orig_Rows();

    elim_mat.Reduce_Rows_By_Frequency(0,fl.num_vertex_blocks,fl.blocks.m);
    elim_mat.Reduce_Rows_By_Frequency(fl.num_vertex_blocks,fl.blocks.m,3);
    elim_mat.Full_Reordered_Elimination();
    elim_mat.Back_Solve();
    elim_mat.Execute_Jobs(threads);

    ARRAY<T> elim_sol;
    elim_mat.Pack_Vector(fl.dof_map,elim_sol,elim_mat.rhs);
    if(!quiet) LOG::printf("ANS DIFF: %g\n",(elim_sol-sol_vector).Max_Abs());
    
    if(!quiet){
        ARRAY<T,FACE_INDEX<TV::m> > face_velocity(grid,1);
        for(FACE_ITERATOR<TV> it(grid,1);it.Valid();it.Next()){
            auto& uf=fl.used_faces(it.Full_Index());
            if(uf.type!=fluid) continue;
            face_velocity(it.Full_Index())=elim_mat.vector_list(elim_mat.rhs(uf.block_id))(uf.block_dof);}
        Flush_Frame(face_velocity,"elim solve");}
}

int main(int argc, char* argv[])
{
    bool use_3d=false;
    bool use_fem=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"use 3D");
    parse_args.Add("-fem",&use_fem,"use FEM");
    parse_args.Parse(true);

    if(use_3d) Run<3>(parse_args);
    else{
        if(use_fem) Run_FEM(parse_args);
        else Run<2>(parse_args);}

    return 0;
}

