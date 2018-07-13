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

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,2> TV;
typedef VECTOR<int,TV::m> TV_INT;

struct MATRIX_ID
{
    enum ID_TYPE {orig_block} type;
    int a,b,dir;

    bool operator<(const MATRIX_ID& m) const
    {
        if(type!=m.type) return type<m.type;
        if(a!=m.a) return a<m.a;
        if(b!=m.b) return b<m.b;
        return dir<m.dir;
    }
};



int main(int argc, char* argv[])
{
    T mu=1;
    T dx=.1;

    std::string pipe_file;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-mu",&mu,"mu","viscosity");
    parse_args.Extra(&pipe_file,"file","file describing pipes");
    parse_args.Parse();

    PARSE_DATA<TV> pd;
    pd.Parse_Input(pipe_file);
    
    GRID<TV> grid(pd.box_size,RANGE<TV>(TV(),TV(pd.box_size)*dx),true);
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE(0.f),grid,"output");
    Flush_Frame<TV>("init");

    FLUID_LAYOUT<TV> fl(grid);
    
    fl.Compute(pd);
    fl.Dump_Layout();
    Flush_Frame<TV>("grid setup");
    fl.Dump_Dofs();
    Flush_Frame<TV>("grid dofs");
    fl.Dump_Blocks();
    Flush_Frame<TV>("grid blocks");
    fl.Dump_Block_Types();
    Flush_Frame<TV>("block types");

    SYSTEM_MATRIX_HELPER<T> MH;
    ARRAY<T> rhs_vector,sol_vector;
    Compute_Full_Matrix(grid,MH,rhs_vector,fl,mu);
    Solve_And_Display_Solution(grid,fl,MH,rhs_vector,&sol_vector);

    CACHED_ELIMINATION_MATRIX<T> elim_mat;
    elim_mat.orig_sizes.Resize(fl.blocks.m);
    for(int i=0;i<fl.blocks.m;i++) elim_mat.orig_sizes(i)=fl.blocks(i).num_dofs;
    Setup_Block_Map(elim_mat,fl);

    elim_mat.Fill_Blocks(fl.dof_map,MH.data,rhs_vector);

    elim_mat.Unpack_Vector(fl.dof_map,elim_mat.test_sol,sol_vector);

    elim_mat.Fill_Orig_Rows();
    elim_mat.Test_State();

    for(int i=2;i<elim_mat.block_list.m;i++)
        OCTAVE_OUTPUT<T>(LOG::sprintf("M-%i.txt",i).c_str()).Write("M",*elim_mat.block_list(i));

    {
        OCTAVE_OUTPUT<T> oo("block.txt");
        ARRAY<T> l(sol_vector.m),r(sol_vector.m);
        ARRAY<ARRAY<T>*> ll,rr;
        int b=l.m;
        oo.Begin_Sparse_Matrix("N",l.m,b);
        l*=0;
        r*=0;
        for(int i=0;i<b;i++){
            r(i)=1;
            l*=0;
            elim_mat.Unpack_Vector(fl.dof_map,ll,l);
            elim_mat.Unpack_Vector(fl.dof_map,rr,r);
            elim_mat.Add_Times(ll,rr);
            elim_mat.Pack_Vector(fl.dof_map,l,ll);
            r(i)=0;
            oo.Append_Sparse_Column(l);}

        oo.End_Sparse_Matrix();
    }

    for(int i=0;i<pd.pts.m;i++)
        elim_mat.Eliminate_Row(i);

    while(1)
    {
        std::map<int,int> counts;
        for(int i=0;i<elim_mat.rows.m;i++)
            if(elim_mat.valid_row(i) && elim_mat.rows(i).m<=3)
                counts[elim_mat.Get_Block_Lazy(i,i)]++;
        int best=-1,best_mat=-1;
        for(auto i:counts)
            if(i.second>best){
                best=i.second;
                best_mat=i.first;}

        if(best<=0) break;

        printf("Eliminate: %i  (%i)\n",best_mat,best);
        for(int i=0;i<elim_mat.rows.m;i++)
            if(elim_mat.valid_row(i) && elim_mat.rows(i).m<=3)
                if(elim_mat.Get_Block_Lazy(i,i)==best_mat)
                    elim_mat.Eliminate_Row(i);
    }

    PHYSBAM_ASSERT(elim_mat.elimination_order.m==elim_mat.orig_sizes.m);
    elim_mat.Back_Solve();

    ARRAY<T,FACE_INDEX<TV::m> > face_velocity(grid,1);
    for(FACE_ITERATOR<TV> it(grid,1);it.Valid();it.Next()){
        auto& uf=fl.used_faces(it.Full_Index());
        if(uf.type!=fluid) continue;
        face_velocity(it.Full_Index())=(*elim_mat.rhs(uf.block_id))(uf.block_dof);}
    Flush_Frame(face_velocity,"elim solve");
    
    return 0;
}

