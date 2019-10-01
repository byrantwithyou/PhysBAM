//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Math_Tools/RANGE_ITERATOR.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_RANGE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <map>
#include "FLUID_LAYOUT.h"

namespace PhysBAM{
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void FLUID_LAYOUT<TV>::
Compute(const PARSE_DATA<TV>& pd)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    used_faces.Resize(grid,1,use_init,{nodof,-1,-1,0});
    used_cells.Resize(grid.Domain_Indices(1),use_init,{nodof,-1,-1,0,DOF_ID(-1)});
    
    for(auto& i:pd.pipes)
    {
        RANGE<TV_INT> box=pd.Pipe_Full_Range(i);
        for(FACE_RANGE_ITERATOR<TV::m> it(box,RF::ghost);it.Valid();it.Next())
            used_faces(it.face).type=wall;
        for(RANGE_ITERATOR<TV::m> it(box,RI::ghost);it.Valid();it.Next())
            used_cells(it.index).type=wall;
    }

    for(auto& i:pd.pts)
    {
        for(FACE_RANGE_ITERATOR<TV::m> it(i.box,RF::ghost);it.Valid();it.Next())
            used_faces(it.face).type=wall;
        for(FACE_RANGE_ITERATOR<TV::m> it(i.box,RF::skip_outer);it.Valid();it.Next())
            used_faces(it.face).type=fluid;
        for(RANGE_ITERATOR<TV::m> it(i.box);it.Valid();it.Next())
            used_cells(it.index).type=fluid;
    }
    
    for(auto& i:pd.pipes)
    {
        RANGE<TV_INT> box=pd.Pipe_Full_Range(i);
        for(FACE_RANGE_ITERATOR<TV::m> it(box,RF::skip_outer);it.Valid();it.Next())
            used_faces(it.face).type=fluid;
        for(RANGE_ITERATOR<TV::m> it(box);it.Valid();it.Next())
            used_cells(it.index).type=fluid;
    }

    for(auto& i:pd.pts)
    {
        if(i.bc_type==dirichlet)
        {
            for(RANGE_ITERATOR<TV::m> it(i.box,1,0,RI::ghost|RI::omit_corners,i.bc_side);it.Valid();it.Next())
            {
                used_cells(it.index)={dirichlet,-1,-1,i.bc_value,DOF_ID(-1)};
            }
            for(FACE_RANGE_ITERATOR<TV::m> it(i.box,0,0,RF::ghost,i.bc_side);it.Valid();it.Next())
            {
                used_faces(it.face)={fluid,-1,-1,0,DOF_ID(-1)};
            }
        }
        if(i.bc_type==wall)
        {
            for(FACE_RANGE_ITERATOR<TV::m> it(i.box,0,0,RF::ghost,i.bc_side);it.Valid();it.Next())
            {
                used_faces(it.face).type=i.bc_type;
                used_faces(it.face).bc_value=pd.Inflow_BC_Value(grid.Face(it.face),i,grid);
            }
        }
        Allocate_Cross_Section_Blocks_Cells(i.box,0);
    }
    num_vertex_blocks=blocks.m;
    if(!quiet) Dump_Blocks();
    if(!quiet) Flush_Frame<TV>("vertex cells");
    for(auto& i:pd.pipes)
    {
        int dir=pd.Pipe_Dir(i);
        RANGE<TV_INT> box=pd.Pipe_Inner_Range(i);
        Allocate_Cross_Section_Blocks_Cells(box,dir);
        Allocate_Cross_Section_Blocks_Faces(box,dir);
    }
    if(!quiet) Dump_Blocks();
    if(!quiet) Flush_Frame<TV>("pipes");

    for(auto& i:pd.pts)
    {
        Allocate_Cross_Section_Blocks_Faces(i.box,0);
    }
    if(!quiet) Dump_Blocks();
    if(!quiet) Flush_Frame<TV>("vertex faces");

    for(FACE_RANGE_ITERATOR<TV::m> it(used_faces.domain_indices);it.Valid();it.Next())
    {
        auto& f=used_faces(it.face);
        if(f.type==fluid)
        {
            f.block_dof=blocks(f.block_id).num_dofs++;
            f.global_id=dof_map.Append({f.block_id,f.block_dof});
        }
    }
    for(RANGE_ITERATOR<TV::m> it(used_cells.domain);it.Valid();it.Next())
    {
        auto& c=used_cells(it.index);
        if(c.type==fluid)
        {
            c.block_dof=blocks(c.block_id).num_dofs++;
            c.global_id=dof_map.Append({c.block_id,c.block_dof});
        }
    }
}
//#####################################################################
// Function Dump_Layout
//#####################################################################
template<class TV> void FLUID_LAYOUT<TV>::
Dump_Layout() const
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    for(FACE_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& f=used_faces(it.face);
        if(f.type==fluid) Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,1));
        else if(f.type==wall) Add_Debug_Particle(it.Location(),VECTOR<T,3>(0,1,1));
        else PHYSBAM_ASSERT(f.type==nodof);
    }
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& c=used_cells(it.index);
        if(c.type==fluid) Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,1,1));
        else if(c.type==dirichlet) Add_Debug_Particle(it.Location(),VECTOR<T,3>(1,0,1));
        else PHYSBAM_ASSERT(c.type==nodof);
    }
}
//#####################################################################
// Function Dump_Dofs
//#####################################################################
template<class TV> void FLUID_LAYOUT<TV>::
Dump_Dofs() const
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    for(FACE_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& f=used_faces(it.face);
        std::string s=LOG::sprintf("%i",f.global_id);
        if(f.type==fluid) Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,1,1));
        else if(f.type==wall) Add_Debug_Text(it.Location(),s,VECTOR<T,3>(0,1,1));
        else Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,0,0));
    }
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& c=used_cells(it.index);
        std::string s=LOG::sprintf("%i",c.global_id);
        if(c.type==fluid) Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,1,1));
        else if(c.type==dirichlet) Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,0,1));
        else Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,0,0));
    }
}
//#####################################################################
// Function Dump_Blocks
//#####################################################################
template<class TV> void FLUID_LAYOUT<TV>::
Dump_Blocks() const
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    for(FACE_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& f=used_faces(it.face);
        if(f.type!=fluid)  continue;
        std::string s=LOG::sprintf("%i",f.block_id);
        Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,1,1));
    }
    for(CELL_ITERATOR<TV> it(grid,1);it.Valid();it.Next())
    {
        auto& c=used_cells(it.index);
        if(c.type!=fluid)  continue;
        std::string s=LOG::sprintf("%i",c.block_id);
        Add_Debug_Text(it.Location(),s,VECTOR<T,3>(1,1,1));
    }
}
//#####################################################################
// Function Assign_Cross_Section_Blocks
//#####################################################################
template<class TV> void FLUID_LAYOUT<TV>::
Allocate_Cross_Section_Blocks_Cells(const RANGE<TV_INT>& box,int dir)
{
    int next_block=blocks.m;
    for(RANGE_ITERATOR<TV::m> it(box);it.Valid();it.Next())
        used_cells(it.index).block_id=it.index(dir)-box.min_corner(dir)+next_block;
    int num_blocks=box.Edge_Lengths()(dir);
    for(int i=0;i<num_blocks;i++) blocks.Append({0});
}
//#####################################################################
// Function Assign_Cross_Section_Blocks
//#####################################################################
template<class TV> void FLUID_LAYOUT<TV>::
Allocate_Cross_Section_Blocks_Faces(const RANGE<TV_INT>& box,int dir)
{
    TV_INT base=box.Center();
    base(dir)=box.min_corner(dir);
    for(FACE_RANGE_ITERATOR<TV::m> it(box);it.Valid();it.Next())
    {
        if(used_faces(it.face).type!=fluid) continue;
        int b0=used_cells(it.face.index).block_id;
        TV_INT cell=it.face.index;
        cell(it.face.axis)--;
        int b1=used_cells(cell).block_id;
        if(b1>=0)
            if(b0<0 || (it.face.axis==dir && cell(!dir)>=base(!dir)))
                b0=b1;
        if(used_faces(it.face).block_id<0)
            used_faces(it.face).block_id=b0;
    }
}

template struct FLUID_LAYOUT<VECTOR<double,2> >;
template struct FLUID_LAYOUT<VECTOR<double,3> >;
}
