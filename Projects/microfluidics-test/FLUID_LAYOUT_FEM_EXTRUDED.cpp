//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/FRAME.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "FLUID_LAYOUT_FEM.h"
#include "FLUID_LAYOUT_FEM_EXTRUDED.h"

namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<typename T> FLUID_LAYOUT_FEM<VECTOR<T,3> >::
FLUID_LAYOUT_FEM(): vol_hidden(*new TETRAHEDRALIZED_VOLUME<T>)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<typename T> FLUID_LAYOUT_FEM<VECTOR<T,3> >::
~FLUID_LAYOUT_FEM()
{
    delete &vol_hidden;
}
//#####################################################################
// Function Extrude
//#####################################################################
template<typename T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Extrude(int layers,T dz)
{
    vol_hidden.particles.Delete_All_Elements();
    vol_hidden.particles.Add_Elements(fl2.area_hidden.particles.number*(layers+1));
    TV3 dx=dz*extrude_dir;
    int m=fl2.area_hidden.particles.number,n=layers*m;
    for(int i=0;i<m;i++) vol_hidden.particles.X(i)=TV3(fl2.area_hidden.particles.X(i));
    for(int i=0;i<n;i++) vol_hidden.particles.X(i+m)=vol_hidden.particles.X(i)+dx;
    //mesh.Initialize_Swept_Mesh(ta.mesh,layers);
    for(int l=0,particle_offset=0;l<layers;l++){
        for(int e=0;e<fl2.area_hidden.mesh.elements.m;e++){
            BLOCK_ID bid=fl2.elem_data(TRIANGLE_ID(e)).block_id;
            VECTOR<int,3> e0=fl2.area_hidden.mesh.elements(e)+particle_offset,e1=e0+fl2.area_hidden.mesh.number_nodes;
            int a=e0.Arg_Min(),c=e0.Arg_Max(),b=3-a-c;
            VECTOR<int,4> t0=e0.Append(e1(a)),t1(t0);
            t1(a)=e1(b);
            VECTOR<int,4> t2=t1;
            t2(b)=e1(c);
            exchange(t1(0),t1(1));
            vol_hidden.mesh.elements.Append(t0);
            vol_hidden.mesh.elements.Append(t1);
            vol_hidden.mesh.elements.Append(t2);
            elem_blocks.Append(bid);
            elem_blocks.Append(bid);
            elem_blocks.Append(bid);}
        particle_offset+=fl2.area_hidden.mesh.number_nodes;}
    vol_hidden.mesh.Set_Number_Nodes(fl2.area_hidden.mesh.number_nodes*(layers+1));
}
//#####################################################################
// Function Update_BC_Planar
//#####################################################################
template<class T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Update_BC_Planar(const PARSE_DATA_FEM<TV,TV3>& pd)
{
    for(EDGE_ID i(0);i<fl2.Number_Edges();i++)
        if(fl2.Edge_Triangles(i).m==1){
            auto key=fl2.Edge(i).Sorted();
            if(!fl2.bc_map.Contains(key)){
                if(BC_ID* pbc=fl2.particle_bc_map.Get_Pointer(key.x)) *pbc=pd.wall_bc;
                else fl2.particle_bc_map.Set(key.x,pd.wall_bc);
                if(BC_ID* pbc=fl2.particle_bc_map.Get_Pointer(key.y)) *pbc=pd.wall_bc;
                else fl2.particle_bc_map.Set(key.y,pd.wall_bc);
                fl2.bc_map.Set({key.x,key.y},pd.wall_bc);}}
}
//#####################################################################
// Function Assign_Edge_Blocks_Planar
//#####################################################################
template<class T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Assign_Edge_Blocks_Planar()
{
    fl2.edge_blocks.Resize(fl2.Number_Edges(),init_all,BLOCK_ID(-1));
    for(EDGE_ID ei(0);ei<fl2.Number_Edges();ei++){
        auto neighbor_tris=fl2.Edge_Triangles(ei);
        PHYSBAM_ASSERT(neighbor_tris.m==1 || neighbor_tris.m==2);
        BLOCK_ID bid=fl2.elem_data(neighbor_tris(0)).block_id;
        if(neighbor_tris.m==2 && fl2.elem_data(neighbor_tris(1)).greed>fl2.elem_data(neighbor_tris(0)).greed)
            bid=fl2.elem_data(neighbor_tris(1)).block_id;
        fl2.edge_blocks(ei)=bid;}
}
//#####################################################################
// Function Allocate_Dofs
//#####################################################################
template<class T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Allocate_Dofs(const PARSE_DATA_FEM<TV,TV3>& pd)
{
    ARRAY<int,BLOCK_ID> block_vel_dofs(blocks.m),block_pressure_dofs(blocks.m);
    for(PARTICLE_ID i(0);i<Number_Particles();i++){
        BLOCK_ID bid=Node_Block(i);
        block_pressure_dofs(bid)++;
        BC_ID bc=Node_BC(i);
        if(bc==BC_ID(-1) || pd.bc(bc).type!=dirichlet_v)
            block_vel_dofs(bid)++;}
    for(EDGE_ID ei(0);ei<Number_Edges();ei++){
        VECTOR<PARTICLE_ID,2> edge=Edge(ei);
        BLOCK_ID bid=Edge_Block(edge(0),edge(1));
        BC_ID bc=Edge_BC(ei);
        if(bc==BC_ID(-1) || pd.bc(bc).type!=dirichlet_v)
            block_vel_dofs(bid)++;}
    for(BLOCK_ID i(1);i<blocks.m;i++){
        block_vel_dofs(i)+=block_vel_dofs(i-1);
        block_pressure_dofs(i)+=block_pressure_dofs(i-1);}
    int num_vel_dofs=block_vel_dofs.Last()*TV3::m;
    num_dofs=DOF_ID(Value(Number_Particles())+num_vel_dofs);

    for(BLOCK_ID i(0);i<blocks.m;i++) blocks(i).num_dofs=0;
    dof_map.Resize(num_dofs,use_init,{-1,-1});
    pressure_dofs.Resize(Number_Particles(),use_init,DOF_ID(-1));
    for(PARTICLE_ID i(0);i<Number_Particles();i++){
        BLOCK_ID bid=Node_Block(i);
        pressure_dofs(i)=DOF_ID(num_vel_dofs+--block_pressure_dofs(bid));
        dof_map(pressure_dofs(i))={Value(bid),blocks(bid).num_dofs++};}

    vel_node_dofs.Resize(Number_Particles(),use_init,DOF_ID(-3));
    for(PARTICLE_ID i(0);i<Number_Particles();i++){
        BC_ID bc=Node_BC(i);
        if(bc!=BC_ID(-1) && pd.bc(bc).type==dirichlet_v) continue;
        BLOCK_ID bid=Node_Block(i);
        vel_node_dofs(i)=DOF_ID(--block_vel_dofs(bid)*TV3::m);
        for(int a=0;a<TV3::m;a++)
            dof_map(vel_node_dofs(i)+a)={Value(bid),blocks(bid).num_dofs++};}

    vel_edge_dofs.Resize(Number_Edges(),use_init,DOF_ID(-3));
    for(EDGE_ID i(0);i<Number_Edges();i++){
        BC_ID bc=Edge_BC(i);
        if(bc!=BC_ID(-1) && pd.bc(bc).type==dirichlet_v) continue;
        VECTOR<PARTICLE_ID,2> edge=Edge(i);
        BLOCK_ID bid=Edge_Block(edge(0),edge(1));
        vel_edge_dofs(i)=DOF_ID(--block_vel_dofs(bid)*TV3::m);
        for(int a=0;a<TV3::m;a++)
            dof_map(vel_edge_dofs(i)+a)={Value(bid),blocks(bid).num_dofs++};}
}
//#####################################################################
// Function Compute
//#####################################################################
template<typename T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Compute(const PARSE_DATA_FEM<TV,TV3>& pd)
{
    typename FLUID_LAYOUT_FEM<TV>::CONNECTION con;
    for(VERTEX_ID i(0);i<pd.pts.m;i++){
        fl2.Generate_Joint(i,pd,con);}
    fl2.pipes.Resize(pd.pipes.m);
    for(PIPE_ID i(0);i<pd.pipes.m;i++){
        fl2.Generate_Pipe(i,pd,con);}
    fl2.area_hidden.mesh.Set_Number_Nodes(fl2.area_hidden.particles.number);
    fl2.area_hidden.mesh.Initialize_Incident_Elements();
    fl2.area_hidden.mesh.Initialize_Element_Edges();
    fl2.area_hidden.mesh.Initialize_Edge_Triangles();
    fl2.area_hidden.mesh.segment_mesh->Initialize_Incident_Elements();
    Assign_Edge_Blocks_Planar();
    Update_BC_Planar(pd);

    layers=pd.height;
    wall_bc=pd.wall_bc;
    Extrude(pd.height,pd.unit_length);
    //vol_hidden.Initialize_Swept_Mesh_And_Particles(fl2.area_hidden,pd.height,FRAME<TV3>(),pd.unit_length*extrude_dir*pd.height);
    vol_hidden.mesh.Initialize_Incident_Elements();
    vol_hidden.mesh.Initialize_Element_Edges();
    Allocate_Dofs(pd);
}
//#####################################################################
// Function Dump_Input
//#####################################################################
template<typename T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Dump_Input(const PARSE_DATA_FEM<TV,TV3>& pd) const
{
    for(auto& i:pd.pts)
        Add_Debug_Particle(TV3(i.pt),VECTOR<T,3>(0,0,1));
    for(auto& i:pd.pipes)
        Add_Debug_Object<TV3,2>({TV3(pd.pts(i.x).pt),TV3(pd.pts(i.y).pt)},VECTOR<T,3>(0,0,1));
}
//#####################################################################
// Function Dump_Mesh
//#####################################################################
template<typename T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Dump_Mesh() const
{
    for(TRIANGLE_ID i(0);i<Number_Tetrahedrons();i++){
        auto tet=TETRAHEDRON<T>(vol_hidden.particles.X.Subset(Tetrahedron(i)));
        for(int j=0;j<4;j++)
            Add_Debug_Object<TV3,3>(tet.triangle(j).X,VECTOR<T,3>(0.5,0.5,0.5));}
}
//#####################################################################
// Function Dump_Layout
//#####################################################################
template<typename T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Dump_Layout() const
{
    for(BLOCK_ID b(0);b<blocks.m;b++){
        for(TRIANGLE_ID e(0);e<Number_Tetrahedrons();e++){
            if(elem_blocks(e)!=b) continue;
            auto tet=TETRAHEDRON<T>(vol_hidden.particles.X.Subset(Tetrahedron(e)));
            for(int j=0;j<4;j++)
                Add_Debug_Object<TV3,3>(tet.triangle(j).X,VECTOR<T,3>(0.5,0.5,0.5));
            std::string s=LOG::sprintf("%i",elem_blocks(e));
            Add_Debug_Text(tet.X.Average(),s,VECTOR<T,3>(1,1,1));}
        Flush_Frame<TV3>("block");}
}
//#####################################################################
// Function Dump_Node_Blocks
//#####################################################################
template<typename T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Dump_Node_Blocks() const
{
    for(BLOCK_ID b(0);b<blocks.m;b++){
        for(TRIANGLE_ID e(0);e<Number_Tetrahedrons();e++){
            if(elem_blocks(e)!=b) continue;
            auto indices=Tetrahedron(e);
            auto tet=vol_hidden.particles.X.Subset(indices);
            auto tet_tri=TETRAHEDRON<T>(tet);
            for(int j=0;j<4;j++){
                Add_Debug_Object<TV3,3>(tet_tri.triangle(j).X,VECTOR<T,3>(0.5,0.5,0.5));
                Add_Debug_Particle(tet(j),VECTOR<T,3>(1,1,1));
                std::string s=LOG::sprintf("%i",Node_Block(PARTICLE_ID(indices(j))));
                Add_Debug_Text(tet(j),s,VECTOR<T,3>(1,1,1));}}
        Flush_Frame<TV3>("node block");}
}
//#####################################################################
// Function Dump_Edge_Blocks
//#####################################################################
template<typename T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Dump_Edge_Blocks() const
{
    for(BLOCK_ID b(0);b<blocks.m;b++){
        for(TRIANGLE_ID e(0);e<Number_Tetrahedrons();e++){
            if(elem_blocks(e)!=b) continue;
            auto indices=Tetrahedron(e);
            auto tet=vol_hidden.particles.X.Subset(indices);
            for(int j=0;j<4;j++) for(int k=j+1;k<j+4;k++){
                Add_Debug_Object<TV3,2>(
                    {tet(j),tet(k%4)},VECTOR<T,3>(0.5,0.5,0.5));
                std::string s=LOG::sprintf("%i",Edge_Block(PARTICLE_ID(indices(j)),PARTICLE_ID(indices(k%4))));
                Add_Debug_Text(0.5*(tet(j)+tet(k%4)),s,VECTOR<T,3>(1,1,1));}}
        Flush_Frame<TV3>("edge block");}
}
//#####################################################################
// Function Dump_Dofs
//#####################################################################
template<typename T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Dump_Dofs() const
{
    for(BLOCK_ID b(0);b<blocks.m;b++){
        for(TRIANGLE_ID e(0);e<Number_Tetrahedrons();e++){
            if(elem_blocks(e)!=b) continue;
            auto indices=Tetrahedron(e);
            auto tet=vol_hidden.particles.X.Subset(indices);
            auto tet_tri=TETRAHEDRON<T>(tet);
            for(int j=0;j<4;j++){
                Add_Debug_Object<TV3,3>(tet_tri.triangle(j).X,VECTOR<T,3>(0.5,0.5,0.5));
                Add_Debug_Particle(tet(j),VECTOR<T,3>(1,1,1));
                DOF_ID dof=vel_node_dofs(PARTICLE_ID(indices(j)));
                if(dof<DOF_ID(0)) continue;
                std::string s=LOG::sprintf("%i",dof);
                Add_Debug_Text(tet(j),s,VECTOR<T,3>(1,1,1));}}
        Flush_Frame<TV3>("vel node dof");}

    for(BLOCK_ID b(0);b<blocks.m;b++){
        for(TRIANGLE_ID e(0);e<Number_Tetrahedrons();e++){
            if(elem_blocks(e)!=b) continue;
            auto indices=Tetrahedron(e);
            auto tet=vol_hidden.particles.X.Subset(indices);
            auto tet_tri=TETRAHEDRON<T>(tet);
            for(int j=0;j<4;j++){
                Add_Debug_Object<TV3,3>(tet_tri.triangle(j).X,VECTOR<T,3>(0.5,0.5,0.5));
                Add_Debug_Particle(tet(j),VECTOR<T,3>(1,1,1));
                DOF_ID dof=pressure_dofs(PARTICLE_ID(indices(j)));
                if(dof<DOF_ID(0)) continue;
                std::string s=LOG::sprintf("%i",dof);
                Add_Debug_Text(tet(j),s,VECTOR<T,3>(1,1,1));}}
        Flush_Frame<TV3>("pressure dof");}

    vol_hidden.mesh.segment_mesh->Initialize_Incident_Elements();
    for(BLOCK_ID b(0);b<blocks.m;b++){
        for(TRIANGLE_ID e(0);e<Number_Tetrahedrons();e++){
            if(elem_blocks(e)!=b) continue;
            auto indices=Tetrahedron(e);
            auto tet=vol_hidden.particles.X.Subset(indices);
            for(int j=0;j<4;j++) for(int k=j+1;k<j+4;k++){
                Add_Debug_Object<TV3,2>(
                    {tet(j),tet(k%4)},VECTOR<T,3>(0.5,0.5,0.5));
                int p0=indices(j),p1=indices(k%4);
                DOF_ID dof=vel_edge_dofs(EDGE_ID(vol_hidden.mesh.segment_mesh->Segment(p0,p1)));
                if(dof<DOF_ID(0)) continue;
                std::string s=LOG::sprintf("%i",dof);
                Add_Debug_Text(0.5*(tet(j)+tet(k%4)),s,VECTOR<T,3>(1,1,1));}}
        Flush_Frame<TV3>("vel edge dof");}

    Dump_Blocks();
}
//#####################################################################
// Function Dump_Blocks
//#####################################################################
template<typename T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Dump_Blocks() const
{
    ARRAY<ARRAY<int>,BLOCK_ID> block_dofs(blocks.m);
    for(BLOCK_ID b(0);b<blocks.m;b++)
        block_dofs(b).Resize(blocks(b).num_dofs,use_init,-1);
    for(PARTICLE_ID i(0);i<pressure_dofs.m;i++){
        DOF_ID dof=pressure_dofs(i);
        VECTOR<int,2> b=dof_map(dof);
        block_dofs(BLOCK_ID(b(0)))(b(1))=Value(i);
        dof=vel_node_dofs(i);
        if(dof<DOF_ID(0)) continue;
        b=dof_map(dof);
        for(int a=0;a<TV3::m;a++)
            block_dofs(BLOCK_ID(b(0)))(b(1)+a)=Value(i);}
    for(int i=0;i<vol_hidden.mesh.segment_mesh->elements.m;i++){
        VECTOR<int,2> p=vol_hidden.mesh.segment_mesh->elements(i).Sorted();
        DOF_ID dof=vel_edge_dofs(EDGE_ID(i));
        if(dof<DOF_ID(0)) continue;
        VECTOR<int,2> b=dof_map(dof);
        block_dofs(BLOCK_ID(b(0)))(b(1))=p(0);
        block_dofs(BLOCK_ID(b(0)))(b(1)+1)=p(1);}

    for(BLOCK_ID b(0);b<blocks.m;b++){
        int m=Value(Number_Particles());
        for(int i=0;i<block_dofs(b).m;i++)
            if(block_dofs(b)(i)!=-1)
                m=std::min(m,block_dofs(b)(i));
        for(int i=0;i<block_dofs(b).m;i++)
            if(block_dofs(b)(i)!=-1)
                block_dofs(b)(i)-=m;}
    ARRAY<bool,BLOCK_ID> visited(blocks.m,use_init,false);
    for(BLOCK_ID b(0);b<blocks.m;b++){
        if(visited(b)) continue;
        for(BLOCK_ID b1(Value(b)+1);b1<blocks.m;b1++){
            if(visited(b1)) continue;
            for(TRIANGLE_ID e(0);e<Number_Tetrahedrons();e++){
                if(elem_blocks(e)!=b && elem_blocks(e)!=b1) continue;
                if(elem_blocks(e)==b1 && block_dofs(b1)!=block_dofs(b)) continue;
                auto indices=Tetrahedron(e);
                auto tet=vol_hidden.particles.X.Subset(indices);
                auto tet_tri=TETRAHEDRON<T>(tet);
                for(int j=0;j<4;j++){
                    Add_Debug_Object<TV3,3>(tet_tri.triangle(j).X,VECTOR<T,3>(1,1,1));}
                if(elem_blocks(e)==b1) visited(b1)=true;}}
        visited(b)=true;
        Flush_Frame<TV3>("check regular blocks");}
}
//#####################################################################
// Function Print_Statistics
//#####################################################################
template<typename T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Print_Statistics() const
{
}

template<typename T>
const VECTOR<T,3> FLUID_LAYOUT_FEM<VECTOR<T,3> >::extrude_dir(0,0,1);

template struct FLUID_LAYOUT_FEM<VECTOR<double,3> >;
}
