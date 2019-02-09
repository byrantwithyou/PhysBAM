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
// Function Allocate_Dofs
//#####################################################################
template<class T> void FLUID_LAYOUT_FEM<VECTOR<T,3> >::
Allocate_Dofs(const PARSE_DATA_FEM<TV,TV3>& pd)
{
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
    Extrude(pd.height,pd.unit_length);
    //vol_hidden.Initialize_Swept_Mesh_And_Particles(fl2.area_hidden,pd.height,FRAME<TV3>(),pd.unit_length*extrude_dir*pd.height);
    vol_hidden.mesh.Initialize_Incident_Elements();
    vol_hidden.mesh.Initialize_Element_Edges();

    fl2.edge_blocks.Resize(fl2.Number_Edges(),init_all,BLOCK_ID(-1));
    for(EDGE_ID i(0);i<fl2.Number_Edges();i++)
        if(fl2.Edge_Triangles(i).m==1){
            auto key=fl2.Edge(i).Sorted();
            if(!fl2.bc_map.Contains(key)){
                if(BC_ID* pbc=fl2.particle_bc_map.Get_Pointer(key.x)) *pbc=pd.wall_bc;
                else fl2.particle_bc_map.Set(key.x,pd.wall_bc);
                if(BC_ID* pbc=fl2.particle_bc_map.Get_Pointer(key.y)) *pbc=pd.wall_bc;
                else fl2.particle_bc_map.Set(key.y,pd.wall_bc);
                fl2.bc_map.Set({key.x,key.y},pd.wall_bc);}}
    for(EDGE_ID ei(0);ei<fl2.Number_Edges();ei++){
        BC_ID* bc_idx=fl2.bc_map.Get_Pointer(fl2.Edge(ei).Sorted());
        if(!bc_idx || pd.bc(*bc_idx).type!=dirichlet_v){
            auto neighbor_tris=fl2.Edge_Triangles(ei);
            PHYSBAM_ASSERT(neighbor_tris.m==1 || neighbor_tris.m==2);
            BLOCK_ID bid=fl2.elem_data(neighbor_tris(0)).block_id;
            if(neighbor_tris.m==2 && fl2.elem_data(neighbor_tris(1)).greed>fl2.elem_data(neighbor_tris(0)).greed)
                bid=fl2.elem_data(neighbor_tris(1)).block_id;
            fl2.edge_blocks(ei)=bid;}}
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
    for(BLOCK_ID b(0);b<fl2.blocks.m;b++){
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
    for(BLOCK_ID b(0);b<fl2.blocks.m;b++){
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
    for(BLOCK_ID b(0);b<fl2.blocks.m;b++){
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
