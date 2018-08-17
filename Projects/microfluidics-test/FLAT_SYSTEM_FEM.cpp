//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Log/LOG.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "FLAT_SYSTEM_FEM.h"
#include "FLUID_LAYOUT_FEM.h"
#if 0
namespace PhysBAM{
template<class T,class TV>
void Compute_Full_Matrix(ARRAY<VECTOR<int,3> >& coded_entries,
    ARRAY<T>& code_values,ARRAY<T>& rhs_vector,FLUID_LAYOUT_FEM<TV>& fl,T mu,T unit_length)
{
    T lcd=6;
    // dim * 6 * 6
    int dd[]={3,1,0,0,0,-4,1,3,0,0,0,-4,0,0,0,0,0,0,0,0,0,8,
        -8,0,0,0,0,-8,8,0,-4,-4,0,0,0,8,3,0,1,0,-4,0,0,
        0,0,0,0,0,1,0,3,0,-4,0,0,0,0,8,0,-8,-4,0,-4,0,8,
        0,0,0,0,-8,0,8};
    // dim * 3 * 6
    int pd[]={-1,0,0,1,-1,1,0,1,0,1,-1,-1,0,0,0,2,-2,0,-1,0,0,
        1,1,-1,0,0,0,2,0,-2,0,0,1,1,-1,-1};

    auto dvdv=[&dd](int r,int a,int b)
    {
        int idx=b+a*6+r*6*6;
        return dd[idx];
    };
    auto pdv=[&pd](int r,int a,int b)
    {
        int idx=b+a*6+r*3*6;
        return pd[idx];
    };

    rhs_vector.Resize(fl.num_vel_dofs+fl.num_pressure_dofs,init_all,0);
    fl.area.mesh.Initialize_Incident_Elements();
    fl.area.mesh.Initialize_Edge_Triangles();
    fl.area.mesh.Initialize_Element_Edges();
    ARRAY<ARRAY<int> >& incident_elements=*fl.area.mesh.incident_elements;
    ARRAY<ARRAY<int> >& edge_triangles=*fl.area.mesh.edge_triangles;
    ARRAY<VECTOR<int,3> >& element_edges=*fl.area.mesh.element_edges;
    SEGMENT_MESH& segment_mesh=*fl.area.mesh.segment_mesh;
    auto node2local=[&fl](int elem,int g)
    {
        for(int i=0;i<3;i++) if(fl.area.mesh.elements(elem)(i)==g) return i;
        PHYSBAM_ASSERT(false);
    };
    auto edge2local=[&fl,&element_edges,&segment_mesh](int elem,int g)
    {
        int p0,p1;
        for(int i=0;i<3;i++) if(element_edges(elem)(i)==g){
            p0=segment_mesh.elements(g)(0);
            p1=segment_mesh.elements(g)(1);
            break;}
        for(int i=0;i<3;i++){
            int i1=(i+1)%3;
            int q0=fl.area.mesh.elements(elem)(i),q1=fl.area.mesh.elements(elem)(i1);
            if((p0==q0 && p1==q1) || (p0==q1 && p1==q0))
                return (i1+1)%3+3;}
    };
    auto nodes2edge=[&element_edegs,&segment_mesh](int elem,int p0,int p1)
    {
        for(int ei:element_edges(elem)){
            int q0=segment_mesh.elements(ei)(0),q1=segment_mesh.elements(ei)(1);
            if((q0==p0 && q1==p1) || (q0==p1 && q1==p0))
                return ei;}
        PHYSBAM_ASSERT(false);
    };
    auto add_entry=[&coded_entries,&code_values](int r,int c,T value)
    {
        int idx=code_values.Append(value);
        coded_entries.Append({r,c,idx});
    };

    for(int i=0;i<fl.area.particles.number;i++){
        int b=fl.vel_node_dofs(i),pb=fl.pressure_dofs(i);
        if(pb<0) continue;
        for(int axis=0;axis<TV::m;axis++){
            T intdvdv=0,intpdv=0;
            for(int inc_tri:incident_elements(i)){
                int l=node2local(inc_tri,i);
                T det=2*fl.area.Area(inc_tri);
                intdvdv+=2*mu*dvdv(axis,l,l)/lcd*det;
                intpdv+=pdv(axis,l,l)/lcd*det;}
            add_entry(b+axis,b+axis,intdvdv);
            add_entry(b+axis,pb,-intpdv);
            add_entry(pb,b+axis,-intpdv);}}

    for(int i=0;i<fl.area.mesh.elements.m;i++){
        for(int j=0;j<3;j++){
            int p0=fl.area.mesh.elements(i)(j),p1=fl.area.mesh.elements(i)((j+1)%3);
            int ei=nodes2edge(i,p0,p1);
            int dof_vedge=fl.vel_edge_dofs(ei);
            if(dof_vedge<0) continue;
            for(int axis=0;axis<TV::m;axis++){
                T intdvdv=0,intdvdv_half=0;
                T intpdv=0,intp0dv_half=0;
                for(int neighbor_tri:edge_triangles(ei)){
                    int n0=node2local(neighbor_tri,p0),n1=node2local(neighbor_tri,p1);
                    int edge=edge2local(neighbor_tri,ei);
                    T det=2*fl.area.Area(neighbor_tri);
                    intdvdv+=2*dvdv(axis,n0,n1)/lcd*det;
                    intdvdv_half+=2*mu*dvdv(axis,n0,edge)/lcd*det;
                    intpdv+=pdv(axis,n0,n1)/lcd*det;
                    intp0dv_half+=pdv(axis,n0,edge)/lcd*det;}}
            int dof_v0=fl.vel_node_dofs(p0),dof_v1=fl.vel_node_dofs(p1);
            int dof_p0=fl.pressure_dofs(p0),dof_p1=fl.pressure_dofs(p1);
            if(dof_v0>=0 && dof_v1>=0){
                add_entry(dof_v0+axis,dof_v1+axis,intdvdv);
                add_entry(dof_p0,dof_v1+axis,-intpdv);}
            if(dof_v0>=0){
                add_entry(dof_vedge+axis,dof_v0+axis,intdvdv_half);
                add_entry(dof_v0+axis,dof_vedge+axis,intdvdv_half);
                add_entry(dof_p0,dof_vedge+axis,-intp0dv_half);
                add_entry(dof_vedge+axis,dof_p0,-intp0dv_half);}}}

    for(int i=0;i<segment_mesh.elements.m;i++){
        int dof_vedge=fl.vel_edge_dofs(i);
        if(dof_vedge<0) continue;
        int p0=segment_mesh.elements(i)(0),p1=segment_mesh.elements(i)(1);
        for(int axis=0;axis<TV::m;axis++){
            T intdvdv=0,intdvdv_self=0;
            T intp0dv=0,intp1dv=0;
            for(int neighbor_tri:edge_triangles(i)){
                int n0=node2local(neighbor_tri,p0),n1=node2local(neighbor_tri,p1);
                int edge=edge2local(neighbor_tri,i);
                T det=2*fl.area.Area(neighbor_tri);
                intdvdv+=2*dvdv(axis,n0,n1)/lcd*det;
                intdvdv_self+=2*mu*dvdv(axis,edge,edge)/lcd*det;
                intp0dv+=pdv(axis,n0,edge)/lcd*det;
                intp1dv+=pdv(axis,n1,edge)/lcd*det;}
            int dof_v0=fl.vel_node_dofs(p0),dof_v1=fl.vel_node_dofs(p1);
            int dof_p0=fl.pressure_dofs(p0),dof_p1=fl.pressure_dofs(p1);
            add_entry(dof_vedge+axis,dof_vedge+axis,intdvdv_self);
            if(dof_v0<0 || dof_v1<0){
                int vdof=dof_v0<0?dof_v1:dof_v0;
                int pdof=dof_p0<0?dof_p1:dof_p0;
                int dof_node=dof_v0<0?p1:p0;
                int boundary_node=dof_v0<0?p0:p1;
                TV bc=fl.bc(fl.particle_bc_map.Get(boundary_node)).bc;
                rhs(dof+axis)-=bc(axis)*intdvdv;
                rhs(pdof)+=bc(axis)*intpdv;
                rhs(dof_vedge+axis)-=bc(axis)*intdvdv_half;

                if(boundary_node==p1){
                    rhs(dof_p0)+=bc(axis)*intpdv;
                    add_entry(dof_p0,dof_vedge+axis,-intp0dv);
                    add_entry(dof_vedge+axis,dof_p0,-intp0dv);}
                else{
                    rhs(dof_p1)+=bc(axis)*intpdv;
                    add_entry(dof_p1,dof_vedge+axis,-intp1dv);
                    add_entry(dof_vedge+axis,dof_p1,-intp1dv);}
            }
        }
    }
}

template void Compute_Full_Matrix<double,VECTOR<double,2> >(ARRAY<VECTOR<int,3>,int>&,
    ARRAY<double>&,ARRAY<double,int>&,FLUID_LAYOUT_FEM<VECTOR<double,2> >&,double,double);
}
#endif
