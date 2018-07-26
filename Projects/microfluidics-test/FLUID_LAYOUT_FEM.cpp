//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <map>
#include "FLUID_LAYOUT_FEM.h"

namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_LAYOUT_FEM<TV>::
FLUID_LAYOUT_FEM(): area(TRIANGULATED_AREA<T>::Create())
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUID_LAYOUT_FEM<TV>::
~FLUID_LAYOUT_FEM()
{
}
//#####################################################################
// Function Generate_Pipe
//#####################################################################
template<class TV> PAIR<ARRAY<int>,ARRAY<int> > FLUID_LAYOUT_FEM<TV>::
Generate_Pipe(const TV& v0,const TV& v1,int half_num_cells,T unit_length,
    const HASHTABLE<PAIR<int,int>,int>& shared_point)
{
    typedef VECTOR<int,3> E;
    GEOMETRY_PARTICLES<TV>& particles=area->particles;
    int base=particles.number;
    //el_base=area->mesh.elements.m;
    TV d=v1-v0;
    T l=d.Normalize();
    TV t(-d.y,d.x);
    int width=2*half_num_cells+1;
    int height=(int)(l/unit_length)+1;
    T min_percent=0.5;
    bool irregular_last=false;
    T rem=l-(height-1)*unit_length;
    if(rem>1e-2*unit_length){
        irregular_last=true;
        if(rem>min_percent*unit_length)
            height++;}
    particles.Add_Elements(width*height);
    auto pid=[base,width,half_num_cells,height,shared_point](int i,int j)
    {
        int ri=i==height-1?-1:i;
        if(const int* r=shared_point.Get_Pointer({ri,j}))
            return *r;
        return base+i*width+j+half_num_cells;
    };
    TV next;
    bool regular=true;
    for(int i=0;i<height;i++){
        for(int j=-half_num_cells;j<=half_num_cells;j++){
            particles.X(pid(i,j))=v0+next+j*unit_length*t;
            if(i==0) continue;
            if(j!=half_num_cells){
                area->mesh.elements.Append(E(pid(i,j),pid(i-1,j+1),pid(i-1,j)));
                blocks.Append({last_block_id,regular});}
            if(j!=-half_num_cells){
                area->mesh.elements.Append(E(pid(i,j),pid(i-1,j),pid(i,j-1)));
                blocks.Append({last_block_id,regular});}}
        if(i==height-2){
            if(irregular_last) regular=false;
            next=l*d;}
        else next=(i+1)*unit_length*d;
        if(i!=0) last_block_id++;}
    PAIR<ARRAY<int>,ARRAY<int> > f;
    for(int j=-half_num_cells;j<=half_num_cells;j++){
        f.x.Append(pid(0,j));
        f.y.Append(pid(height-1,j));}
    return f;
}
//#####################################################################
// Function Weld_Arc
//#####################################################################
template<class TV> ARRAY<int> FLUID_LAYOUT_FEM<TV>::
Weld_Arc(int p0,int p1,const ARRAY<int>& side,const TV& c,T unit_length)
{
    typedef VECTOR<int,3> E;
    GEOMETRY_PARTICLES<TV>& particles=area->particles;
    T polyline_arc_ratio=0.999;
    TV v0=particles.X(p0)-c,v1=particles.X(p1)-c;
    TV v=v0.Normalized();
    T r=v0.Magnitude();
    T unit_rad=unit_length/polyline_arc_ratio/r;
    T total_rad=acos(v1.Normalized().Dot(v));
    int num_segs=(int)(total_rad/unit_rad);
    T min_percent=0.5;
    T rem=total_rad-num_segs*unit_rad;
    int base=particles.number;
    if(rem>min_percent*unit_rad)
        num_segs++;
    particles.Add_Elements(num_segs-1);

    ARRAY<int> new_side;
    auto next=[v1,num_segs,p1,unit_rad,v0,&particles,base,c,&new_side](int i,TV& g)
    {
        int next_p;
        if(i==num_segs-1){
            g=v1;
            next_p=p1;}
        else{
            T rad=-(i+1)*unit_rad;
            g=TV(cos(rad)*v0(0)-sin(rad)*v0(1),sin(rad)*v0(0)+cos(rad)*v0(1));
            particles.X(base+i)=c+g;
            next_p=base+i;}
        new_side.Append(next_p);
        return next_p;
    };

    int p=p0,j=0,i=1,pnext;
    new_side.Append(p0);
    TV u;
    pnext=next(0,u);
    while(p!=p1){
        if(j==side.m-1 || u.Normalized().Dot(v)>(particles.X(side(j+1))-c).Normalized().Dot(v)){
            area->mesh.elements.Append(E(pnext,p,side(j)));
            blocks.Append({last_block_id,false});
            p=pnext;
            if(i<num_segs) pnext=next(i,u);
            i++;}
        else{
            area->mesh.elements.Append(E(side(j+1),p,side(j)));
            blocks.Append({last_block_id,false});
            j++;}}
    if(j==side.m-2){
        area->mesh.elements.Append(E(side(j+1),p,side(j)));
        blocks.Append({last_block_id,false});}
    last_block_id++;
    return new_side;
}
//#####################################################################
// Function Mark_BC
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Mark_BC(const ARRAY<int>& pindices,BC_TYPE bc_type)
{
    for(int i=1;i<pindices.m;i++){
        bc.Set({pindices(i),pindices(i-1)},{bc_type});
        bc.Set({pindices(i-1),pindices(i)},{bc_type});}
}
//#####################################################################
// Function Pipe_Joint_Connection
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Pipe_Joint_Connection(const TV& joint,const TV& p0,const TV& p1,int half_width,T unit_length,TV& q0,TV& q1) const
{
    TV e0=(p0-joint).Normalized(),e1=(p1-joint).Normalized();
    TV m=(e0+e1).Normalized();
    T s=half_width/abs(e0.Cross(m)(0));
    q0=joint+s*m.Dot(e0)*e0;
    q1=joint+s*m.Dot(e1)*e1;
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Compute(const PARSE_DATA_FEM<TV>& pd)
{
    GEOMETRY_PARTICLES<TV>& particles=area->particles;
    for(auto& i:pd.pipes){
        if(pd.pts(i.x).joint_type==end_vertex && pd.pts(i.y).joint_type==end_vertex){
            auto f=Generate_Pipe(pd.pts(i.x).pt,pd.pts(i.y).pt,pd.half_num_cells,pd.unit_length);
            if(pd.pts(i.x).bc_type!=nobc)
                Mark_BC(f.x,pd.pts(i.x).bc_type);
            if(pd.pts(i.y).bc_type!=nobc)
                Mark_BC(f.y,pd.pts(i.y).bc_type);}}

    for(int i=0;i<pd.pts.m;i++){
        if(pd.pts(i).joint_type==arc){
            const ARRAY<int>* pipes=pd.joints.Get_Pointer(i);
            PHYSBAM_ASSERT(pipes);
            PHYSBAM_ASSERT(pipes->m==2);
            int p0=(*pipes)(0),p1=(*pipes)(1);
            int j0=pd.pipes(p0)(0),j1=pd.pipes(p1)(0);
            if(j0==i) j0=pd.pipes(p0)(1);
            if(j1==i) j1=pd.pipes(p1)(1);
            if(pd.pts(j0).joint_type!=end_vertex || pd.pts(j1).joint_type!=end_vertex) continue;
            TV joint=pd.pts(i).pt;
            if((pd.pts(j0).pt-joint).Cross(pd.pts(j1).pt-joint)(0)<0){
                std::swap(p0,p1);
                std::swap(j0,j1);}
            TV q0,q1;
            Pipe_Joint_Connection(joint,pd.pts(j0).pt,pd.pts(j1).pt,pd.half_width,pd.unit_length,q0,q1);
            //Add_Debug_Particle(q0,VECTOR<T,3>(1,0,0));
            //Add_Debug_Particle(q1,VECTOR<T,3>(1,0,0));
            auto f0=Generate_Pipe(q0,pd.pts(j0).pt,pd.half_num_cells,pd.unit_length);
            int center_pid=f0.x(2*pd.half_num_cells);
            HASHTABLE<PAIR<int,int>,int> shared_point;
            shared_point.Set({0,-pd.half_num_cells},center_pid);
            auto f1=Generate_Pipe(q1,pd.pts(j1).pt,pd.half_num_cells,pd.unit_length,shared_point);

            TV center=particles.X(center_pid);
            ARRAY<int> side;
            side.Append(center_pid);
            for(int j=1;j<=2*pd.half_num_cells;j++){
                ARRAY<int> new_side=Weld_Arc(f0.x(2*pd.half_num_cells-j),f1.x(j),side,center,pd.unit_length);
                side=new_side;}

            if(pd.pts(j0).bc_type!=nobc)
                Mark_BC(f0.y,pd.pts(j0).bc_type);
            if(pd.pts(j1).bc_type!=nobc)
                Mark_BC(f1.y,pd.pts(j1).bc_type);
        }
    }
    Dump_Input(pd);
}
//#####################################################################
// Function Dump_Mesh
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Dump_Mesh() const
{
    const GEOMETRY_PARTICLES<TV>& particles=area->particles;
    //for(int i=0;i<particles.number;i++)
    //    Add_Debug_Particle(particles.X(i),VECTOR<T,3>(1,1,1));
    VECTOR<T,3> white(1,1,1),dirichlet_bc(1,0,1),traction_bc(0,1,1);
    for(int i=0;i<area->mesh.elements.m;i++){
        //Add_Debug_Object(VECTOR<TV,3>(particles.X.Subset(area->mesh.elements(i))),VECTOR<T,3>(1,1,1));
        for(int j=0;j<3;j++){
            VECTOR<T,3> color=white;
            int v0=area->mesh.elements(i)(j),v1=area->mesh.elements(i)((j+1)%3);
            if(const BC_DATA* bc_data=bc.Get_Pointer({v0,v1})){
                if(bc_data->bc_type==dirichlet_v)
                    color=dirichlet_bc;
                else if(bc_data->bc_type==traction)
                    color=traction_bc;}
            Add_Debug_Object<TV,2>(VECTOR<TV,2>(particles.X(v0),particles.X(v1)),color);}
    }
}
//#####################################################################
// Function Dump_Layout
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Dump_Layout() const
{
    const GEOMETRY_PARTICLES<TV>& particles=area->particles;
    Dump_Mesh();
    VECTOR<T,3> reg(0,1,0),irreg(1,0,0);
    for(int i=0;i<area->mesh.elements.m;i++){
        std::string s=LOG::sprintf("%i",blocks(i).block_id);
        Add_Debug_Text(particles.X.Subset(area->mesh.elements(i)).Average(),s,
            blocks(i).regular?reg:irreg);}
}
//#####################################################################
// Function Dump_Input
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Dump_Input(const PARSE_DATA_FEM<TV>& pd) const
{
    for(auto& i:pd.pts)
        Add_Debug_Particle(i.pt,VECTOR<T,3>(0,0,1));
    for(auto& i:pd.pipes)
        Add_Debug_Object<TV,2>(VECTOR<TV,2>(pd.pts(i.x).pt,pd.pts(i.y).pt),VECTOR<T,3>(0,0,1));
}
template struct FLUID_LAYOUT_FEM<VECTOR<double,2> >;
}
