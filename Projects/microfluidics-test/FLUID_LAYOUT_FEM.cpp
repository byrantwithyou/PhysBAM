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
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_Pipe(int pipe,const PARSE_DATA_FEM<TV>& pd,const CONNECTION& con)
{
    typedef VECTOR<int,3> E;
    GEOMETRY_PARTICLES<TV>& particles=area->particles;
    int base=particles.number;
    const ARRAY<int>* f0=con.Get_Pointer({pd.pipes(pipe).x,pipe}),*f1=con.Get_Pointer({pd.pipes(pipe).y,pipe});
    PHYSBAM_ASSERT(f0 && f1);
    TV v0=particles.X((*f0)(pd.half_width)),v1=particles.X((*f1)(pd.half_width));
    TV d=v1-v0;
    T l=d.Normalize();
    TV t(-d.y,d.x);
    int width=2*pd.half_width+1;
    int height=(int)(l/pd.unit_length)+1;
    T min_percent=0.5;
    bool irregular_last=false;
    T rem=l-(height-1)*pd.unit_length;
    if(rem>1e-2*pd.unit_length){
        irregular_last=true;
        if(rem>min_percent*pd.unit_length)
            height++;}
    if(height>2)
        particles.Add_Elements(width*(height-2));
    auto pid=[base,width,height,&pd,f0,f1](int i,int j)
    {
        if(i==0) return (*f0)(pd.half_width+j);
        else if(i==height-1) return (*f1)(pd.half_width+j);
        else return base+(i-1)*width+j+pd.half_width;
    };
    TV next=pd.unit_length*d;;
    bool regular=true;
    for(int i=1;i<height;i++){
        for(int j=-pd.half_width;j<=pd.half_width;j++){
            if(i!=height-1)
                particles.X(pid(i,j))=v0+next+j*pd.unit_length*t;
            if(j!=pd.half_width){
                area->mesh.elements.Append(E(pid(i,j),pid(i-1,j+1),pid(i-1,j)));
                elem_data.Append({blocks.m});}
            if(j!=-pd.half_width){
                area->mesh.elements.Append(E(pid(i,j),pid(i-1,j),pid(i,j-1)));
                elem_data.Append({blocks.m});}}
        if(i!=0) blocks.Append({regular});
        if(i==height-2){
            if(irregular_last) regular=false;
            next=l*d;}
        else next=(i+1)*pd.unit_length*d;}
}
//#####################################################################
// Function March_Arc
//#####################################################################
template<class TV> ARRAY<int> FLUID_LAYOUT_FEM<TV>::
March_Arc(int p0,int p1,const ARRAY<int>& side,const TV& c,T unit_length)
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
    if(rem>min_percent*unit_rad || num_segs==0)
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
        if(j==side.m-1 || u.Normalized().Dot(v)>=(particles.X(side(j+1))-c).Normalized().Dot(v)){
            area->mesh.elements.Append(E(pnext,p,side(j)));
            p=pnext;
            if(i<num_segs) pnext=next(i,u);
            i++;}
        else{
            area->mesh.elements.Append(E(side(j+1),p,side(j)));
            j++;}
        elem_data.Append({blocks.m});}
    if(j==side.m-2){
        area->mesh.elements.Append(E(side(j+1),p,side(j)));
        elem_data.Append({blocks.m});}
    blocks.Append({false});
    return new_side;
}
//#####################################################################
// Function Mark_BC
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Mark_BC(const ARRAY<int>& pindices,BC_TYPE bc_type)
{
    int bc_index=bc.Append({bc_type});
    for(int i=0;i<pindices.m;i++)
        bc_map.Set(pindices(i),bc_index);
}
//#####################################################################
// Function Generate_End
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_End(int i,int pipe,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con)
{
    int j=pd.pipes(pipe).x;
    if(j==i) j=pd.pipes(pipe).y;
    TV p=pd.pts(i).pt;
    TV v=pd.pts(j).pt-pd.pts(i).pt;
    TV u=TV(-v.y,v.x).Normalized();
    GEOMETRY_PARTICLES<TV>& particles=area->particles;
    int base=particles.Add_Elements(2*pd.half_width+1);
    ARRAY<int>& indices=con.Get_Or_Insert({i,pipe});
    for(int k=-pd.half_width;k<=pd.half_width;k++){
        int pid=base+k+pd.half_width;
        indices.Append(pid);
        particles.X(pid)=p+k*pd.unit_length*u;}
    if(pd.pipes(pipe).x!=i) indices.Reverse();
    if(pd.pts(i).bc_type!=nobc) Mark_BC(indices,pd.pts(i).bc_type);
}
//#####################################################################
// Function Arc
//#####################################################################
template<class TV> VECTOR<TV,3> FLUID_LAYOUT_FEM<TV>::
Arc(const TV& joint,const TV& p0,const TV& p1,int half_width,T unit_length) const
{
    TV e0=(p0-joint).Normalized(),e1=(p1-joint).Normalized();
    TV m=(e0+e1).Normalized();
    T s=half_width*unit_length/abs(e0.Cross(m)(0));
    return VECTOR<TV,3>(joint+s*m,(s*m.Dot(e0)*e0-s*m).Normalized(),(s*m.Dot(e1)*e1-s*m).Normalized());
}
//#####################################################################
// Function Generate_Arc
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_Arc(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con)
{
    GEOMETRY_PARTICLES<TV>& particles=area->particles;
    const ARRAY<int>* pipes=pd.joints.Get_Pointer(i);
    PHYSBAM_ASSERT(pipes);
    PHYSBAM_ASSERT(pipes->m==2);
    int p0=(*pipes)(0),p1=(*pipes)(1);
    int j0=pd.pipes(p0)(0),j1=pd.pipes(p1)(0);
    if(j0==i) j0=pd.pipes(p0)(1);
    if(j1==i) j1=pd.pipes(p1)(1);
    TV joint=pd.pts(i).pt;
    if((pd.pts(j0).pt-joint).Cross(pd.pts(j1).pt-joint)(0)<0){
        std::swap(p0,p1);
        std::swap(j0,j1);}
    VECTOR<TV,3> arc=Arc(joint,pd.pts(j0).pt,pd.pts(j1).pt,pd.half_width,pd.unit_length);
    int center_pid=particles.Add_Element();
    particles.X(center_pid)=arc(0);
    ARRAY<int>& f0=con.Get_Or_Insert({i,p0}),&f1=con.Get_Or_Insert({i,p1});
    f0.Append(center_pid);
    f1.Append(center_pid);
    ARRAY<int> side;
    side.Append(center_pid);
    for(int j=0;j<2*pd.half_width;j++){
        int s0=particles.Add_Element();
        f0.Append(s0);
        particles.X(s0)=arc(0)+arc(1)*pd.unit_length*(j+1);
        int s1=particles.Add_Element();
        f1.Append(s1);
        particles.X(s1)=arc(0)+arc(2)*pd.unit_length*(j+1);
        ARRAY<int> new_side=March_Arc(s0,s1,side,arc(0),pd.unit_length);
        side=new_side;}
    if(pd.pipes(p0).x==i) f0.Reverse();
    if(pd.pipes(p1).x!=i) f1.Reverse();
}
//#####################################################################
// Function Generate_Joint
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Generate_Joint(int i,const PARSE_DATA_FEM<TV>& pd,CONNECTION& con)
{
    if(pd.pts(i).joint_type==end_vertex){
        const ARRAY<int>* p=pd.joints.Get_Pointer(i);
        PHYSBAM_ASSERT(p);
        PHYSBAM_ASSERT(p->m==1);
        Generate_End(i,(*p)(0),pd,con);
        return;}

    if(pd.pts(i).joint_type==arc_joint){
        Generate_Arc(i,pd,con);
        return;}

    PHYSBAM_FATAL_ERROR("joint type not handled");
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Compute(const PARSE_DATA_FEM<TV>& pd)
{
    CONNECTION con;
    for(int i=0;i<pd.pts.m;i++){
        Generate_Joint(i,pd,con);}
    for(int i=0;i<pd.pipes.m;i++){
        Generate_Pipe(i,pd,con);}
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
            if(const int* bc0=bc_map.Get_Pointer(v0)){
                const int* bc1=bc_map.Get_Pointer(v1);
                if(bc1 && *bc1==*bc0)
                    color=bc(*bc0).bc_type==dirichlet_v?dirichlet_bc:traction_bc;}
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
        std::string s=LOG::sprintf("%i",elem_data(i).block_id);
        Add_Debug_Text(particles.X.Subset(area->mesh.elements(i)).Average(),s,
            blocks(elem_data(i).block_id).regular?reg:irreg);}
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
