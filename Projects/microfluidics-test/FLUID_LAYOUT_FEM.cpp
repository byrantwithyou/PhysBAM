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
Generate_Pipe(const TV& v0,const TV& v1,int half_width,T unit_length)
{
    typedef VECTOR<int,3> E;
    GEOMETRY_PARTICLES<TV>& particles=area->particles;
    int base=particles.number;
    //el_base=area->mesh.elements.m;
    TV d=v1-v0;
    T l=d.Normalize();
    TV t(-d.y,d.x);
    int width=2*half_width+1;
    int height=(int)(l/unit_length)+1;
    T avg_step=0;
    T rem=l-(height-1)*unit_length;
    if(rem>1e-2*unit_length){
        avg_step=(rem+unit_length)*0.5;
        height++;}
    particles.Add_Elements(width*height);
    auto pid=[base,width,half_width](int i,int j){return base+i*width+j+half_width;};
    TV next;
    int block=-1;
    bool regular=true;
    for(int i=0;i<height;i++){
        for(int j=-half_width;j<=half_width;j++){
            particles.X(pid(i,j))=v0+next+j*t;
            if(i==0) continue;
            if(j!=half_width){
                area->mesh.elements.Append(E(pid(i,j),pid(i-1,j+1),pid(i-1,j)));
                blocks.Append({block,regular});}
            if(j!=-half_width){
                area->mesh.elements.Append(E(pid(i,j),pid(i-1,j),pid(i,j-1)));
                blocks.Append({block,regular});}}
        if(avg_step && i==height-3)
            regular=false;
        if(regular) next=(i+1)*unit_length*d;
        else next=next+avg_step*d;
        block++;}
    PAIR<ARRAY<int>,ARRAY<int> > f;
    for(int j=-half_width;j<=half_width;j++){
        f.x.Append(pid(0,j));
        f.y.Append(pid(height-1,j));}
    return f;
}
//#####################################################################
// Function Mark_BC
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Mark_BC(const ARRAY<int>& pindices,BC_TYPE bc_type)
{
    for(int i=1;i<pindices.m;i++){
        bc.Set(PAIR<int,int>(pindices(i),pindices(i-1)),{bc_type});
        bc.Set(PAIR<int,int>(pindices(i-1),pindices(i)),{bc_type});}
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void FLUID_LAYOUT_FEM<TV>::
Compute(const PARSE_DATA_FEM<TV>& pd)
{
    for(auto& i:pd.pipes){
        auto f=Generate_Pipe(pd.pts(i.x).pt,pd.pts(i.y).pt,pd.half_width,pd.unit_length);
        if(pd.pts(i.x).bc_type!=nobc)
            Mark_BC(f.x,pd.pts(i.x).bc_type);
        if(pd.pts(i.y).bc_type!=nobc)
            Mark_BC(f.y,pd.pts(i.y).bc_type);}
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
            if(const BC_DATA* bc_data=bc.Get_Pointer(PAIR<int,int>(v0,v1))){
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
