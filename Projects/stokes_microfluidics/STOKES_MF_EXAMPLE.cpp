//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/RANGE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include "STOKES_MF_EXAMPLE.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> STOKES_MF_EXAMPLE<TV>::
STOKES_MF_EXAMPLE(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
    :stream_type(stream_type),
    output_directory("output"),data_directory("../../Public_Data"),
    debug_particles(new DEBUG_PARTICLES<TV>)
{
    parse_args.Extra(&test_number,"example number","example number to run");
    parse_args.Add("-resolution",&resolution,"resolution","grid resolution");
    parse_args.Add("-mu",&viscosity,"mu","viscosity");
    parse_args.Parse();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STOKES_MF_EXAMPLE<TV>::
~STOKES_MF_EXAMPLE()
{
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void STOKES_MF_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    Write_To_File(stream_type,output_directory+"/common/grid",grid);
    debug_particles->Write_Debug_Particles(stream_type,output_directory,frame);
#if 0
    Write_To_File(stream_type,LOG::sprintf("%s/%d/mac_velocities",output_directory.c_str(),frame),phases(PHASE_ID()).velocity);
#endif
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void STOKES_MF_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
}
//#####################################################################
// Function Add_Vertex
//#####################################################################
template<class TV> int STOKES_MF_EXAMPLE<TV>::
Add_Vertex(const TV_INT& x)
{
    iverts.Append(x);
    return iverts.m-1;
}
//#####################################################################
// Function Add_Edge
//#####################################################################
template<class TV> int STOKES_MF_EXAMPLE<TV>::
Add_Edge(int v,int dir,int len)
{
    iverts.Append(iverts(v));
    iverts.Last()(dir)+=len;
    edges.Append(PAIR<int,int>(v,iverts.m-1));
    return iverts.m-1;
}
//#####################################################################
// Function Add_Edge
//#####################################################################
template<class TV> void STOKES_MF_EXAMPLE<TV>::
Add_Edge(int v0,int v1)
{
    edges.Append(PAIR<int,int>(v0,v1));
}
//#####################################################################
// Function Build_Grid
//#####################################################################
static int gcd(int a,int b)
{
    if(b==0) return a;
    else return gcd(b,a%b);
}
template<class TV> void STOKES_MF_EXAMPLE<TV>::
Build_Grid(T len_scale,int cross_section_radius)
{
    TV_INT min_corner(iverts(0)),max_corner(iverts(0));
    for(int i=1;i<iverts.m;i++){
        min_corner=TV_INT::Componentwise_Min(min_corner,iverts(i));
        max_corner=TV_INT::Componentwise_Max(max_corner,iverts(i));}
    TV_INT extent=max_corner-min_corner;
    int d=gcd(cross_section_radius,extent(0));
    for(int i=1;i<TV_INT::m;i++) d=gcd(d,extent(i));
    TV_INT aspect=extent/=d;
    RANGE<TV> domain=RANGE<TV>(TV(),TV(extent)*len_scale).Thickened(cross_section_radius*len_scale);
    grid.Initialize((aspect+2)*resolution,domain,true);
    for(int i=0;i<iverts.m;i++){
        vertices.Append(resolution*((iverts(i)-min_corner)/d+1));
        Add_Debug_Particle(grid.Node(vertices(i)),VECTOR<T,3>(1,0,0));}

    domain_cells.Remove_All();
    boundary_faces.Remove_All();
    ARRAY<VECTOR<bool,2*TV::m> > vert_connect(vertices.m); // -x +x -y +y -z +z
    for(int i=0;i<edges.m;i++){
        Add_Debug_Object(VECTOR<TV,2>(grid.Node(vertices(edges(i).x)),grid.Node(vertices(edges(i).y))),VECTOR<T,3>(0,1,0));
        int v0=edges(i).x,v1=edges(i).y;
        TV_INT first=vertices(v0),second=vertices(v1);
        int d;
        for(d=0;d<TV::m;d++)
            if(first(d)!=second(d)){
                if(first(d)>second(d)){
                    std::swap(first,second);
                    std::swap(v0,v1);}
                first(d)+=resolution;
                second(d)-=resolution;
                break;}
        for(int j=0;j<TV::m;j++)
            if(d!=j){
                first(j)-=resolution;
                second(j)+=resolution;}
        vert_connect(v0)(2*d+1)=true;
        vert_connect(v1)(2*d)=true;
        RANGE<TV_INT> tube(first,second);
        for(RANGE_ITERATOR<TV::m> it(tube);it.Valid();it.Next()){
            bool strict_inside=it.index.All_Greater(first) && it.index.All_Less(second-1);
            if(tube.Lazy_Inside_Half_Open(it.index) && !strict_inside){
                for(int j=0;j<TV::m;j++){
                    if(j==d) continue;
                    TV_INT cell=it.index;
                    cell(j)+=1;
                    if(tube.Lazy_Outside_Half_Open(cell))
                        boundary_faces.Insert(FACE_INDEX<TV::m>(j,cell));
                    cell(j)-=2;
                    if(tube.Lazy_Outside_Half_Open(cell))
                        boundary_faces.Insert(FACE_INDEX<TV::m>(j,it.index));}}
            domain_cells.Append(it.index);}}
    for(int i=0;i<vertices.m;i++){
        TV_INT v=vertices(i);
        RANGE<TV_INT> cube(v-resolution,v+resolution);
        for(RANGE_ITERATOR<TV::m> it(cube);it.Valid();it.Next()){
            bool strict_inside=it.index.All_Greater(v-resolution) && it.index.All_Less(v+resolution-1);
            if(cube.Lazy_Inside_Half_Open(it.index) && !strict_inside){
                for(int j=0;j<TV::m;j++){
                    TV_INT cell=it.index;
                    cell(j)+=1;
                    if(cube.Lazy_Outside_Half_Open(cell) && !vert_connect(i)(2*j+1))
                        boundary_faces.Insert(FACE_INDEX<TV::m>(j,cell));
                    cell(j)-=2;
                    if(cube.Lazy_Outside_Half_Open(cell) && !vert_connect(i)(2*j))
                        boundary_faces.Insert(FACE_INDEX<TV::m>(j,it.index));}}
            domain_cells.Append(it.index);}}
    for(typename HASHTABLE<FACE_INDEX<TV::m> >::ITERATOR it(boundary_faces);it.Valid();it.Next()){
        Add_Debug_Particle(grid.Face(it.Key()),VECTOR<T,3>(1,1,0));}
}
//#####################################################################
namespace PhysBAM{
template class STOKES_MF_EXAMPLE<VECTOR<float,2> >;
template class STOKES_MF_EXAMPLE<VECTOR<float,3> >;
template class STOKES_MF_EXAMPLE<VECTOR<double,2> >;
template class STOKES_MF_EXAMPLE<VECTOR<double,3> >;
}
