//#####################################################################
// Copyright 2007-2008, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include "PHYSBAM_INTERFACE.h"
using namespace PhysBAM;

int main(int argc,char* argv[])
{
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    typedef double T;
#else
    typedef float T;
#endif
    typedef VECTOR<T,3> TV;
    typedef double RW;

    PARSE_ARGS parse_args;
    parse_args.Add_Double_Argument("-tolerance",(T)0);
    parse_args.Parse(argc,argv);
    T tolerance=(T)parse_args.Get_Double_Value("-tolerance");

//##########################  INITIALIZATION  #########################
    GEOMETRY_PARTICLES<TV> particles;
    particles.array_collection->Resize(8);
    particles.X(1)=TV(-40,-1,-1);particles.X(2)=TV(-40,-1,1);
    particles.X(3)=TV(40,-1,-1);particles.X(4)=TV(40,-1,1);
    particles.X(5)=TV(-40,1,-1);particles.X(6)=TV(-40,1,1);
    particles.X(7)=TV(40,1,-1);particles.X(8)=TV(40,1,1);

    ARRAY<VECTOR<int,3>,int> triangle_list;
    triangle_list.Append(VECTOR<int,3>(0,1,2));triangle_list.Append(VECTOR<int,3>(1,3,2));
    triangle_list.Append(VECTOR<int,3>(4,6,5));triangle_list.Append(VECTOR<int,3>(5,6,7));
    triangle_list.Append(VECTOR<int,3>(0,4,1));triangle_list.Append(VECTOR<int,3>(1,4,5));
    triangle_list.Append(VECTOR<int,3>(2,3,6));triangle_list.Append(VECTOR<int,3>(3,7,6));
    triangle_list.Append(VECTOR<int,3>(0,2,4));triangle_list.Append(VECTOR<int,3>(2,6,4));
    triangle_list.Append(VECTOR<int,3>(1,5,3));triangle_list.Append(VECTOR<int,3>(3,5,7));

    TRIANGLE_MESH mesh(8,triangle_list);mesh.Initialize_Adjacent_Elements();
    TRIANGULATED_SURFACE<T> triangulated_surface(mesh,particles);triangulated_surface.Update_Triangle_List();

    PhysBAMInterface<T> interface(triangulated_surface);interface.Update(true);
//#####################################################################

    ARRAY<TV> nodes_for_rays;nodes_for_rays.Append(triangulated_surface.Get_Element(0).Center());
    nodes_for_rays.Append(triangulated_surface.Get_Element(0).Center()+TV(0,10,0));
    nodes_for_rays.Append(triangulated_surface.Get_Element(0).Center()+TV(0,-10,0));

    ARRAY<PAIR<VECTOR<int,2>,IntersectionResult<T> > > rays_to_intersect(3);
    rays_to_intersect(0).x=VECTOR<int,2>(0,1);
    rays_to_intersect(1).x=VECTOR<int,2>(1,0);
    rays_to_intersect(2).x=VECTOR<int,2>(0,2);
    rays_to_intersect(2).x=VECTOR<int,2>(2,0);

    interface.Intersect(nodes_for_rays,rays_to_intersect,tolerance);
    for(int i=0;i<rays_to_intersect.Size();i++)
        LOG::cout<<rays_to_intersect(i).y.triangleID<<": "
            <<10*rays_to_intersect(i).y.alpha<<" ("<<rays_to_intersect(i).y.zeta[0]<<", "
            <<rays_to_intersect(i).y.zeta[1]<<", "<<rays_to_intersect(i).y.zeta[2]<<")"<<std::endl;

    return 0;
}
//#####################################################################
