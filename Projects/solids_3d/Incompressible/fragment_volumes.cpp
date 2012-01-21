//#####################################################################
// Copyright 2007, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Deformable_Objects/DEFORMABLE_OBJECT.h>
using namespace PhysBAM;

typedef float T;
typedef float RW;
typedef VECTOR<T,3> TV;

int main(int argc,char* argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Set_Extra_Arguments(-1,"<frame>");
    parse_args.Parse(argc,argv);

    int start_frame=0,last_frame=0;
    if(parse_args.Num_Extra_Args()>=1) start_frame=last_frame=atoi(parse_args.Extra_Arg(0).c_str());
    else FILE_UTILITIES::Read_From_Text_File("last_frame",last_frame);

    RANGE<VECTOR<T,1> > volume_box=RANGE<VECTOR<T,1> >(FLT_MAX,-FLT_MAX);T rest_volume=(T)4.85889;
    for(int frame=start_frame;frame<=last_frame;frame++){
        LOG::cout<<"frame "<<frame<<": ";
        // read
        std::string f=STRING_UTILITIES::string_sprintf(".%d",frame);
        COLLISION_BODY_LIST<TV> collision_body_list;
        DEFORMABLE_OBJECT<TV> deformable_object(collision_body_list);
        deformable_object.Read(STREAM_TYPE(RW()),"",frame,-1,true);
        TETRAHEDRALIZED_VOLUME<T>& volume=deformable_object.Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        PARTICLES<TV>& particles=volume.particles;
        // find components
        UNION_FIND<> union_find(particles.array_collection->Size());
        for(int t=0;t<volume.mesh.elements.m;t++) union_find.Union(volume.mesh.elements(t));
        // compute volumes
        ARRAY<T> volumes(volume.particles.array_collection->Size());
        for(int t=0;t<volume.mesh.elements.m;t++){VECTOR<int,4>& nodes=volume.mesh.elements(t);
            int root=union_find.Find(nodes[1]);
            volumes(root)+=TETRAHEDRON<T>::Signed_Volume(particles.X(nodes[1]),particles.X(nodes[2]),particles.X(nodes[3]),particles.X(nodes[4]));}
        // print volumes
        RANGE<VECTOR<T,1> > fragment_volume_box=RANGE<VECTOR<T,1> >(FLT_MAX,-FLT_MAX);
        for(int p=0;p<particles.array_collection->Size();p++)if(union_find.Is_Root(p)){
            if(frame==0) rest_volume=volumes(p);
            fragment_volume_box.Enlarge_To_Include_Point(VECTOR<T,1>(volumes(p)/rest_volume-1));
            volume_box.Enlarge_To_Include_Point(VECTOR<T,1>(volumes(p)/rest_volume-1));}
        LOG::cout<<fragment_volume_box.min_corner.x<<" "<<fragment_volume_box.max_corner.x<<std::endl;}
    LOG::cout<<"rest volume = "<<rest_volume<<std::endl;
    LOG::cout<<"error = "<<volume_box.min_corner.x<<" "<<volume_box.max_corner.x<<std::endl;

    return 0;
}
//#####################################################################
