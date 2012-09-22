//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Particles/PARTICLES_FORWARD.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <fstream>
#include <iostream>
#include <ostream>

using namespace PhysBAM;

template<class TV>
void Dump_Rigid_Baody_Simplicial_Object_Particles(const STREAM_TYPE& stream_type,std::ostream* output_stream,const std::string& basedir,const int frame,
    const int rigid_body_id)
{
    ARRAY<int> needs_init,needs_destroy;
    RIGID_BODY_COLLECTION<TV> rigid_body_collection(0,0);
    rigid_body_collection.Read(stream_type,basedir,0,&needs_init);
    rigid_body_collection.Read(stream_type,basedir,frame,&needs_init);
    
    RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection=rigid_body_collection.rigid_geometry_collection;
    RIGID_GEOMETRY<TV>& rigid_geometry=rigid_geometry_collection.Rigid_Geometry(rigid_body_id);

    *output_stream<<"Number of Particles="<<rigid_geometry.simplicial_object->particles.Size()<<std::endl;
    for(int t=0;t<rigid_geometry.simplicial_object->particles.Size();t++){
        TV X=rigid_geometry.simplicial_object->particles.X(t);
        *output_stream<<"X("<<t<<")="<<rigid_geometry.Frame()*X<<std::endl;}
}
int main(int argc,char *argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    bool type_double=false;
    int dimension=1,frame=0,rigid_body_id=1;
    std::string output_file;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-d",&dimension,"dim","dimension");
    parse_args.Add("-frame",&frame,"frame","frame number");
    parse_args.Add("-rbid",&rigid_body_id,"id","rigid body id");
    parse_args.Add("-o",&output_file,"file","output file");
    parse_args.Set_Extra_Arguments(1,"<basedir>","base directory");
    parse_args.Parse();

    std::ostream* output_stream=&(std::cout);
    std::ofstream output_file_stream;
    if(output_file.size()){
        output_file_stream.open(output_file.c_str());
        output_stream=&output_file_stream;}

    std::string basedir=parse_args.Extra_Arg(0);

#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
        if(type_double) PHYSBAM_FATAL_ERROR("No double support");
#endif

    STREAM_TYPE stream_type(type_double?STREAM_TYPE(0.0):STREAM_TYPE(0.0f));

    if(!type_double){
        switch(dimension){
            case 1:
                Dump_Rigid_Baody_Simplicial_Object_Particles<VECTOR<float,1> >(stream_type,output_stream,basedir,frame,rigid_body_id);
                break;
            case 2:
                Dump_Rigid_Baody_Simplicial_Object_Particles<VECTOR<float,2> >(stream_type,output_stream,basedir,frame,rigid_body_id);
                break;
            case 3:
                Dump_Rigid_Baody_Simplicial_Object_Particles<VECTOR<float,3> >(stream_type,output_stream,basedir,frame,rigid_body_id);
                break;}}
    else{
        switch(dimension){
            case 1:
                Dump_Rigid_Baody_Simplicial_Object_Particles<VECTOR<double,1> >(stream_type,output_stream,basedir,frame,rigid_body_id);
                break;
            case 2:
                Dump_Rigid_Baody_Simplicial_Object_Particles<VECTOR<double,2> >(stream_type,output_stream,basedir,frame,rigid_body_id);
                break;
            case 3:
                Dump_Rigid_Baody_Simplicial_Object_Particles<VECTOR<double,3> >(stream_type,output_stream,basedir,frame,rigid_body_id);
                break;}}
    output_file_stream.close();
}
