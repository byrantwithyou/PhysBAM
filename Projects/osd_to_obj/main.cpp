#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Topology_Based_Geometry/B_SPLINE_PATCH.h>
#include <Geometry/Topology_Based_Geometry/OPENSUBDIV_SURFACE.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;
typedef double T;
typedef VECTOR<T,3> TV;

HASHTABLE<PAIR<std::string,int>,DEFORMABLE_BODY_COLLECTION<TV>*> deformable_geometry_collection_cache;

DEFORMABLE_BODY_COLLECTION<TV>& Load_Deformable_Geometry_Collection(const std::string& location,int frame,bool use_doubles)
{
    DEFORMABLE_BODY_COLLECTION<TV>*& deformable_geometry_collection=deformable_geometry_collection_cache.Get_Or_Insert(PAIR<std::string,int>(location,frame));
    if(deformable_geometry_collection) return *deformable_geometry_collection;
    deformable_geometry_collection=new DEFORMABLE_BODY_COLLECTION<TV>(0,0);
    if(use_doubles) deformable_geometry_collection->Read(STREAM_TYPE(double()),location,location,frame,-1,true,true);
    else deformable_geometry_collection->Read(STREAM_TYPE(float()),location,location,frame,-1,true,true);
    return *deformable_geometry_collection;
}

void Emit_OpenSubdiv_Surface(std::ostream& output,OPENSUBDIV_SURFACE<TV,3>& surf)
{
    output<<"g\n";
    for(int i=0;i<surf.m;i++)
        output<<LOG::sprintf("f %d %d %d %d\n",surf.control_points(surf.mesh(i)(0))+1,surf.control_points(surf.mesh(i)(1))+1,surf.control_points(surf.mesh(i)(2))+1,surf.control_points(surf.mesh(i)(3))+1);
}

void Emit_Deformable_Bodies(std::ostream& output,DEFORMABLE_BODY_COLLECTION<TV>& collection)
{
    for(int p=0;p<collection.particles.Size();p++)
        output<<LOG::sprintf("v %lg %lg %lg\n",collection.particles.X(p)[0],collection.particles.X(p)[1],collection.particles.X(p)[2]);
    
    for(int i=0;i<collection.structures.m;i++){
        if(OPENSUBDIV_SURFACE<TV,3>* surf=dynamic_cast<OPENSUBDIV_SURFACE<TV,3>*>(collection.structures(i)))
            Emit_OpenSubdiv_Surface(output,*surf);
    }
}

int main(int argc,char *argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    int frame_number=0;
    std::string output_filename;
    std::string input_folder;
    bool use_doubles=false;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Extra(&frame_number,"frame number","frame number");
    parse_args.Extra(&output_filename,"obj file","output obj file name");
    parse_args.Extra(&input_folder,"input_folder","input folder name");
    parse_args.Add("-double",&use_doubles,"read doubles instead of floats");
    parse_args.Parse();
    if(parse_args.unclaimed_arguments){parse_args.Print_Usage();exit(0);}

    std::ostream* output=Safe_Open_Output(output_filename,false);
    DEFORMABLE_BODY_COLLECTION<TV>& collection=Load_Deformable_Geometry_Collection(input_folder,frame_number,use_doubles);
    Emit_Deformable_Bodies(*output,collection);

    delete output;
}
