#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
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

HASHTABLE<std::string,DEFORMABLE_BODY_COLLECTION<TV>*> deformable_body_collection_cache;

DEFORMABLE_BODY_COLLECTION<TV>& Load_Deformable_Body_Collection(const VIEWER_DIR& viewer_dir)
{
    DEFORMABLE_BODY_COLLECTION<TV>*& deformable_body_collection=deformable_body_collection_cache.Get_Or_Insert(viewer_dir.current_directory);
    if(deformable_body_collection) return *deformable_body_collection;
    deformable_body_collection=new DEFORMABLE_BODY_COLLECTION<TV>(0,0);
    deformable_body_collection->Read(viewer_dir,true);
    return *deformable_body_collection;
}

void Emit_OpenSubdiv_Surface(std::ostream& output,OPENSUBDIV_SURFACE<TV,3>& surf)
{
    output<<"g\n";
    for(int i=0;i<surf.m;i++)
        LOG::fprintf(output,"f %d %d %d %d\n",surf.control_points(surf.mesh(i)(0))+1,surf.control_points(surf.mesh(i)(1))+1,surf.control_points(surf.mesh(i)(2))+1,surf.control_points(surf.mesh(i)(3))+1);
}

void Emit_Deformable_Bodies(std::ostream& output,DEFORMABLE_BODY_COLLECTION<TV>& collection)
{
    for(int p=0;p<collection.particles.Size();p++)
        LOG::fprintf(output,"v %lg %lg %lg\n",collection.particles.X(p)[0],collection.particles.X(p)[1],collection.particles.X(p)[2]);
    
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
    std::string input_directory;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Extra(&frame_number,"frame number","frame number");
    parse_args.Extra(&output_filename,"obj file","output obj file name");
    parse_args.Extra(&input_directory,"input_directory","input directory name");
    parse_args.Parse();
    if(parse_args.unclaimed_arguments){parse_args.Print_Usage();exit(0);}

    std::ostream* output=Safe_Open_Output_Raw(output_filename,false);
    VIEWER_DIR viewer_dir(input_directory);
    viewer_dir.Set(frame_number);
    DEFORMABLE_BODY_COLLECTION<TV>& collection=Load_Deformable_Body_Collection(viewer_dir);
    Emit_Deformable_Bodies(*output,collection);

    delete output;
}
