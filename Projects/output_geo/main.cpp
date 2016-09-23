#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Topology_Based_Geometry/B_SPLINE_PATCH.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;
typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;

HASHTABLE<PAIR<std::string,int>,DEFORMABLE_BODY_COLLECTION<TV>*> deformable_geometry_collection_cache;

DEFORMABLE_BODY_COLLECTION<TV>& Load_Deformable_Geometry_Collection(const std::string& location,int frame)
{
    DEFORMABLE_BODY_COLLECTION<TV>*& deformable_geometry_collection=deformable_geometry_collection_cache.Get_Or_Insert(PAIR<std::string,int>(location,frame));
    if(deformable_geometry_collection) return *deformable_geometry_collection;
    deformable_geometry_collection=new DEFORMABLE_BODY_COLLECTION<TV>(0,0);
    deformable_geometry_collection->Read(STREAM_TYPE(RW()),location,location,frame,-1,true,true);
    return *deformable_geometry_collection;
}

void Emit_Spline_Patch(std::ostream& output,B_SPLINE_PATCH<TV,3>& spline)
{
    output<<"[ [ [ \"type\",\"NURBMesh\" ], [ \"vertex\",["<<std::endl;
    int M=spline.control_points.domain.max_corner(0);
    if(spline.loop_s) M-=3;
    int N=spline.control_points.domain.max_corner(1);
    if(spline.loop_t) N-=3;

    for(int j=0;j<N;j++){
        output<<"[";
        for(int i=0;i<M;i++)
            output<<spline.control_points(i,j)<<",";
        output<<"],";}
    output<<"], \"surface\",\"quads\",";
    output<<" \"uwrap\","<<(spline.loop_s?"true":"false")<<",";
    output<<" \"vwrap\","<<(spline.loop_t?"true":"false")<<",";
    output<<" \"ubasis\",\n[ \"type\",\"NURBS\", \"order\",4,";
    output<<" \"endinterpolation\","<<(spline.loop_s?"false":"true")<<", \"knots\",\n[";
    output<<spline.knots_s(0)<<",";
    for(int i=0;i<spline.knots_s.m;i++)
        output<<spline.knots_s(i)<<",";
    output<<spline.knots_s(spline.knots_s.m-1);
    output<<"] ],";
    output<<"\"vbasis\",[ \"type\",\"NURBS\", \"order\",4,";
    output<<" \"endinterpolation\","<<(spline.loop_t?"false":"true")<<", \"knots\",\n[";
    output<<spline.knots_t(0)<<",";
    for(int i=0;i<spline.knots_t.m;i++)
        output<<spline.knots_t(i)<<",";
    output<<spline.knots_t(spline.knots_t.m-1);
    output<<"] ] ] ] ]"<<std::endl;
}


void Emit_Deformable_Bodies(std::ostream& output,DEFORMABLE_BODY_COLLECTION<TV>& collection)
{
    output<<"[\n";
    output<<"\"pointcount\","<<collection.particles.X.m<<','<<std::endl;
    output<<"\"vertexcount\","<<collection.particles.X.m<<','<<std::endl;
    output<<"\"primitivecount\",1,\n";
    output<<"\"topology\",[\n";
    output<<"\"pointref\",[\n";
    output<<"\"indices\",[\n";
    for(int p=0;p<collection.particles.Size();p++)
        output<<LOG::sprintf("%d,",p);
    output<<"]]],";
    output<<"\"attributes\",[\n";
    output<<"\"pointattributes\",[\n";
    output<<"[                        "<<std::endl;
    output<<"[                        "<<std::endl;
    output<<"\"scope\",\"public\",    "<<std::endl;
    output<<"\"type\",\"numeric\",    "<<std::endl;
    output<<"\"name\",\"P\",          "<<std::endl;
    output<<"],                       "<<std::endl;
    output<<"[                        "<<std::endl;
    output<<"\"size\",4,              "<<std::endl;
    output<<"\"storage\",\"fpreal32\","<<std::endl;
    output<<"\"values\",["<<std::endl;
    output<<"\"size\",4,"<<std::endl;
    output<<"\"storage\",\"fpreal32\","<<std::endl;
    output<<"\"tuples\",["<<std::endl;

    for(int p=0;p<collection.particles.Size();p++)
        output<<LOG::sprintf("[%lg,%lg,%lg,1],\n",collection.particles.X(p)[0],collection.particles.X(p)[1],collection.particles.X(p)[2]);

    output<<"] ] ] ] ] ],"<<std::endl;
    output<<"\"primitives\","<<std::endl;

    for(int i=0;i<collection.structures.m;i++){
        if(B_SPLINE_PATCH<TV,3>* spline=dynamic_cast<B_SPLINE_PATCH<TV,3>*>(collection.structures(i)))
            Emit_Spline_Patch(output,*spline);}
    output<<"]"<<std::endl;
}

int main(int argc,char *argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    int frame_number=0;
    std::string output_filename;
    std::string input_folder;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Extra(&frame_number,"frame number","frame number");
    parse_args.Extra(&output_filename,"geo file","output geo file name");
    parse_args.Extra(&input_folder,"input_folder","input folder name");
    parse_args.Parse();
    if(parse_args.unclaimed_arguments){parse_args.Print_Usage();exit(0);}

    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(output_filename,false);
    DEFORMABLE_BODY_COLLECTION<TV>& collection=Load_Deformable_Geometry_Collection(input_folder,frame_number);
    Emit_Deformable_Bodies(*output,collection);

    delete output;
}
