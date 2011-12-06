//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDING.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
#include <PhysBAM_Dynamics/Read_Write/Particles/READ_WRITE_PARTICLES.h>
#include <fstream>
using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;

HASHTABLE<PAIR<std::string,int>,RIGID_BODY_COLLECTION<TV>*> rigid_body_collection_cache;
HASHTABLE<PAIR<std::string,int>,DEFORMABLE_GEOMETRY_COLLECTION<TV>*> deformable_geometry_collection_cache;

RIGID_BODY_COLLECTION<TV>& Load_Rigid_Body_Collection(const std::string& location,int frame)
{
    RIGID_BODY_COLLECTION<TV>*& rigid_body_collection=rigid_body_collection_cache.Get_Or_Insert(PAIR<std::string,int>(location,frame));
    if(rigid_body_collection) return *rigid_body_collection;
    rigid_body_collection=new RIGID_BODY_COLLECTION<TV>(0,0);
    rigid_body_collection->Read(STREAM_TYPE(RW()),location,frame); // TODO: only load implicit surfaces if load_implicit_surfaces
    return *rigid_body_collection;
}

DEFORMABLE_GEOMETRY_COLLECTION<TV>& Load_Deformable_Geometry_Collection(const std::string& location,int frame)
{
    DEFORMABLE_GEOMETRY_COLLECTION<TV>*& deformable_geometry_collection=deformable_geometry_collection_cache.Get_Or_Insert(PAIR<std::string,int>(location,frame));
    if(deformable_geometry_collection) return *deformable_geometry_collection;
    deformable_geometry_collection=new DEFORMABLE_GEOMETRY_COLLECTION<TV>(*new GEOMETRY_PARTICLES<TV>);
    deformable_geometry_collection->Read(STREAM_TYPE(RW()),location,location,frame,-1,true);
    return *deformable_geometry_collection;
}

void Apply_Options(TRIANGULATED_SURFACE<T>* ts, const HASHTABLE<std::string,std::string>& options)
{
    if(const std::string* value=options.Get_Pointer("subdivide")){
        int n=atoi(value->c_str());
        if(!n) n=1;
        for(int i=1;i<=n;i++) ts->Loop_Subdivide();
        ts->Refresh_Auxiliary_Structures();}

    if(const std::string* value=options.Get_Pointer("scale")){
        T s=atof(value->c_str());
        if(s!=1) ts->Rescale(s);}

    if(options.Contains("preserve_creases")){
        ts->avoid_normal_interpolation_across_sharp_edges=true;}
}

void Emit_Vector(std::ofstream& fout,const TV& v, const char* str="")
{
    fout<<"< "<<v.x<<" , "<<v.y<<" , "<<v.z<<" > "<<str;
}

void Emit_Smooth_Surface(std::ofstream& fout,TRIANGULATED_SURFACE<T>* ts)
{
    ts->Update_Vertex_Normals();
    for(int i=1;i<=ts->mesh.elements.m;i++){
        fout<<"smooth_triangle { ";
        VECTOR<TV,3> X(ts->particles.X.Subset(ts->mesh.elements(i)));
        VECTOR<TV,3> N;
        if(ts->face_vertex_normals) N=(*ts->face_vertex_normals)(i);
        else N=VECTOR<TV,3>(ts->vertex_normals->Subset(ts->mesh.elements(i)));

        Emit_Vector(fout,X(1),",");
        Emit_Vector(fout,N(1),",");
        Emit_Vector(fout,X(2),",");
        Emit_Vector(fout,N(2),",");
        Emit_Vector(fout,X(3),",");
        Emit_Vector(fout,N(3),"}\n");}
}

void Emit_Rigid_Body(std::ofstream& fout,const HASHTABLE<std::string,std::string>& options,int frame)
{
    RIGID_BODY_COLLECTION<TV>& collection=Load_Rigid_Body_Collection(options.Get("location"),frame);
    RIGID_BODY<TV>& rigid_body=collection.Rigid_Body(atoi(options.Get("index").c_str()));

    TRIANGULATED_SURFACE<T>* ts=rigid_body.simplicial_object->Create_Compact_Copy();
    for(int i=1;i<=ts->particles.X.m;i++) ts->particles.X(i)=rigid_body.Frame()*ts->particles.X(i);

    Apply_Options(ts,options);
    Emit_Smooth_Surface(fout,ts);
}

void Emit_Deformable_Body(std::ofstream& fout,const HASHTABLE<std::string,std::string>& options,int frame)
{
    DEFORMABLE_GEOMETRY_COLLECTION<TV>& collection=Load_Deformable_Geometry_Collection(options.Get("location"),frame);
    STRUCTURE<TV>* structure=collection.structures(atoi(options.Get("index").c_str()));
    TRIANGULATED_SURFACE<T>* ts=0;

    if(TRIANGULATED_SURFACE<T>* triangulated_surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(structure)){
        ts=triangulated_surface->Create_Compact_Copy();}
    else if(TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(structure)){
        tetrahedralized_volume->Update_Number_Nodes();
        tetrahedralized_volume->Initialize_Triangulated_Surface();
        ts=tetrahedralized_volume->triangulated_surface->Create_Compact_Copy();}
    else if(HEXAHEDRALIZED_VOLUME<T>* hexahedralized_volume=dynamic_cast<HEXAHEDRALIZED_VOLUME<T>*>(structure)){
        hexahedralized_volume->Initialize_Triangulated_Surface();
        ts=hexahedralized_volume->triangulated_surface->Create_Compact_Copy();}
    else if(EMBEDDING<TV>* embedding=dynamic_cast<EMBEDDING<TV>*>(structure)){
        ts=embedding->material_surface.Create_Compact_Copy();}
    else PHYSBAM_FATAL_ERROR(std::string("Unhandled object type: ")+typeid(*structure).name());

    Apply_Options(ts,options);
    Emit_Smooth_Surface(fout,ts);
}

bool Parse_Pair(const char*& str,std::string& key,std::string& value)
{
    int a=strcspn(str, " \t=");
    if(a==0) return false;
    key=std::string(str,a);
    str+=a+1;
    if(str[-1]!='=')
    {
        value="";
        str+=strspn(str, " \t");
        return true;
    }
    char ec=' ';
    if(*str=='\'' || *str=='"')
    {
        ec=*str;
        str++;
    }

    const char* end=strchr(str, ec);
    PHYSBAM_ASSERT(end);
    value=std::string(str, end-str);
    str=end+1;
    str+=strspn(str, " \t");
    return true;
}

int main(int argc, char *argv[]) 
{  
    PROCESS_UTILITIES::Set_Backtrace(true);
    Initialize_Particles();
    Initialize_Read_Write_General_Structures();

    PARSE_ARGS parse_args;
    parse_args.Set_Extra_Arguments(-1, "<scene file> <output scene file> <frame number>");
    parse_args.Parse(argc,argv);
    if(parse_args.Num_Extra_Args() != 3){parse_args.Print_Usage();exit(0);}
    std::string scene_filename=parse_args.Extra_Arg(1);
    std::string output_filename=parse_args.Extra_Arg(2);
    int frame_number=atoi(parse_args.Extra_Arg(3).c_str());

    std::ifstream fin(scene_filename.c_str());
    std::ofstream fout(output_filename.c_str());
    std::string line;
    while(getline(fin,line))
    {
        const char* str=line.c_str();
        int sp=strspn(str, " \t");
        
        if((int)line.size()<sp+5 || std::string(str+sp,5)!="#emit")
        {
            fout<<line<<std::endl;
            continue;
        }

        str+=sp+5;
        str+=strspn(str, " \t");
        sp=strcspn(str, " \t");
        std::string type(str, sp);
        str+=sp;
        str+=strspn(str, " \t");

        HASHTABLE<std::string,std::string> options;
        std::string key,value;
        while(Parse_Pair(str,key,value))
            options.Set(key,value);

        if(type=="rigid_body")
            Emit_Rigid_Body(fout,options,frame_number);
        else if(type=="deformable_body")
            Emit_Deformable_Body(fout,options,frame_number);
        else PHYSBAM_FATAL_ERROR("unexpected replacement type: "+type);
    }


    return 0;
}
//#####################################################################
