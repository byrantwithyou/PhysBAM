//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Images/PNG_FILE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Symbolics/PROGRAM.h>
#include <Tools/Symbolics/PROGRAM_CONTEXT.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Fracture/EMBEDDING.h>
#include <fstream>
using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<T,2> TV2;

HASHTABLE<PAIR<std::string,int>,RIGID_BODY_COLLECTION<TV>*> rigid_body_collection_cache;
HASHTABLE<PAIR<std::string,int>,DEFORMABLE_BODY_COLLECTION<TV>*> deformable_geometry_collection_cache;
HASHTABLE<std::string,TRIANGULATED_SURFACE<T>*> saved_surface;
struct TEXTURE
{
    ARRAY<TV2> coords;
    ARRAY<VECTOR<int,3> > map;
};
HASHTABLE<std::string,TEXTURE*> saved_texture;

RIGID_BODY_COLLECTION<TV>& Load_Rigid_Body_Collection(const std::string& location,int frame)
{
    RIGID_BODY_COLLECTION<TV>*& rigid_body_collection=rigid_body_collection_cache.Get_Or_Insert(PAIR<std::string,int>(location,frame));
    if(rigid_body_collection) return *rigid_body_collection;
    rigid_body_collection=new RIGID_BODY_COLLECTION<TV>(0);
    rigid_body_collection->Read(STREAM_TYPE(RW()),location,frame); // TODO: only load implicit surfaces if load_implicit_surfaces
    return *rigid_body_collection;
}

DEFORMABLE_BODY_COLLECTION<TV>& Load_Deformable_Geometry_Collection(const std::string& location,int frame)
{
    DEFORMABLE_BODY_COLLECTION<TV>*& deformable_geometry_collection=deformable_geometry_collection_cache.Get_Or_Insert(PAIR<std::string,int>(location,frame));
    if(deformable_geometry_collection) return *deformable_geometry_collection;
    deformable_geometry_collection=new DEFORMABLE_BODY_COLLECTION<TV>(0,0);
    deformable_geometry_collection->Read(STREAM_TYPE(RW()),location,location,frame,-1,true,true);
    return *deformable_geometry_collection;
}

void Apply_Options(TRIANGULATED_SURFACE<T>* ts, const HASHTABLE<std::string,std::string>& options)
{
    if(options.Contains("preserve_creases")){
        ts->avoid_normal_interpolation_across_sharp_edges=true;}

    if(const std::string* value=options.Get_Pointer("subdivide")){
        int n=atoi(value->c_str());
        if(!n) n=1;
        for(int i=0;i<n;i++) ts->Loop_Subdivide();
        ts->Refresh_Auxiliary_Structures();}

    if(const std::string* value=options.Get_Pointer("scale")){
        T s=atof(value->c_str());
        if(s!=1) ts->Rescale(s);}

    if(const std::string* value=options.Get_Pointer("jitter")){
        T e=atof(value->c_str());
        if(e>0){
            ARRAY<TV> dX(ts->particles.X.m);
            RANDOM_NUMBERS<T> random;
            random.Fill_Uniform(dX,-e,e);
            ts->particles.X+=dX;}}
}

template<class T,int d>
void Emit_Vector(std::ofstream& fout,const VECTOR<T,d>& v, const char* str="")
{
    fout<<"< "<<v(0);
    for(int i=1;i<d;i++) fout<<" , "<<v(i);
    fout<<" > "<<str;
}

void Emit_Smooth_Surface(std::ofstream& fout,TRIANGULATED_SURFACE<T>* ts, const HASHTABLE<std::string,std::string>& options)
{
    ts->Update_Vertex_Normals();

    ARRAY<TV2> coords;
    ARRAY<VECTOR<int,3> > map;
    if(const std::string* texture_map_file=options.Get_Pointer("texture_map")){
        int ignore;
        FILE_UTILITIES::Read_From_File(STREAM_TYPE((RW)0),texture_map_file->c_str(),coords,ignore,map);
        LOG::cout<<"Texture mapping file data:  "<<texture_map_file<<"  "<<coords.m<<"  "<<ignore<<"  "<<map.m<<"  "<<ts->mesh.elements.m<<std::endl;}
    else if(const std::string* str=options.Get_Pointer("saved_texture_map")){
        TEXTURE* tex=saved_texture.Get(*str);
        coords=tex->coords;
        map=tex->map;}

    for(int i=0;i<ts->mesh.elements.m;i++){
        fout<<"smooth_triangle { ";
        VECTOR<TV,3> X(ts->particles.X.Subset(ts->mesh.elements(i)));
        VECTOR<TV,3> N;
        if(ts->face_vertex_normals) N=(*ts->face_vertex_normals)(i);
        else N=VECTOR<TV,3>(ts->vertex_normals->Subset(ts->mesh.elements(i)));

        Emit_Vector(fout,X(0),",");
        Emit_Vector(fout,N(0),",");
        Emit_Vector(fout,X(1),",");
        Emit_Vector(fout,N(1),",");
        Emit_Vector(fout,X(2),",");
        Emit_Vector(fout,N(2),"");

        if(coords.m){
            fout<<"uv_vectors ";
            int uv1,uv2,uv3;map(i).Get(uv1,uv2,uv3);
            Emit_Vector(fout, coords(uv1), " , ");
            Emit_Vector(fout, coords(uv2), " , ");
            Emit_Vector(fout, coords(uv3), "");}

        fout<<"}\n";}
}

void Emit_Rigid_Body(std::ofstream& fout,const HASHTABLE<std::string,std::string>& options,int frame)
{
    RIGID_BODY_COLLECTION<TV>& collection=Load_Rigid_Body_Collection(options.Get("location"),frame);
    RIGID_BODY<TV>& rigid_body=collection.Rigid_Body(atoi(options.Get("index").c_str()));

    TRIANGULATED_SURFACE<T>* ts=rigid_body.simplicial_object->Create_Compact_Copy();
    for(int i=0;i<ts->particles.X.m;i++) ts->particles.X(i)=rigid_body.Frame()*ts->particles.X(i);

    Apply_Options(ts,options);
    if(const std::string* value=options.Get_Pointer("save"))
        saved_surface.Get_Or_Insert(*value)=ts;
    else Emit_Smooth_Surface(fout,ts,options);
}

void Emit_Rigid_Body_Frame(std::ofstream& fout,const HASHTABLE<std::string,std::string>& options,int frame)
{
    RIGID_BODY_COLLECTION<TV>& collection=Load_Rigid_Body_Collection(options.Get("location"),frame);
    RIGID_BODY<TV>& rigid_body=collection.Rigid_Body(atoi(options.Get("index").c_str()));
    MATRIX<T,3> rot=rigid_body.Frame().r.Rotation_Matrix();
    TV X=rigid_body.Frame().t;
    fout<<"matrix < ";
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            fout<<rot(j,i)<<" , ";
    fout<<X.x<<" , "<<X.y<<" , "<<X.z<<" >\n";
}

void Emit_Deformable_Body(std::ofstream& fout,const HASHTABLE<std::string,std::string>& options,int frame)
{
    if(const std::string* value=options.Get_Pointer("frame")) frame=atoi(value->c_str());
    DEFORMABLE_BODY_COLLECTION<TV>& collection=Load_Deformable_Geometry_Collection(options.Get("location"),frame);
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
    if(const std::string* value=options.Get_Pointer("save"))
        saved_surface.Get_Or_Insert(*value)=ts;
    else Emit_Smooth_Surface(fout,ts,options);
}

struct PATCH
{
    PROGRAM<T> map_prog,color_prog;
    PROGRAM_CONTEXT<T> map_context,color_context;
    ARRAY<int> patch_elements;
    ARRAY<int> patch_vertices;
    ARRAY<TV2> coords;
    RANGE<TV2> uv_range;
    int x_start;
    TV2 offset;

    PATCH(){}

    void Set_Map_Program(const char* map_str)
    {
        map_prog.var_in.Append("x");
        map_prog.var_in.Append("y");
        map_prog.var_in.Append("z");
        map_prog.var_in.Append("patch");
        map_prog.var_out.Append("u");
        map_prog.var_out.Append("v");
        map_prog.Parse(map_str,false);
        map_prog.Optimize();
        map_prog.Finalize();
        map_context.Initialize(map_prog);
    }

    void Set_Color_Program(const char* color_str)
    {
        color_prog.var_in.Append("u");
        color_prog.var_in.Append("v");
        color_prog.var_in.Append("patch");
        color_prog.var_out.Append("r");
        color_prog.var_out.Append("g");
        color_prog.var_out.Append("b");
        color_prog.var_out.Append("a");
        color_prog.Parse(color_str,false);
        color_prog.Optimize();
        color_prog.Finalize();
        color_context.Initialize(color_prog);
    }
};

void Create_Texture_Map(std::ofstream& fout,const HASHTABLE<std::string,std::string>& options,int frame)
{
    TRIANGULATED_SURFACE<T>* ts=saved_surface.Get(options.Get("surface"));

    int samples=250000;
    if(const std::string* sample_str=options.Get_Pointer("samples")) samples=atoi(sample_str->c_str());
    bool use_alpha=options.Contains("alpha"); 

    int num_patches=atoi(options.Get("num_patches").c_str());
    ARRAY<PATCH> patches(num_patches);

    std::stringstream mss(options.Get("map_funcs").c_str());
    std::stringstream css(options.Get("color_funcs").c_str());
    for(int i=0;i<num_patches;i++){
        PATCH& p=patches(i);
        std::string m,c;
        mss>>m;
        css>>c;
        p.Set_Map_Program(options.Get(m).c_str());
        p.Set_Color_Program(options.Get(c).c_str());}

    std::string patch_str=options.Get("patch");
    PROGRAM<T> patch_program;
    patch_program.var_in.Append("x@");
    patch_program.var_in.Append("y@");
    patch_program.var_in.Append("z@");
    patch_program.var_in.Append("x0");
    patch_program.var_in.Append("y0");
    patch_program.var_in.Append("z0");
    patch_program.var_in.Append("x1");
    patch_program.var_in.Append("y1");
    patch_program.var_in.Append("z1");
    patch_program.var_in.Append("nx");
    patch_program.var_in.Append("ny");
    patch_program.var_in.Append("nz");
    patch_program.var_in.Append("index");
    patch_program.var_out.Append("patch");
    patch_program.Parse(patch_str.c_str(),false);
    patch_program.Optimize();
    patch_program.Finalize();

    TEXTURE* tex=new TEXTURE;
    PROGRAM_CONTEXT<T> context(patch_program);
    ARRAY<int> element_patch,element_map,element_color;
    for(int i=0;i<ts->mesh.elements.m;i++){
        TRIANGLE_3D<T> tri(ts->Get_Element(i));
        TV n=tri.Normal();
        for(int j=0;j<TV::m;j++) for(int k=0;k<TV::m;k++) context.data_in(TV::m*j+k)=tri.X(j)(k);
        for(int j=0;j<TV::m;j++) context.data_in(TV::m*TV::m+j)=n(j);
        context.data_in(TV::m*(TV::m+1))=i;
        patch_program.Execute(context);
        int patch=context.data_out(0);
        tex->map.Append(ts->mesh.elements(i)+ts->particles.X.m*patch);
        patches(patch).patch_elements.Append(i);}

    T uv_area=0;
    for(int i=0;i<num_patches;i++){
        PATCH& p=patches(i);
        for(int j=0;j<p.patch_elements.m;j++)
            p.patch_vertices.Append_Elements(ts->mesh.elements(p.patch_elements(j)));
        Prune_Duplicates(p.patch_vertices);
        p.coords.Resize(ts->particles.X.m);
        for(int j=0;j<p.patch_vertices.m;j++){
            int q=p.patch_vertices(j);
            p.map_context.data_in=ts->particles.X(q).Append(i);
            p.map_prog.Execute(p.map_context);
            p.coords(q)=TV2(p.map_context.data_out);}
        for(int j=0;j<p.patch_elements.m;j++){
            T area=TRIANGLE_2D<T>::Size(p.coords.Subset(ts->mesh.elements(p.patch_elements(j))));
            uv_area+=area;}
        p.uv_range=RANGE<TV2>::Bounding_Box(p.coords.Subset(p.patch_vertices));}

    T scale=sqrt(samples/uv_area);

    int margin=2;
    int start=0,height=0;
    for(int i=0;i<num_patches;i++){
        PATCH& p=patches(i);
        p.x_start=start;
        p.offset=TV2(p.x_start,0)+margin-scale*p.uv_range.min_corner;
        TV2 end=scale*p.uv_range.max_corner+p.offset+margin;
        start=end.x;
        height=std::max(height,(int)end.y);}

    ARRAY<VECTOR<T,4>,VECTOR<int,2> > image(VECTOR<int,2>(start,height)+1);

    for(int i=0;i<num_patches;i++){
        PATCH& p=patches(i);
        RANGE<VECTOR<int,2> > rast_box(scale*p.uv_range+p.offset);
        for(RANGE_ITERATOR<2> it(rast_box.To_Half_Opened().Thickened(margin));it.Valid();it.Next()){
            TV2 uv=(TV2(it.index)-p.offset)/scale;
            p.color_context.data_in=uv.Append(i);
            p.color_prog.Execute(p.color_context);
            VECTOR<T,4> c(p.color_context.data_out);
            if(!use_alpha) c(3)=1;
            image(it.index)=c;}
        for(int j=0;j<p.coords.m;j++)
            tex->coords.Append((p.coords(j)*scale+p.offset)/TV2(start,height));}

    char buff[1000];
    sprintf(buff,options.Get("png_file").c_str(),frame);
    PNG_FILE<T>::Write(buff,image);

    saved_texture.Get_Or_Insert(options.Get("uv_save"))=tex;
    if(const std::string* str=options.Get_Pointer("uv_file")){
        sprintf(buff,str->c_str(),frame);
        FILE_UTILITIES::Write_To_File(STREAM_TYPE((RW)0),buff,tex->coords,0,tex->map);}
}

bool Parse_Pair(const char*& str,std::string& key,std::string& value)
{
    int a=strcspn(str, " \t=");
    if(a==0) return false;
    key=std::string(str,a);
    str+=a;
    if(*str!='=')
    {
        value="";
        str+=strspn(str, " \t");
        return true;
    }
    str++;
    char ec[]=" \t\n";
    if(*str=='\'' || *str=='"')
    {
        ec[0]=*str;
        ec[1]=0;
        str++;
    }

    int len=strcspn(str, ec);
    value=std::string(str, len);
    str+=len;
    if(*str) str++;
    str+=strspn(str, " \t");
    return true;
}

void Emit_Camera(std::ofstream& fout,const HASHTABLE<std::string,std::string>& options)
{
    TV loc, at, up;
    T angle;
    std::string file=options.Get("location"), line;
    std::ifstream fin(file.c_str());
    while(getline(fin, line))
    {
        VECTOR<double,TV::m> dloc,dup,dat;
        double dangle;
        if(sscanf(line.c_str(), " Location = ( %lg %lg %lg ) \n", &dloc.x, &dloc.y, &dloc.z)==3){loc=TV(dloc);}
        else if(sscanf(line.c_str(), " Pseudo_Up = ( %lg %lg %lg ) \n", &dup.x, &dup.y, &dup.z)==3){up=TV(dup);}
        else if(sscanf(line.c_str(), " Look_At = ( %lg %lg %lg ) \n", &dat.x, &dat.y, &dat.z)==3){at=TV(dat);}
        else if(sscanf(line.c_str(), " Field_Of_View = %lg \n", &dangle)==1){angle=dangle;}
    }
    fout<<"camera { location ";
    Emit_Vector(fout, loc, " look_at ");
    Emit_Vector(fout, at, " up ");
    Emit_Vector(fout, up, " right ");
    Emit_Vector(fout, TV(-4./3,0,0), " angle ");
    fout<<angle<<" }"<<std::endl;
}

int main(int argc, char *argv[]) 
{  
    PROCESS_UTILITIES::Set_Backtrace(true);

    std::string scene_filename,output_filename;
    int frame_number=0;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Extra(&scene_filename,"scene file","scene file");
    parse_args.Extra(&output_filename,"output scene file","output scene file");
    parse_args.Extra(&frame_number,"frame number","frame number");
    parse_args.Parse();
    if(parse_args.unclaimed_arguments){parse_args.Print_Usage();exit(0);}

    std::ifstream fin(scene_filename.c_str());
    std::ofstream fout(output_filename.c_str());
    std::string line,more;
    while(getline(fin,line))
    {
        while(line[line.size()-1]=='\\' && getline(fin,more)){
            line.resize(line.size()-1);
            line+=more;}

        const char* str=line.c_str();
        int sp=strspn(str, " \t");

        if((int)line.size()<sp+5 || std::string(str+sp,5)!="#emit"){
            fout<<line<<std::endl;
            continue;}

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
        else if(type=="camera")
            Emit_Camera(fout,options);
        else if(type=="rigid_body_frame")
            Emit_Rigid_Body_Frame(fout,options,frame_number);
        else if(type=="triangulated_surface")
            Emit_Rigid_Body_Frame(fout,options,frame_number);
        else if(type=="texture_map")
            Create_Texture_Map(fout,options,frame_number);
        else PHYSBAM_FATAL_ERROR("unexpected replacement type: '"+type+"'.");
    }


    return 0;
}
//#####################################################################
