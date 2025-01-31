//#####################################################################
// Copyright 2002-2012, Robert Bridson, Eran Guendelman, Eilene Hao, Geoffrey Irving, Cynthia Lau, Neil Molino, Avi Robinson-Mosher, Andrew Selle, Eftychios Sifakis, Alexey Stomakhin, Jerry Talton, Joseph Teran, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Log/SCOPE.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Particles/PARTICLES_FORWARD.h>
#include <Geometry/Grids_Uniform_Computations/DUALCONTOUR_2D.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <OpenGL/OpenGL/BASIC_VISUALIZATION.h>
#include <OpenGL/OpenGL/OPENGL_AXES.h>
#include <OpenGL/OpenGL/OPENGL_BOX_3D.h>
#include <OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <OpenGL/OpenGL/OPENGL_HEXAHEDRALIZED_VOLUME.h>
#include <OpenGL/OpenGL/OPENGL_LEVELSET_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_LEVELSET_MULTIVIEW.h>
#include <OpenGL/OpenGL/OPENGL_LIGHT.h>
#include <OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_3D.h>
#include <OpenGL/OpenGL/OPENGL_TETRAHEDRALIZED_VOLUME.h>
#include <OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <OpenGL/OpenGL/OPENGL_WINDOW.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <cstring>
#include <fstream>
using namespace PhysBAM;
template<class T> void Add_File(OPENGL_WORLD<T>& opengl_world,const std::string& filename,int number);
template<class T> void Add_Tri2D_File(const std::string& filename,OPENGL_WORLD<T>& world,int number);
template<class T> void Add_Tri_File(const std::string& filename,OPENGL_WORLD<T>& world,int number);
template<class T> void Add_Obj_File(const std::string& filename,OPENGL_WORLD<T>& world,int number);
template<class T> void Add_Phi_File(const std::string& filename,OPENGL_WORLD<T>& world,int number);
template<class T> void Add_Phi2D_File(const std::string& filename,OPENGL_WORLD<T>& world,int number);
template<class T> void Add_Curve_File(const std::string& filename,OPENGL_WORLD<T>& world,int number);
template<class T> void Add_Curve2D_File(const std::string& filename,OPENGL_WORLD<T>& world,int number);
template<class T> void Add_Tet_File(const std::string& filename,OPENGL_WORLD<T>& world,int number);
template<class T> void Add_Hex_File(const std::string& filename,OPENGL_WORLD<T>& world,int number);
template<class T> void Add_Box_File(const std::string& filename,OPENGL_WORLD<T>& world,int number);
//#################################################################
static bool triangulated_surface_highlight_boundary=false;
static bool print_statistics=true;
//#################################################################
// Function main
//#################################################################

template<class T>
class VIEW_GEOMETRY_VISUALIZATION:public BASIC_VISUALIZATION<T>
{
  public:
    ARRAY<std::string> files;

    VIEW_GEOMETRY_VISUALIZATION()
    {
    }

    void Setup_And_Run(PARSE_ARGS &parse_args)
    {
        parse_args.Extra_Optional(&files,0,"files","Files to be viewed");
        this->Initialize(parse_args);
        this->opengl_window_title="PhysBAM geometry viewer";
        this->Run();
        
        // world.Add_Light(new OPENGL_LIGHT(VECTOR<double,3>(1,1,1),.8));
        // world.Add_Light(new OPENGL_LIGHT(VECTOR<double,3>(.3,-2,.1),.4));
        // world.Add_Light(new OPENGL_LIGHT(VECTOR<double,3>(-2,.1,.5),.2));
        // world.Set_Ambient_Light(.2);
        // world.Initialize("PhysBAM geometry viewer");

        // world.Bind_Key("^c",{[&world](){world.Save_View("camera_script",true);},"Save view"});
        // world.Bind_Key('c',{[&world](){world.Load_View("camera_script",true);},"Load view"});
    }

    virtual void Initialize_Components_And_Key_Bindings()
    {
        BASIC_VISUALIZATION<T>::Initialize_Components_And_Key_Bindings();
        std::map<std::string,std::function<void(const std::string&,OPENGL_WORLD<T>&,int)> > func_map;
        func_map["tri"]=Add_Tri_File<T>;
        func_map["obj"]=Add_Obj_File<T>;
        func_map["tri2d"]=Add_Tri2D_File<T>;
        func_map["phi"]=Add_Phi_File<T>;
        func_map["phi2d"]=Add_Phi2D_File<T>;
        func_map["curve"]=Add_Curve_File<T>;
        func_map["curve2d"]=Add_Curve2D_File<T>;
        func_map["tet"]=Add_Tet_File<T>;
        func_map["hex"]=Add_Hex_File<T>;
        func_map["box"]=Add_Box_File<T>;

        for(int i=0;i<files.m;i++)
        {
            auto it=func_map.find(Get_File_Extension(files(i)));
            if(it!=func_map.end()) it->second(files(i),this->opengl_world,i);
            else LOG::cerr<<"Unrecognized file "<<files(i)<<std::endl;
        }
        LOG::cout<<std::flush;
    }
};

int main(int argc,char *argv[])
{
    bool type_double=false;    
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Compute in floats");
    parse_args.Add("-double",&type_double,"Compute in doubles");
    parse_args.Add("-show_boundary",&triangulated_surface_highlight_boundary,"Highlight boundary");
    parse_args.Add_Not("-nostat",&print_statistics,"Disable statistics");
    parse_args.Parse(true);
    
    if(type_double){
        VIEW_GEOMETRY_VISUALIZATION<double> v;
        v.Setup_And_Run(parse_args);}
    else{
        VIEW_GEOMETRY_VISUALIZATION<float> v;
        v.Setup_And_Run(parse_args);}
    
    return 0;
}
//#################################################################
// Function Add_Box_File
//#################################################################
template<class T> void Add_Box_File(const std::string& filename,OPENGL_WORLD<T>& world,int number)
{
//    typedef VECTOR<T,3> TV;
    try{
        // BOX<TV>* box=new BOX<TV>;
        // Read_From_File<T>(filename,*box);
        // OPENGL_BOX_3D<T>* ob=new OPENGL_BOX_3D<T>(*box);
        // world.Add_Object(ob,true,true);
    }
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Hex_File
//#################################################################
template<class T> void Add_Hex_File(const std::string& filename,OPENGL_WORLD<T>& world,int number)
{
    try{
        HEXAHEDRALIZED_VOLUME<T>* hex_vol;
        Create_From_File(filename,hex_vol);
        OPENGL_HEXAHEDRALIZED_VOLUME<T>* ohv=new OPENGL_HEXAHEDRALIZED_VOLUME<T>(&hex_vol->mesh,&hex_vol->particles,OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.7),float(1),float(.8))),OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.7),float(8),float(.1))));
        ohv->selectable=true;
        world.Add_Object(ohv,true,true);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Tri2D_File
//#################################################################
template<class T> void Add_Tri2D_File(const std::string& filename,OPENGL_WORLD<T>& world,int number)
{
    try{
        TRIANGULATED_AREA<T>* area;
        Create_From_File(filename,area);
        TRIANGULATED_SURFACE<T>* surface=new TRIANGULATED_SURFACE<T>(area->mesh,*new DEFORMABLE_PARTICLES<VECTOR<T,3> >);
        surface->particles.Add_Elements(area->particles.Size());
        for(int p=0;p<area->particles.Size();p++)surface->particles.X(p)=VECTOR<T,3>(area->particles.X(p));
        area->Update_Bounding_Box();
        LOG::cout<<"bounding box: "<<*area->bounding_box<<std::endl;
        OPENGL_TRIANGULATED_SURFACE<T>* opengl_triangulated_surface=new OPENGL_TRIANGULATED_SURFACE<T>(*surface,false,OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Cyan()));
        world.Bind_Key('0'+number,{[opengl_triangulated_surface](){static bool is_two_sided=false;is_two_sided=!is_two_sided;opengl_triangulated_surface->Set_Two_Sided(is_two_sided);},"Toggle Two Sided"});
        opengl_triangulated_surface->selectable=true;
        world.Add_Object(opengl_triangulated_surface,true,true);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Tri_File
//#################################################################
template<class T> void Add_Tri_File(const std::string& filename,OPENGL_WORLD<T>& world,int number)
{
    try{
        TRIANGULATED_SURFACE<T>* surface;
        Create_From_File(filename,surface);
        {LOG::SCOPE scope("mesh statistics","mesh statistics");
            LOG::cout<<"filename = "<<filename<<std::endl;
            if(print_statistics) surface->Print_Statistics(LOG::cout);}
        OPENGL_TRIANGULATED_SURFACE<T>* opengl_triangulated_surface=new OPENGL_TRIANGULATED_SURFACE<T>(*surface,false,
            OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Red()),OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Blue()));
        if(triangulated_surface_highlight_boundary) opengl_triangulated_surface->highlight_boundary=true;
        world.Bind_Key('0'+number,{[opengl_triangulated_surface](){static bool is_two_sided=false;is_two_sided=!is_two_sided;opengl_triangulated_surface->Set_Two_Sided(is_two_sided);},"Toggle Two Sided"});
        opengl_triangulated_surface->selectable=true;
        world.Add_Object(opengl_triangulated_surface,true,true);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Tri_File
//#################################################################
template<class T> void Add_Obj_File(const std::string& filename,OPENGL_WORLD<T>& world,int number)
{
    try{
        TRIANGULATED_SURFACE<T>* surface=new TRIANGULATED_SURFACE<T>;
        surface->Read_Obj(filename);
        {LOG::SCOPE scope("mesh statistics","mesh statistics");
            LOG::cout<<"filename = "<<filename<<std::endl;
            if(print_statistics) surface->Print_Statistics(LOG::cout);}
        OPENGL_TRIANGULATED_SURFACE<T>* opengl_triangulated_surface=new OPENGL_TRIANGULATED_SURFACE<T>(*surface,false,
            OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Red()),OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Blue()));
        if(triangulated_surface_highlight_boundary) opengl_triangulated_surface->highlight_boundary=true;
        world.Bind_Key('0'+number,{[opengl_triangulated_surface](){static bool is_two_sided=false;is_two_sided=!is_two_sided;opengl_triangulated_surface->Set_Two_Sided(is_two_sided);},"Toggle Two Sided"});
        opengl_triangulated_surface->selectable=true;
        world.Add_Object(opengl_triangulated_surface,true,true);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Phi_File
//#################################################################
template<class T> void Add_Phi_File(const std::string& filename,OPENGL_WORLD<T>& world,int number)
{
    try{
        LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* surface;
        Create_From_File(filename,surface);
        ARRAY<T,VECTOR<int,3> >& phi=surface->levelset.phi;
        LOG::cout<<filename<<" statistics:"<<std::endl;
        LOG::cout<<"  grid = "<<surface->levelset.grid<<std::endl;
        LOG::cout<<"  phi array bounds = "<<phi.domain.min_corner.x<<" "<<phi.domain.max_corner.x<<", "<<phi.domain.min_corner.y<<" "<<phi.domain.max_corner.y<<", "<<phi.domain.min_corner.z<<" "<<phi.domain.max_corner.z<<std::endl;
        for(int i=phi.domain.min_corner.x;i<phi.domain.max_corner.x;i++)
            for(int j=phi.domain.min_corner.y;j<phi.domain.max_corner.y;j++)
                if(phi(i,j,phi.domain.min_corner.z)<=0){
                    LOG::cout<<"  phi<=0 on domain.min_corner.z"<<std::endl;
                    goto check1_end;}
      check1_end:
        for(int i=phi.domain.min_corner.x;i<phi.domain.max_corner.x;i++)
            for(int j=phi.domain.min_corner.y;j<phi.domain.max_corner.y;j++)
                if(phi(i,j,phi.domain.max_corner.z-1)<=0){
                    LOG::cout<<"  phi<=0 on domain.max_corner.z"<<std::endl;
                    goto check2_end;}
      check2_end:
        for(int i=phi.domain.min_corner.x;i<phi.domain.max_corner.x;i++)
            for(int k=phi.domain.min_corner.z;k<phi.domain.max_corner.z;k++)
                if(phi(i,phi.domain.min_corner.y,k)<=0){
                    LOG::cout<<"  phi<=0 on domain.min_corner.y"<<std::endl;
                    goto check3_end;}
      check3_end:
        for(int i=phi.domain.min_corner.x;i<phi.domain.max_corner.x;i++)
            for(int k=phi.domain.min_corner.z;k<phi.domain.max_corner.z;k++)
                if(phi(i,phi.domain.max_corner.y-1,k)<=0){
                    LOG::cout<<"  phi<=0 on domain.max_corner.y"<<std::endl;
                    goto check4_end;}
      check4_end:
        for(int j=phi.domain.min_corner.y;j<phi.domain.max_corner.y;j++)
            for(int k=phi.domain.min_corner.z;k<phi.domain.max_corner.z;k++)
                if(phi(phi.domain.min_corner.x,j,k)<=0){
                    LOG::cout<<"  phi<=0 on domain.min_corner.x"<<std::endl;
                    goto check5_end;}
      check5_end:
        for(int j=phi.domain.min_corner.y;j<phi.domain.max_corner.y;j++)
            for(int k=phi.domain.min_corner.z;k<phi.domain.max_corner.z;k++)
                if(phi(phi.domain.max_corner.x-1,j,k)<=0){
                    LOG::cout<<"  phi<=0 on domain.max_corner.x"<<std::endl;
                    goto check6_end;}
      check6_end:
        OPENGL_LEVELSET_MULTIVIEW<T>* opengl_surface=new OPENGL_LEVELSET_MULTIVIEW<T>;
        opengl_surface->Set_Levelset(surface->levelset);
        opengl_surface->Generate_Triangulated_Surface(false,"");
        opengl_surface->Set_Surface_Material(OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.8),float(.7),float(1))),OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.8),float(.8),float(.1))));
        opengl_surface->selectable=true;
        world.Add_Object(opengl_surface);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Phi2D_File
//#################################################################
template<class T> void Add_Phi2D_File(const std::string& filename,OPENGL_WORLD<T>& world,int number)
{
    try{
        LEVELSET_IMPLICIT_OBJECT<VECTOR<T,2> >* area;
        Create_From_File(filename,area);
        ARRAY<T,VECTOR<int,2> >& phi=area->levelset.phi;
        LOG::cout<<filename<<" statistics:"<<std::endl;
        LOG::cout<<"  grid = "<<area->levelset.grid<<std::endl;
        LOG::cout<<"  phi array bounds = "<<phi.domain.min_corner.x<<" "<<phi.domain.max_corner.x<<", "<<phi.domain.min_corner.y<<" "<<phi.domain.max_corner.y<<std::endl;
        for(int i=phi.domain.min_corner.x;i<phi.domain.max_corner.x;i++)
            if(phi(i,phi.domain.min_corner.y)<=0){
                LOG::cout<<"  phi<=0 on domain.min_corner.y"<<std::endl;
                goto check1_end;}
      check1_end:
        for(int i=phi.domain.min_corner.x;i<phi.domain.max_corner.x;i++)
            if(phi(i,phi.domain.max_corner.y-1)<=0){
                LOG::cout<<"  phi<=0 on domain.max_corner.y"<<std::endl;
                goto check2_end;}
      check2_end:
        for(int j=phi.domain.min_corner.y;j<phi.domain.max_corner.y;j++)
            if(phi(phi.domain.min_corner.x,j)<=0){
                LOG::cout<<"  phi<=0 on domain.min_corner.x"<<std::endl;
                goto check5_end;}
      check5_end:
        for(int j=phi.domain.min_corner.y;j<phi.domain.max_corner.y;j++)
            if (phi(phi.domain.max_corner.x-1,j)<=0){
                LOG::cout<<"  phi<=0 on domain.max_corner.x"<<std::endl;
                goto check6_end;}
      check6_end:
        SEGMENTED_CURVE_2D<T>* curve=DUALCONTOUR_2D<T>::Create_Segmented_Curve_From_Levelset(area->levelset);
        curve->Update_Bounding_Box();
        LOG::cout<<"bounding box: "<<*curve->bounding_box<<std::endl;
        OPENGL_SEGMENTED_CURVE_2D<T>* og_curve=new OPENGL_SEGMENTED_CURVE_2D<T>(*curve,OPENGL_COLOR::Yellow());
        og_curve->selectable=true;
        world.Add_Object(og_curve);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Curve_File
//#################################################################
template<class T> void Add_Curve_File(const std::string& filename,OPENGL_WORLD<T>& world,int number)
{
    try{
        SEGMENTED_CURVE<VECTOR<T,3> >* curve;
        Create_From_File(filename,curve);
        curve->Update_Bounding_Box();
        LOG::cout<<"bounding box: "<<*curve->bounding_box<<std::endl;
        OPENGL_SEGMENTED_CURVE_3D<T>* og_curve=new OPENGL_SEGMENTED_CURVE_3D<T>(*curve,OPENGL_COLOR::Yellow());
        og_curve->selectable=true;
        world.Add_Object(og_curve);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Curve2D_File
//#################################################################
template<class T> void Add_Curve2D_File(const std::string& filename,OPENGL_WORLD<T>& world,int number)
{
    try{
        SEGMENTED_CURVE_2D<T>* curve;
        Create_From_File(filename,curve);
        curve->Update_Bounding_Box();
        LOG::cout<<"bounding box: "<<*curve->bounding_box<<std::endl;
        OPENGL_SEGMENTED_CURVE_2D<T>* og_curve=new OPENGL_SEGMENTED_CURVE_2D<T>(*curve,OPENGL_COLOR::Yellow());
        og_curve->selectable=true;
        world.Add_Object(og_curve);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
// Function Add_Tet_File
//#################################################################
template<class T> void Add_Tet_File(const std::string& filename,OPENGL_WORLD<T>& world,int number)
{
    try{
        TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume;
        Create_From_File(filename,tetrahedralized_volume);
        {LOG::SCOPE scope("mesh statistics","mesh statistics");
            LOG::cout<<"filename = "<<filename<<std::endl;
            if(print_statistics) tetrahedralized_volume->Print_Statistics(LOG::cout);}
        OPENGL_TETRAHEDRALIZED_VOLUME<T>* tets=new OPENGL_TETRAHEDRALIZED_VOLUME<T>(&(tetrahedralized_volume->mesh),&(tetrahedralized_volume->particles),
            OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.9),float(.1),float(.1))),OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.1),float(.9),float(.1))));
        T ymin,ymax;
        int range=100;
        tets->subset.Preallocate(tets->mesh->elements.m);
        ymin=FLT_MAX;ymax=-FLT_MAX;
        for(int p=0;p<tets->particles->Size();p++){
            ymin=PhysBAM::min(ymin,tets->particles->X(p).y);
            ymax=PhysBAM::max(ymax,tets->particles->X(p).y);}
        world.Bind_Key('c',{[tets,ymin,ymax,range](){
                    static int cut=range;
                    tets->boundary_only=false;
                    tets->subset.Remove_All();
                    cut--;
                    if(cut<0)cut=range;
                    T y=ymin+(ymax-ymin)*cut/range;
                    for(int t=0;t<tets->mesh->elements.m;t++){
                        int i,j,k,l;
                        tets->mesh->elements(t).Get(i,j,k,l);
                        if(tets->particles->X(i).y<y || tets->particles->X(j).y<y || tets->particles->X(k).y<y || tets->particles->X(l).y<y)
                            tets->subset.Append(t);}
                },"Cross section"});
        tets->selectable=true;
        world.Add_Object(tets,true,true);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
