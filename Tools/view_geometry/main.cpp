//#####################################################################
// Copyright 2002-2012, Robert Bridson, Eran Guendelman, Eilene Hao, Geoffrey Irving, Cynthia Lau, Neil Molino, Avi Robinson-Mosher, Andrew Selle, Eftychios Sifakis, Alexey Stomakhin, Jerry Talton, Joseph Teran, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Particles/PARTICLES_FORWARD.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Grids_Uniform_Computations/DUALCONTOUR_2D.h>
#include <Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_AXES.h>
#include <OpenGL/OpenGL/OPENGL_BASIC_CALLBACKS.h>
#include <OpenGL/OpenGL/OPENGL_BOX_3D.h>
#include <OpenGL/OpenGL/OPENGL_CALLBACK.h>
#include <OpenGL/OpenGL/OPENGL_HEXAHEDRALIZED_VOLUME.h>
#include <OpenGL/OpenGL/OPENGL_IMPLICIT_SURFACE.h>
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
using namespace FILE_UTILITIES;
template<class T> void Add_File(const std::string& filename,int number);
template<class T> void Add_Tri2D_File(const std::string& filename,OPENGL_WORLD<T>& world,int number);
template<class T> void Add_Tri_File(const std::string& filename,OPENGL_WORLD<T>& world,int number);
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
// Class OPENGL_CALLBACK_TOGGLE_TWO_SIDED
//#################################################################
template<class T>
class OPENGL_CALLBACK_TOGGLE_TWO_SIDED:public OPENGL_CALLBACK
{
public:
    OPENGL_TRIANGULATED_SURFACE<T>& surface;
    bool is_two_sided;

    explicit OPENGL_CALLBACK_TOGGLE_TWO_SIDED(OPENGL_TRIANGULATED_SURFACE<T>& surface_input)
        :surface(surface_input),is_two_sided(false)
    {}

    void operator()()
    {is_two_sided=!is_two_sided;
    surface.Set_Two_Sided(is_two_sided);}

    void Print(std::ostream& out){out<<"Toggle two sided";}
};
//#################################################################
// Class OPENGL_CALLBACK_CROSS_SECTION
//#################################################################
template<class T>
class OPENGL_CALLBACK_CROSS_SECTION:public OPENGL_CALLBACK
{
public:
    OPENGL_TETRAHEDRALIZED_VOLUME<T>& tets;
    T ymin,ymax;
    int cut,range;

    explicit OPENGL_CALLBACK_CROSS_SECTION(OPENGL_TETRAHEDRALIZED_VOLUME<T>& tets_input)
        :tets(tets_input),range(100)
    {
        cut=range;
        tets.subset.Preallocate(tets.mesh->elements.m);
        ymin=FLT_MAX;ymax=-FLT_MAX;
        for(int p=0;p<tets.particles->Size();p++){
            ymin=PhysBAM::min(ymin,tets.particles->X(p).y);
            ymax=PhysBAM::max(ymax,tets.particles->X(p).y);}
    }

    void operator()()
    {tets.boundary_only=false;
    tets.subset.Remove_All();
    cut--;if(cut<0)cut=range;
    T y=ymin+(ymax-ymin)*cut/range;
    for(int t=0;t<tets.mesh->elements.m;t++){
        int i,j,k,l;tets.mesh->elements(t).Get(i,j,k,l);
        if(tets.particles->X(i).y<y || tets.particles->X(j).y<y || tets.particles->X(k).y<y || tets.particles->X(l).y<y)
            tets.subset.Append(t);}}

    void Print(std::ostream& out){out<<"Cross section";}
};
//#################################################################
// Function main
//#################################################################
int main(int argc,char *argv[])
{
    bool type_double=false;    

    ARRAY<std::string> files;
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"-float")) type_double=false;
        else if (!strcmp(argv[i],"-double")) type_double=true;
        else if (!strcmp(argv[i],"-show_boundary")) triangulated_surface_highlight_boundary = true;
        else if (!strcmp(argv[i],"-nostat")) print_statistics = false;
        else files.Append(argv[i]);}

    for(int i=0;i<files.m;i++){
        if(!type_double) Add_File<float>(files(i),i);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        else Add_File<double>(files(i),i);
#else
        else{LOG::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
        }
    return 0;
}
//#################################################################
// Function Add_File
//#################################################################
template<class T> void Add_File(const std::string& filename,int number)
{
    OPENGL_WORLD<T> world;
    world.Add_Light(new OPENGL_LIGHT(VECTOR<double,3>(1,1,1),.8));
    world.Add_Light(new OPENGL_LIGHT(VECTOR<double,3>(.3,-2,.1),.4));
    world.Add_Light(new OPENGL_LIGHT(VECTOR<double,3>(-2,.1,.5),.2));
    world.Set_Ambient_Light(.2);
    world.Initialize("PhysBAM geometry viewer");

    world.Bind_Key("^c", new OPENGL_CALLBACK_SAVE_VIEW<T>(world, "camera_script", true));
    world.Bind_Key('c', new OPENGL_CALLBACK_LOAD_VIEW<T>(world, "camera_script", true));

    FILE_TYPE type=Get_File_Type(filename);
    switch(type){
        case TRI_FILE: Add_Tri_File<T>(filename,world,number);break;
        case TRI2D_FILE: Add_Tri2D_File<T>(filename,world,number);break;
        case PHI_FILE: Add_Phi_File<T>(filename,world,number);break;
        case PHI2D_FILE: Add_Phi2D_File<T>(filename,world,number);break;
        case CURVE_FILE: Add_Curve_File<T>(filename,world,number);break;
        case CURVE2D_FILE: Add_Curve2D_File<T>(filename,world,number);break;
        case TET_FILE: Add_Tet_File<T>(filename,world,number);break;
        case HEX_FILE: Add_Hex_File<T>(filename,world,number);break;
        case BOX_FILE: Add_Box_File<T>(filename,world,number);break;
        default: LOG::cerr<<"Unrecognized file "<<filename<<std::endl;}
    LOG::cout<<std::flush;
    world.Center_Camera_On_Scene();
    world.window->Main_Loop();
}
//#################################################################
// Function Add_Box_File
//#################################################################
template<class T> void Add_Box_File(const std::string& filename,OPENGL_WORLD<T>& world,int number)
{
//    typedef VECTOR<T,3> TV;
    try{
        // BOX<TV>* box=new BOX<TV>;
        // FILE_UTILITIES::Read_From_File<T>(filename,*box);
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
        FILE_UTILITIES::Create_From_File<T>(filename,hex_vol);
        OPENGL_HEXAHEDRALIZED_VOLUME<T>* ohv=new OPENGL_HEXAHEDRALIZED_VOLUME<T>(&hex_vol->mesh,&hex_vol->particles,OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.7),float(1),float(.8))),OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.7),float(8),float(.1))));
        //world.Bind_Key('0'+number,new OPENGL_CALLBACK_TOGGLE_TWO_SIDED(*ohv));
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
        FILE_UTILITIES::Create_From_File<T>(filename,area);
        TRIANGULATED_SURFACE<T>* surface=new TRIANGULATED_SURFACE<T>(area->mesh,*new DEFORMABLE_PARTICLES<VECTOR<T,3> >);
        surface->particles.Add_Elements(area->particles.Size());
        for(int p=0;p<area->particles.Size();p++)surface->particles.X(p)=VECTOR<T,3>(area->particles.X(p));
        area->Update_Bounding_Box();
        LOG::cout<<"bounding box: "<<*area->bounding_box<<std::endl;
        OPENGL_TRIANGULATED_SURFACE<T>* opengl_triangulated_surface=new OPENGL_TRIANGULATED_SURFACE<T>(*surface,false);
        world.Bind_Key('0'+number,new OPENGL_CALLBACK_TOGGLE_TWO_SIDED<T>(*opengl_triangulated_surface));
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
        FILE_UTILITIES::Create_From_File<T>(filename,surface);
        {LOG::SCOPE scope("mesh statistics","mesh statistics");
        LOG::cout<<"filename = "<<filename<<std::endl;
        if(print_statistics) surface->Print_Statistics(LOG::cout);}
        OPENGL_TRIANGULATED_SURFACE<T>* opengl_triangulated_surface=new OPENGL_TRIANGULATED_SURFACE<T>(*surface,false,
            OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Red()),OPENGL_MATERIAL::Plastic(OPENGL_COLOR::Blue()));
        if(triangulated_surface_highlight_boundary) opengl_triangulated_surface->highlight_boundary=true;
        world.Bind_Key('0'+number,new OPENGL_CALLBACK_TOGGLE_TWO_SIDED<T>(*opengl_triangulated_surface));
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
        FILE_UTILITIES::Create_From_File<T>(filename,surface);
        ARRAY<T,VECTOR<int,3> >& phi=surface->levelset.phi;
        LOG::cout<<filename<<" statistics:"<<std::endl;
        LOG::cout<<"  grid = "<<surface->levelset.grid<<std::endl;
        LOG::cout<<"  phi array bounds = "<<phi.domain.min_corner.x<<" "<<phi.domain.max_corner.x<<", "<<phi.domain.min_corner.y<<" "<<phi.domain.max_corner.y<<", "<<phi.domain.min_corner.z<<" "<<phi.domain.max_corner.z<<std::endl;
        for(int i=phi.domain.min_corner.x;i<phi.domain.max_corner.x;i++)for(int j=phi.domain.min_corner.y;j<phi.domain.max_corner.y;j++)if(phi(i,j,phi.domain.min_corner.z)<=0){LOG::cout<<"  phi<=0 on domain.min_corner.z"<<std::endl;goto check1_end;}check1_end:
        for(int i=phi.domain.min_corner.x;i<phi.domain.max_corner.x;i++)for(int j=phi.domain.min_corner.y;j<phi.domain.max_corner.y;j++)if(phi(i,j,phi.domain.max_corner.z-1)<=0){LOG::cout<<"  phi<=0 on domain.max_corner.z"<<std::endl;goto check2_end;}check2_end:
        for(int i=phi.domain.min_corner.x;i<phi.domain.max_corner.x;i++)for(int k=phi.domain.min_corner.z;k<phi.domain.max_corner.z;k++)if(phi(i,phi.domain.min_corner.y,k)<=0){LOG::cout<<"  phi<=0 on domain.min_corner.y"<<std::endl;goto check3_end;}check3_end:
        for(int i=phi.domain.min_corner.x;i<phi.domain.max_corner.x;i++)for(int k=phi.domain.min_corner.z;k<phi.domain.max_corner.z;k++)if(phi(i,phi.domain.max_corner.y-1,k)<=0){LOG::cout<<"  phi<=0 on domain.max_corner.y"<<std::endl;goto check4_end;}check4_end:
        for(int j=phi.domain.min_corner.y;j<phi.domain.max_corner.y;j++)for(int k=phi.domain.min_corner.z;k<phi.domain.max_corner.z;k++)if(phi(phi.domain.min_corner.x,j,k)<=0){LOG::cout<<"  phi<=0 on domain.min_corner.x"<<std::endl;goto check5_end;}check5_end:
        for(int j=phi.domain.min_corner.y;j<phi.domain.max_corner.y;j++)for(int k=phi.domain.min_corner.z;k<phi.domain.max_corner.z;k++)if(phi(phi.domain.max_corner.x-1,j,k)<=0){LOG::cout<<"  phi<=0 on domain.max_corner.x"<<std::endl;goto check6_end;}check6_end:
        OPENGL_LEVELSET_MULTIVIEW<T>* opengl_surface=new OPENGL_LEVELSET_MULTIVIEW<T>;
        opengl_surface->Set_Levelset(surface->levelset);
        opengl_surface->Generate_Triangulated_Surface(false,"");
        opengl_surface->Set_Surface_Material(OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.8),float(.7),float(1))),OPENGL_MATERIAL::Matte(OPENGL_COLOR(float(.8),float(.8),float(.1))));
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
        FILE_UTILITIES::Create_From_File<T>(filename,area);
        ARRAY<T,VECTOR<int,2> >& phi=area->levelset.phi;
        LOG::cout<<filename<<" statistics:"<<std::endl;
        LOG::cout<<"  grid = "<<area->levelset.grid<<std::endl;
        LOG::cout<<"  phi array bounds = "<<phi.domain.min_corner.x<<" "<<phi.domain.max_corner.x<<", "<<phi.domain.min_corner.y<<" "<<phi.domain.max_corner.y<<std::endl;
        for(int i=phi.domain.min_corner.x;i<phi.domain.max_corner.x;i++)if(phi(i,phi.domain.min_corner.y)<=0){LOG::cout<<"  phi<=0 on domain.min_corner.y"<<std::endl;goto check1_end;}check1_end:
        for(int i=phi.domain.min_corner.x;i<phi.domain.max_corner.x;i++)if(phi(i,phi.domain.max_corner.y-1)<=0){LOG::cout<<"  phi<=0 on domain.max_corner.y"<<std::endl;goto check2_end;}check2_end:
        for(int j=phi.domain.min_corner.y;j<phi.domain.max_corner.y;j++)if(phi(phi.domain.min_corner.x,j)<=0){LOG::cout<<"  phi<=0 on domain.min_corner.x"<<std::endl;goto check5_end;}check5_end:
        for(int j=phi.domain.min_corner.y;j<phi.domain.max_corner.y;j++)if(phi(phi.domain.max_corner.x-1,j)<=0){LOG::cout<<"  phi<=0 on domain.max_corner.x"<<std::endl;goto check6_end;}check6_end:
        SEGMENTED_CURVE_2D<T>* curve=DUALCONTOUR_2D<T>::Create_Segmented_Curve_From_Levelset(area->levelset);
        curve->Update_Bounding_Box();
        LOG::cout<<"bounding box: "<<*curve->bounding_box<<std::endl;
        OPENGL_SEGMENTED_CURVE_2D<T>* og_curve=new OPENGL_SEGMENTED_CURVE_2D<T>(*curve,OPENGL_COLOR::Yellow());
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
        FILE_UTILITIES::Create_From_File<T>(filename,curve);
        curve->Update_Bounding_Box();
        LOG::cout<<"bounding box: "<<*curve->bounding_box<<std::endl;
        OPENGL_SEGMENTED_CURVE_3D<T>* og_curve=new OPENGL_SEGMENTED_CURVE_3D<T>(*curve,OPENGL_COLOR::Yellow());
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
        FILE_UTILITIES::Create_From_File<T>(filename,curve);
        curve->Update_Bounding_Box();
        LOG::cout<<"bounding box: "<<*curve->bounding_box<<std::endl;
        OPENGL_SEGMENTED_CURVE_2D<T>* og_curve=new OPENGL_SEGMENTED_CURVE_2D<T>(*curve,OPENGL_COLOR::Yellow());
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
        FILE_UTILITIES::Create_From_File<T>(filename,tetrahedralized_volume);
        {LOG::SCOPE scope("mesh statistics","mesh statistics");
        LOG::cout<<"filename = "<<filename<<std::endl;
        if(print_statistics) tetrahedralized_volume->Print_Statistics(LOG::cout);}
        OPENGL_TETRAHEDRALIZED_VOLUME<T>* tets=new OPENGL_TETRAHEDRALIZED_VOLUME<T>(&(tetrahedralized_volume->mesh),&(tetrahedralized_volume->particles),
            OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.9),float(.1),float(.1))),OPENGL_MATERIAL::Plastic(OPENGL_COLOR(float(.1),float(.9),float(.1))));
        world.Bind_Key('c',new OPENGL_CALLBACK_CROSS_SECTION<T>(*tets));
        world.Add_Object(tets,true,true);}
    catch(FILESYSTEM_ERROR&){}
}
//#################################################################
