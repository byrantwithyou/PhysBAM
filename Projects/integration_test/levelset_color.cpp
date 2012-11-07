//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

using namespace PhysBAM;

typedef float RW;
typedef HASHTABLE<VECTOR<int,2>,INTERVAL<int> > HASH_INTERFACE;
typedef HASHTABLE<int,INTERVAL<int> > HASH_BOUNDARY;

std::string output_directory="levelset";

template<class TV>
GRID<TV>* Global_Grid(GRID<TV>* grid_in=0)
{
    static GRID<TV>* grid=0;
    GRID<TV>* old_grid=grid;
    if(grid_in) grid=grid_in;
    return old_grid;
}

typedef VECTOR<double,3> TV3;
TV3 color_map[4]={TV3(0,0.7,0),TV3(0.8,0.8,0),TV3(0,0.4,1),TV3(0.8,0.2,0)};

template<class TV> 
DEBUG_PARTICLES<TV>& Get_Debug_Particles()
{
    static DEBUG_PARTICLES<TV> debug_particles;
    return debug_particles;
}

template<class T, class TV>
void Dump_Frame(const ARRAY<T,FACE_INDEX<TV::m> >& u,const char* title)
{
    static int frame=0;
    char buff[100];
    sprintf(buff, "%s/%i", output_directory.c_str(), frame);
    FILE_UTILITIES::Create_Directory(buff);
    FILE_UTILITIES::Write_To_File<RW>((std::string)buff+"/mac_velocities.gz",u);
    if(title) FILE_UTILITIES::Write_To_Text_File((std::string)buff+"/frame_title",title);
    Get_Debug_Particles<TV>().Write_Debug_Particles(STREAM_TYPE((RW())),output_directory,frame);
    frame++;
}

template<class T,class TV>
void Flush_Frame(const char* title)
{
    Dump_Frame<T,TV>(ARRAY<T,FACE_INDEX<TV::m> >(*Global_Grid<TV>()),title);
}

template<class T,class TV>
void Dump_Element(const VECTOR<int,2>& color_pair,SEGMENT_2D<T>& x)
{
    Add_Debug_Object(x.X+x.Normal()*(T).005*Global_Grid<TV>()->dX.Min(),color_map[color_pair.x]);
    Add_Debug_Object(x.X-x.Normal()*(T).005*Global_Grid<TV>()->dX.Min(),color_map[color_pair.y]);
}

template<class T,class TV>
void Dump_Element(const VECTOR<int,2>& color_pair,TRIANGLE_3D<T>& x)
{
    Add_Debug_Object(x.X,color_map[color_pair.x],color_map[color_pair.y]);
}

template<class T,class TV,class T_SURFACE,class T_FACE>
void Dump_Interface(const HASHTABLE<VECTOR<int,2>,T_SURFACE*>& interface)
{
    for(typename HASHTABLE<VECTOR<int,2>,T_SURFACE*>::CONST_ITERATOR it(interface);it.Valid();it.Next()){
        const VECTOR<int,2>& color_pair=it.Key();
        const T_SURFACE& surf=*it.Data();
        for(int i=0;i<surf.mesh.elements.m;i++){
            T_FACE x=surf.Get_Element(i);
            Dump_Element<T,TV>(color_pair,x);}}
}

template<class T,class TV,class T_SURFACE,class T_FACE>
void Dump_Boundary(const HASHTABLE<int,T_SURFACE*>& boundary)
{
    for(typename HASHTABLE<int,T_SURFACE*>::CONST_ITERATOR it(boundary);it.Valid();it.Next()){
        VECTOR<int,2> color_pair(it.Key(),it.Key());
        const T_SURFACE& surf=*it.Data();
        for(int i=0;i<surf.mesh.elements.m;i++){
            T_FACE x=surf.Get_Element(i);
            Dump_Element<T,TV>(color_pair,x);}}
}

template<class T,class TV>
struct ANALYTIC_LS
{
    virtual ~ANALYTIC_LS() {}
    virtual void Initialize()=0;
    virtual T Phi_Value(const TV& X)const=0;
    virtual int Phi_Color(const TV& X)const=0;
};

template<class T,class TV>
void Build_Surface(int argc,char* argv[],PARSE_ARGS& parse_args)
{
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;

    int resolution=32;
    int iterations=100;
    int test_number=0;
    bool opt_arg=false;
    bool verbose=false;
    parse_args.Extra_Optional(&test_number,&opt_arg,"example number","example number to run");
    parse_args.Add("-o",&output_directory,"output","output directory");
    parse_args.Add("-resolution",&resolution,"res","resolution");
    parse_args.Add("-iterations",&iterations,"iter","iterations");
    parse_args.Add("-verbose",&verbose,"verbose");
    parse_args.Parse();

    TV_INT counts=TV_INT()+resolution;
    GRID<TV> grid(counts,RANGE<TV>(TV(),TV()+1),true);

    Global_Grid<TV>(&grid);

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",grid);

    Get_Debug_Particles<TV>().debug_particles.template Add_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);

    ARRAY<T,TV_INT> phi_value(grid.Node_Indices());
    ARRAY<int,TV_INT> phi_color(grid.Node_Indices());
    
    ANALYTIC_LS<T,TV>* als=0;

    switch(test_number){
        case 0:{
            struct ANALYTIC_LS_0:public ANALYTIC_LS<T,TV>
            {
                TV n;
                T r;
                void Initialize(){for(int i=0;i<TV::m;i++) n(i)=i+M_PI/(i+M_PI);n.Normalize();r=(T)1/M_PI;}
                T Phi_Value(const TV& X) const {TV x=X-0.5;T s=x.Magnitude()-r;return (s<0)?min(abs(s),abs(x.Dot(n))):abs(s);}
                int Phi_Color(const TV& X) const {TV x=X-0.5;return (x.Magnitude()-r)<0?((x.Dot(n)<0)?2:1):0;}
            };
            als=new ANALYTIC_LS_0;
            break;}
        case 1:{
            struct ANALYTIC_LS_1:public ANALYTIC_LS<T,TV>
            {
                TV n1,n2;
                void Initialize(){
                    for(int i=0;i<TV::m;i++)
                        n1(i)=i+M_PI/(i+M_PI);
                    n1.Normalize();
                    n2=n1.Orthogonal_Vector();}
                T Phi_Value(const TV& X) const {TV x=X-0.5223;return n1.Dot(x)>0?n1.Dot(x):(min(abs(n1.Dot(x)),abs(n2.Dot(x))));}
                int Phi_Color(const TV& X) const {TV x=X-0.5223;return n1.Dot(x)>0?0:(n2.Dot(x)>0?1:2);}
            };
            als=new ANALYTIC_LS_1;
            break;}
        default: PHYSBAM_FATAL_ERROR();}

    als->Initialize();
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next()){
        phi_value(it.index)=als->Phi_Value(it.Location());
        phi_color(it.index)=als->Phi_Color(it.Location());}

    MARCHING_CUBES_COLOR<TV>::Initialize_Case_Table();

    for(int i=0;i<=iterations;i++){
        HASHTABLE<TV_INT,PAIR<HASH_INTERFACE,HASH_BOUNDARY> > cell_to_element;
        HASHTABLE<VECTOR<int,2>,T_SURFACE*> surface;
        HASHTABLE<int,T_SURFACE*> boundary;
        MARCHING_CUBES_COLOR<TV>::Get_Elements(grid,surface,boundary,cell_to_element,phi_color,phi_value,i,verbose);
        Dump_Interface<T,TV,T_SURFACE,T_FACE>(surface);
        Dump_Boundary<T,TV,T_SURFACE,T_FACE>(boundary);
        char buffer[100];
        sprintf(buffer, "newton step %i",i);
        Flush_Frame<T,TV>(buffer);}

    delete als;
}

int main(int argc, char* argv[])
{
    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"Use 3D");
    parse_args.Parse(true);

    if(use_3d) Build_Surface<double,VECTOR<double,3> >(argc,argv,parse_args);
    else       Build_Surface<double,VECTOR<double,2> >(argc,argv,parse_args);

    return 0;
}
