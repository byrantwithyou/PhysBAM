//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <Tools/Data_Structures/PAIR.h>
#include <Tools/Data_Structures/TRIPLE.h>
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <iomanip>

using namespace PhysBAM;

typedef float RW;

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
void Dump_Element(const SEGMENT_2D<T>& x,const int& color0,const int& color1)
{
    static RANDOM_NUMBERS<T> rand;
    Add_Debug_Object(x.X+x.Normal()*(T).01*rand.Get_Uniform_Number(0,1)*Global_Grid<TV>()->dX.Min(),color_map[color0]);
    Add_Debug_Object(x.X-x.Normal()*(T).01*rand.Get_Uniform_Number(0,1)*Global_Grid<TV>()->dX.Min(),color_map[color1]);
}

template<class T,class TV>
void Dump_Element(const SEGMENT_2D<T>& x,const int& color0)
{
    static RANDOM_NUMBERS<T> rand;
    Add_Debug_Object(x.X+x.Normal()*(T).005*rand.Get_Uniform_Number(-1,1)*Global_Grid<TV>()->dX.Min(),color_map[color0]);
}

template<class T,class TV>
void Dump_Element(const TRIANGLE_3D<T>& x,const int& color0,const int& color1)
{
    Add_Debug_Object(x.X,color_map[color0],color_map[color1]);
}

template<class T,class TV>
void Dump_Element(const TRIANGLE_3D<T>& x,const int& color0)
{
    Add_Debug_Object(x.X,color_map[color0],color_map[color0]);
}

template<class T,class TV,class T_FACE,class T_ELEMENT>
void Dump_Interface(const ARRAY<T_ELEMENT>& interface)
{
    for(int e=0;e<interface.m;e++)
        Dump_Element<T,TV>(interface(e).face,interface(e).color_pair.x,interface(e).color_pair.y);
}

template<class T,class TV,class T_FACE,class T_ELEMENT>
void Dump_Boundary(const ARRAY<T_ELEMENT>& boundary)
{
    for(int e=0;e<boundary.m;e++)
        Dump_Element<T,TV>(boundary(e).face,boundary(e).color);
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
    typedef typename MARCHING_CUBES_COLOR<TV>::CELL_ELEMENTS CELL_ELEMENTS;

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
                    // for(int i=0;i<TV::m;i++)
                        // n1(i)=i+M_PI/(i+M_PI);
                    // n1.Normalize();
                    n1(0)=1;
                    n2=n1.Orthogonal_Vector();}
                T Phi_Value(const TV& X) const {TV x=X-0.5223;return n1.Dot(x)>0?n1.Dot(x):(min(abs(n1.Dot(x)),abs(n2.Dot(x))));}
                int Phi_Color(const TV& X) const {TV x=X-0.5223;return n1.Dot(x)>0?0:(n2.Dot(x)>0?1:2);}
            };
            als=new ANALYTIC_LS_1;
            break;}
        default: PHYSBAM_FATAL_ERROR();}

    als->Initialize();
    for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        phi_value(it.index)=als->Phi_Value(it.Location());
        phi_color(it.index)=als->Phi_Color(it.Location());}

    MARCHING_CUBES_COLOR<TV>::Initialize_Case_Table();

    for(int i=0;i<=iterations;i++){
        HASHTABLE<TV_INT,CELL_ELEMENTS> index_to_cell_elements;
        MARCHING_CUBES_COLOR<TV>::Get_Elements(index_to_cell_elements,grid,phi_color,phi_value,i,verbose);
        for(typename HASHTABLE<TV_INT,CELL_ELEMENTS>::CONST_ITERATOR it(index_to_cell_elements);it.Valid();it.Next()){
            const CELL_ELEMENTS& cell_elements=it.Data();
            Dump_Interface<T,TV,T_FACE>(cell_elements.interface);}
            // Dump_Boundary<T,TV,T_FACE>(cell_elements.boundary);}
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
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);
    LOG::cout<<std::setprecision(16);

    if(use_3d) Build_Surface<double,VECTOR<double,3> >(argc,argv,parse_args);
    else       Build_Surface<double,VECTOR<double,2> >(argc,argv,parse_args);

    return 0;
}
