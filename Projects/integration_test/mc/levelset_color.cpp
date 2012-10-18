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

using namespace PhysBAM;

typedef float RW;
std::string output_directory="levelset";

typedef double T;
typedef VECTOR<T,2> TV;
typedef VECTOR<int,2> TV_INT;

GRID<TV>* Global_Grid(GRID<TV>* grid_in=0)
{
    static GRID<TV>* grid=0;
    GRID<TV>* old_grid=grid;
    if(grid_in) grid=grid_in;
    return old_grid;
}

typedef VECTOR<double,3> TV3;
TV3 color_map[4]={TV3(0,0.7,0),TV3(0.8,0.8,0),TV3(0,0.4,1),TV3(0.8,0.2,0)};

DEBUG_PARTICLES<TV>& Get_Debug_Particles()
{
    static DEBUG_PARTICLES<TV> debug_particles;
    return debug_particles;
}

void Dump_Frame(const ARRAY<T,FACE_INDEX<TV::m> >& u,const char* title)
{
    static int frame=0;
    char buff[100];
    sprintf(buff, "%s/%i", output_directory.c_str(), frame);
    FILE_UTILITIES::Create_Directory(buff);
    FILE_UTILITIES::Write_To_File<RW>((std::string)buff+"/mac_velocities.gz",u);
    if(title) FILE_UTILITIES::Write_To_Text_File((std::string)buff+"/frame_title",title);
    Get_Debug_Particles().Write_Debug_Particles(STREAM_TYPE((RW())),output_directory,frame);
    frame++;
}

void Flush_Frame(const char* title)
{
    Dump_Frame(ARRAY<T,FACE_INDEX<TV::m> >(*Global_Grid()),title);
}

void Dump_Interface(const HASHTABLE<VECTOR<int,2>,SEGMENTED_CURVE_2D<T>*>& surface)
{
    for(typename HASHTABLE<VECTOR<int,2>,SEGMENTED_CURVE_2D<T>*>::CONST_ITERATOR it(surface);it.Valid();it.Next()){
        const VECTOR<int,2>& color_pair=it.Key();
        const SEGMENTED_CURVE_2D<T>& surf=*it.Data();
        for(int i=0;i<surf.mesh.elements.m;i++){
            const SEGMENT_2D<T>& x=surf.Get_Element(i);
            Add_Debug_Object(x.X+x.Normal()*(T).03*Global_Grid()->dX.Min(),color_map[color_pair.x]);
            Add_Debug_Object(x.X-x.Normal()*(T).03*Global_Grid()->dX.Min(),color_map[color_pair.y]);}}
}

struct ANALYTIC_LS
{
    TV n;
    T r;
    void Initialize(){for(int i=0;i<TV::m;i++) n(i)=i+M_PI/(i+M_PI);n.Normalize();r=(T)1/M_PI;}
    T Phi_Value(const TV& X){TV x=X-0.5;T s=x.Magnitude()-r;return (s<0)?min(abs(s),abs(x.Dot(n))):abs(s);}
    int Phi_Color(const TV& X){TV x=X-0.5;return (x.Magnitude()-r)<0?((x.Dot(n)<0)?2:1):0;}
};

int main(int argc, char* argv[])
{
    int resolution=32;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-o",&output_directory,"output","output directory");
    parse_args.Add("-resolution",&resolution,"res","resolution");
    parse_args.Parse();

    TV_INT counts=TV_INT()+resolution;
    GRID<TV> grid(counts,RANGE<TV>(TV(),TV()+1),true);

    Global_Grid(&grid);

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",grid);

    Get_Debug_Particles().debug_particles.Add_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);

    ARRAY<T,TV_INT> phi_value(grid.Node_Indices());
    ARRAY<int,TV_INT> phi_color(grid.Node_Indices());
    
    ANALYTIC_LS als;
    als.Initialize();

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next()){
        phi_value(it.index)=als.Phi_Value(it.Location());
        phi_color(it.index)=als.Phi_Color(it.Location());}

    MARCHING_CUBES_COLOR<TV>::Initialize_Case_Table();

    for(int i=0;i<100;i++){
        HASHTABLE<VECTOR<int,2>,SEGMENTED_CURVE_2D<T>*> surface;
        HASHTABLE<int,SEGMENTED_CURVE_2D<T>*> boundary;
        MARCHING_CUBES_COLOR<TV>::Get_Elements(grid,surface,boundary,phi_color,phi_value,i);
        Dump_Interface(surface);
        char buffer[100];
        sprintf(buffer, "newton step %i",i);
        Flush_Frame(buffer);}
    
    return 0;
}
