//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_TETRAHEDRA.h>
#include <iostream>

using namespace PhysBAM;

typedef float RW;
std::string output_directory="mt";

template<class TV>
GRID<TV>* Global_Grid(GRID<TV>* grid_in=0)
{
    static GRID<TV>* grid=0;
    GRID<TV>* old_grid=grid;
    if(grid_in) grid=grid_in;
    return old_grid;
}

//#################################################################################################################################################
// Debug Particles ################################################################################################################################
//#################################################################################################################################################

template<class TV> DEBUG_PARTICLES<TV>& Get_Debug_Particles()
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

template<class TV>
void Dump_Element(VECTOR<TV,2> X)
{
    typedef typename TV::SCALAR T;
    Add_Debug_Object(X-(X(1)-X(0)).Orthogonal_Vector().Normalized()*0.01,VECTOR<T,3>(1,0,0));
    Add_Debug_Object(X+(X(1)-X(0)).Orthogonal_Vector().Normalized()*0.01,VECTOR<T,3>(0,1,0));
}

template<class TV>
void Dump_Element(VECTOR<TV,3> X)
{
    typedef typename TV::SCALAR T;
    Add_Debug_Object(X,VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,1,0));
}

template<class TV>
void Test()
{
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename TV::SCALAR T;
    typedef typename MARCHING_TETRAHEDRA<TV>::T_FACE T_FACE;

    Get_Debug_Particles<TV>().debug_particles.template Add_Array<T>(ATTRIBUTE_ID_DISPLAY_SIZE);
    GRID<TV> grid(TV_INT()+1,RANGE<TV>(TV(),TV()+1),true);
    Global_Grid(&grid);

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);
    FILE_UTILITIES::Write_To_File<RW>(output_directory+"/common/grid.gz",grid);

    VECTOR<TV,TV::m+1> X;
    VECTOR<T,TV::m+1> phi;
    for(int i=0;i<TV::m;i++) X(i)(i)=1;

    for(int s=0;s<2;s++)
    for(int v=0;v<TV::m+1;v++)
    for(int c=0;c<(1<<(TV::m+1));c++){
        for(int i=0;i<TV::m+1;i++)
            phi(i)=((c>>i)&1)?-1:1;

        for(int i=0;i<TV::m+1;i++) Add_Debug_Particle(X(i),phi(i)<0?VECTOR<T,3>(1,0,0):VECTOR<T,3>(0,1,0));

        ARRAY<T_FACE> surface;
        VECTOR<VECTOR<ARRAY<T_FACE>*,2>,TV::m+1> boundary;
        boundary(v)(s)=&surface;

        MARCHING_TETRAHEDRA<TV>::Get_Elements_For_Tetrahedron(surface,boundary,phi,X);
        
        for(int i=0;i<surface.m;i++)
            Dump_Element(surface(i).X);
        char buff[100];
        sprintf(buff,"sign: %i vertex: %i case: %i",s,v,c);
        Flush_Frame<T,TV>(buff);
    }

    LOG::Finish_Logging();
}

//#################################################################################################################################################
// Main ###########################################################################################################################################
//#################################################################################################################################################

int main(int argc,char* argv[])
{
    bool use_3d=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-3d",&use_3d,"Use 3D");
    parse_args.Parse(true);
    
    if(use_3d) Test<VECTOR<double,3> >();
    else Test<VECTOR<double,2> >();

    return 0;
}
