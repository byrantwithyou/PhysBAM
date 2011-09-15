#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_Rendering/Rendering_Objects/RENDERING_IMPLICIT_SURFACE.h>
#include <fstream>
#include <iostream>
#include "../../Public_Library/Deformable_Objects/DEFORMABLE_OBJECT_LIST_3D.h"
#include "../../Public_Library/Rigid_Bodies/RIGID_BODY_LIST_3D.h"
#include <PhysBAM_Geometry/Level_Sets/IMPLICIT_SURFACE_MAKER.h>

using namespace PhysBAM;

typedef float T;
typedef float RW;

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-o", "output", "directory", "output directory");
    parse_args.Set_Extra_Arguments(1, "<frame number>");
    parse_args.Parse(argc,argv);
    if (parse_args.Num_Extra_Args() != 1) 
    {
        parse_args.Print_Usage();
        return 1;
    }
    std::string output_directory=parse_args.Get_String_Value("-o");
    int frame=atoi(parse_args.Extra_Arg(1).c_str());


    std::cout<<"Reading levelset: "<<std::endl;
    ARRAYS<VECTOR<T,3> > levelset_phi;GRID_3D<T> levelset_grid;
    LEVELSET_3D<T> levelset(levelset_grid,levelset_phi);
    FILE_UTILITIES::Read_From_File<RW>(STRING_UTILITIES::string_sprintf("levelset_input/levelset.%d",frame),levelset);
    std::cout<<levelset_grid<<std::endl;

    std::cout<<"Reading deep water texture: ";
    ARRAYS<VECTOR<T,2> > deep_water_height;GRID_2D<T> deep_water_grid;
    FILE_UTILITIES::Read_From_File<RW>("deep_water_input/grid",deep_water_grid);
    std::cout<<deep_water_grid<<std::endl;
    FILE_UTILITIES::Read_From_File<RW>(STRING_UTILITIES::string_sprintf("deep_water_input/eta.%d",frame),deep_water_height);

    GRID_3D<T> grid(600,80,400,-3,4.5,0,1,-2,3);BOX_3D<T> levelset_domain=levelset_grid.Domain();
    ARRAYS<VECTOR<T,3> > phi(grid,3);
    LINEAR_INTERPOLATION<T,T> interpolation;
    for(int i=1;i<=grid.m;i++)for(int j=1;j<=grid.n;j++)for(int ij=1;ij<=grid.mn;ij++){
        T x=grid.x(i),y=grid.y(j),z=grid.z(ij);
        T distance=levelset_domain.Signed_Distance(VECTOR_3D<T>(x,y,z));
        T clamped_phi=levelset.Phi(levelset_grid.Clamp(VECTOR_3D<T>(x,y,z)));
        T flat_phi=y-.35;
        if(distance<0)phi(i,j,ij)=clamped_phi;
        else if(distance<1)phi(i,j,ij)=distance*flat_phi+(1-distance)*clamped_phi;
        else phi(i,j,ij)=flat_phi;
    }
    std::cout<<"Fast marching"<<std::endl;
    LEVELSET_3D<T> output_levelset(grid,phi);
    output_levelset.Fast_Marching_Method();
    for(int i=1;i<=grid.m;i++)for(int j=1;j<=grid.n;j++)for(int ij=1;ij<=grid.mn;ij++){
        T x=grid.x(i)*300;
        T z=grid.z(ij)*300;
        int x_node=(int)((x-deep_water_grid.xmin)*deep_water_grid.one_over_dx);
        int y_node=(int)((z-deep_water_grid.ymin)*deep_water_grid.one_over_dy);
        while(x_node<0)x_node+=deep_water_grid.m;x_node=x_node%(deep_water_grid.m-1)+1;
        while(y_node<0)y_node+=deep_water_grid.n;y_node=y_node%(deep_water_grid.n-1)+1;
        VECTOR_2D<T> tiled_point(x_node*deep_water_grid.dx+deep_water_grid.xmin,y_node*deep_water_grid.dy+deep_water_grid.ymin);
        T tiled_deep_water=interpolation.From_Base_Node(deep_water_grid,deep_water_height,tiled_point,x_node,y_node);
        phi(i,j,ij)-=tiled_deep_water/700;
    }
    std::cout<<"Fast marching"<<std::endl;
    output_levelset.Fast_Marching_Method();

    std::cout<<"Writing merged levelset file"<<std::endl;
    FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/levelset.%d",output_directory.c_str(),frame).c_str(),output_levelset);

    return 0;
}

