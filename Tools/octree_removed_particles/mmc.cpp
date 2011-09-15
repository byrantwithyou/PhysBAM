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
#include "../../Public_Library/Collisions_And_Interactions/COLLISION_BODY_LIST_3D.h"
#include "../../Public_Library/Deformable_Objects/DEFORMABLE_OBJECT_LIST_3D.h"
#include "../../Public_Library/Rigid_Bodies/RIGID_BODY_LIST_3D.h"
#include <Level_Sets/OCTREE_REMOVED_PARTICLE_PROCESSING.h>

using namespace PhysBAM;

typedef float T;
typedef float RW;

void Motion_By_Mean_Curvature(OCTREE_LEVELSET<T>& levelset,T fmm_band,const T cfl_number,T& time)
{
    T minimum_cell_size=levelset.grid.Get_Minimum_Cell_Size();
    T dt_curvature=levelset.sigma*(6/sqr(minimum_cell_size));
    T dt=cfl_number/dt_curvature;

    std::cout << "motion by mean curvature: stepping by dt=" << dt << std::endl;
    ARRAY<T> phi_ghost(levelset.grid.number_of_nodes);levelset.boundary->Fill_Ghost_Cells(levelset.grid,levelset.phi,phi_ghost,time); 
    ARRAY<VECTOR<T,3> >& node_locations=levelset.grid.Node_Locations();
    levelset.Compute_Curvature(time);
    T one_over_two_dx=1/(2*minimum_cell_size),one_over_two_dy=1/(2*minimum_cell_size),one_over_two_dz=1/(2*minimum_cell_size);
    VECTOR<T,3> x_offset(minimum_cell_size,0,0),y_offset(0,minimum_cell_size,0),z_offset(0,0,minimum_cell_size);
    for(int i=1;i<=levelset.grid.number_of_nodes;i++){
        T phix=(levelset.grid.Interpolate_Nodes(phi_ghost,node_locations(i)+x_offset)-levelset.grid.Interpolate_Nodes(phi_ghost,node_locations(i)-x_offset))*one_over_two_dx;
        T phiy=(levelset.grid.Interpolate_Nodes(phi_ghost,node_locations(i)+y_offset)-levelset.grid.Interpolate_Nodes(phi_ghost,node_locations(i)-y_offset))*one_over_two_dy;
        T phiz=(levelset.grid.Interpolate_Nodes(phi_ghost,node_locations(i)+z_offset)-levelset.grid.Interpolate_Nodes(phi_ghost,node_locations(i)-z_offset))*one_over_two_dz;
        T potential=levelset.phi(i)-dt*levelset.sigma*(*levelset.curvature)(i)*sqrt(sqr(phix)+sqr(phiy)+sqr(phiz));
        levelset.phi(i)=potential;
    }

    std::cout << "motion by mean curvature: calling fmm with bandwidth=" << fmm_band << " cells" << std::endl;
    levelset.Fast_Marching_Method(time,fmm_band*minimum_cell_size);

    time+=dt;
}

int main(int argc,char *argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_Double_Argument("-cfl",1,"cfl");
    parse_args.Add_Integer_Argument("-steps",10,"steps");
    parse_args.Add_Double_Argument("-band",5,"fmm band");
    parse_args.Add_Double_Argument("-sigma",.0001,"fmm band");
    parse_args.Add_String_Argument("-o","output.phi");
    parse_args.Add_Option_Argument("-substeps");
    parse_args.Set_Extra_Arguments(2, "<grid filename> <phi filename>");
    parse_args.Parse(argc,argv);
    if (parse_args.Num_Extra_Args() != 2) 
    {
        parse_args.Print_Usage();
        return 1;
    }

    std::string grid_filename=parse_args.Extra_Arg(1);
    std::string phi_filename=parse_args.Extra_Arg(2);
    std::string output_filename=parse_args.Get_String_Value("-o");
    bool write_substeps=parse_args.Get_Option_Value("-substeps");

    std::cout << "Reading grid from " << grid_filename << std::endl;
    OCTREE_GRID<T> grid;
    FILE_UTILITIES::Read_From_File<T>(grid_filename,grid);

    std::cout << "Reading phi from " << phi_filename << std::endl;
    ARRAY<T> phi;
    FILE_UTILITIES::Read_From_File<T>(phi_filename,phi);

    OCTREE_LEVELSET<T> levelset(grid,phi);
    levelset.Set_Curvature_Motion((T)parse_args.Get_Double_Value("-sigma"));

    T fmm_band=(T)parse_args.Get_Double_Value("-band");
    T cfl=(T)parse_args.Get_Double_Value("-cfl");
    T time=0;
    int total_steps=parse_args.Get_Integer_Value("-steps");
    for(int substep=1;substep<=total_steps;substep++){
        std::cout << "time = " << time << std::endl;
        Motion_By_Mean_Curvature(levelset,fmm_band,cfl,time); 
        if(write_substeps){
            std::string filename;
            if(FILE_UTILITIES::Is_Animated(output_filename)) filename=FILE_UTILITIES::Get_Frame_Filename(output_filename,substep);
            else filename=STRING_UTILITIES::string_sprintf("%s_%d",output_filename.c_str(),substep);
            std::cout << "Writing phi to " << filename << std::endl;
            FILE_UTILITIES::Write_To_File<T>(filename,phi);}
    }

    if(!write_substeps){
        std::cout << "Writing phi to " << output_filename << std::endl;
        FILE_UTILITIES::Write_To_File<T>(output_filename,phi);}
}
