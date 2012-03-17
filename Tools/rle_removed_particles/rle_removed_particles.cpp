//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <PhysBAM_Dynamics/Level_Sets/RLE_REMOVED_PARTICLES_PROCESSING.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Implicit_Objects_RLE/RLE_IMPLICIT_OBJECT.h>
using namespace PhysBAM;
//#####################################################################
// Class PROCESSOR
//#####################################################################
template<class T>
class PROCESSOR
{
public:
    typedef RLE_GRID_3D<T> T_GRID;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;

    PARSE_ARGS& parse_args;
    T sigma;

    PROCESSOR(PARSE_ARGS& parse_args_input)
        :parse_args(parse_args_input)
    {
        sigma=(T)parse_args.Get_Double_Value("-sigma");
    }

//#####################################################################
    void Fast_March(T_GRID& grid,ARRAY<T>& phi);
    void Curvature_Motion(T_GRID& grid,ARRAY<T>& phi);
    void Rebuild(T_GRID& grid,ARRAY<T>& phi,const int extra_cells);
    template<class RW> void Process(const int frame);
//#####################################################################
};
//#####################################################################
// Function Fast_March
//#####################################################################
template<class T> void PROCESSOR<T>::
Fast_March(T_GRID& grid,ARRAY<T>& phi)
{
    LEVELSET_RLE<T_GRID>(grid,phi).Fast_Marching_Method(0);
}
//#####################################################################
// Function Curvature_Motion
//#####################################################################
template<class T> void PROCESSOR<T>::
Curvature_Motion(T_GRID& grid,ARRAY<T>& phi)
{
    LOG::Time("curvature motion");
    LEVELSET_RLE<T_GRID> levelset(grid,phi);
    levelset.Set_Curvature_Motion(sigma);
    T cfl=(T)parse_args.Get_Double_Value("-curvature_cfl");
    T dt_cfl=cfl*levelset.Curvature_CFL();
    int substeps=(int)ceil(1/dt_cfl);
    LOG::cout<<"substeps = "<<substeps<<std::endl;
    T dt=1/(T)substeps;
    for(int step=0;step<substeps;step++){
        Rebuild(grid,phi,1);
        LOG::cout<<"number of cells = "<<grid.number_of_cells<<std::endl;
        levelset.Curvature_Motion(dt,0);
        levelset.Fast_Marching_Method(0);}
    Rebuild(grid,phi,1);
    LOG::Stop_Time();
}
//#####################################################################
// Function Rebuild
//#####################################################################
template<class T> void PROCESSOR<T>::
Rebuild(T_GRID& grid,ARRAY<T>& phi,const int extra_cells)
{
    LEVELSET_RLE<T_GRID>(grid,phi).Rebuild_Grid(extra_cells,0);
    for(CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++)if(!cell.Short()){int c=cell.Cell();phi(c)=phi(c+1)=grid.positive_bandwidth;}
}
//#####################################################################
// Function Process
//#####################################################################
template<class T> template<class RW> void PROCESSOR<T>::
Process(const int frame)
{
    typedef VECTOR<T,3> TV;

    LOG::SCOPE scope("PROCESSING","Processing frame %d",frame);
    bool extrapolate_into_objects=parse_args.Get_Option_Value("-extrapolate_into_objects");
    T refinement=(T)parse_args.Get_Double_Value("-refine");
    T bandwidth=(T)parse_args.Get_Double_Value("-bandwidth");
    std::string f=STRING_UTILITIES::string_sprintf(".%d",frame);

    LOG::Time("reading grid");
    T_GRID water_grid;FILE_UTILITIES::Read_From_File<RW>("rle_grid"+f,water_grid);
    LOG::Time("reading levelset");
    ARRAY<T> water_phi;FILE_UTILITIES::Read_From_File<RW>("rle_levelset"+f,water_phi);
    LOG::cout<<"water uniform grid = "<<water_grid.uniform_grid<<std::endl;

    LOG::Time("reading particles");
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> particles;
    {ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*> particles_array;
    FILE_UTILITIES::Read_From_File<RW>("removed_negative_particles"+f,particles_array);
    int number_of_particles=0;
    for(int b=0;b<particles_array.m;b++)if(particles_array(b)) number_of_particles+=particles_array(b)->number;
    LOG::cout<<"particles = "<<number_of_particles<<std::endl;
    particles.Preallocate(number_of_particles);
    for(int b=0;b<particles_array.m;b++)if(particles_array(b)) particles.Take(*particles_array(b));
    particles_array.Delete_Pointers_And_Clean_Memory();}

    RLE_REMOVED_PARTICLE_PROCESSING<T> removed_particle_processing(particles);
    removed_particle_processing.blending_parameter=(T)parse_args.Get_Double_Value("-b");
    removed_particle_processing.scale=(T)parse_args.Get_Double_Value("-scale");
    removed_particle_processing.particle_power=(T)parse_args.Get_Double_Value("-power");
    removed_particle_processing.use_velocity_scaling=parse_args.Get_Option_Value("-use_velocity_scaling");
    removed_particle_processing.preserve_volume=parse_args.Get_Option_Value("-preserve_volume");
    removed_particle_processing.dt=(T)parse_args.Get_Double_Value("-dt");

    T_GRID& particle_grid=removed_particle_processing.grid;
    ARRAY<T>& particle_phi=removed_particle_processing.phi;

    LOG::Time("determining particle grid");
    {BOX<TV> particle_box=BOX<TV>::Empty_Box();
    for(int p=0;p<particles.number;p++)particle_box.Enlarge_To_Include_Point(particles.X(p));
    particle_box.Change_Size(2*water_grid.positive_bandwidth);
    GRID_3D<T> uniform_grid=GRID_3D<T>::Create_Grid_Given_Cell_Size(particle_box,water_grid.uniform_grid.dx/refinement,false);
    // shift grid to hit origin exactly (to prevent frame incoherent jittering)
    TV shift=water_grid.uniform_grid.Domain().Minimum_Corner()-uniform_grid.Domain().Minimum_Corner();
    shift.x=fmod(shift.x,uniform_grid.dx);shift.y=fmod(shift.y,uniform_grid.dx);shift.z=fmod(shift.z,uniform_grid.dx);
    uniform_grid=GRID_3D<T>(uniform_grid.m,uniform_grid.n,uniform_grid.mn,uniform_grid.Domain()+shift);
    LOG::cout<<"particle uniform grid = "<<uniform_grid<<std::endl;
    particle_grid.Set_Uniform_Grid(uniform_grid);
    particle_grid.Set_Uniform_Grid(uniform_grid);
    particle_grid.Set_Linear_Pressure_And_Linear_Horizontal_Velocity_Profile();
    particle_grid.Set_Positive_Bandwidth_In_Cells(bandwidth);
    particle_grid.Set_Negative_Bandwidth_In_Cells(bandwidth);}

    if(extrapolate_into_objects) PHYSBAM_NOT_IMPLEMENTED();

    removed_particle_processing.Build_Grid_And_Rasterize_Particles();
    particles.Clean_Memory();

/*
    LOG::Time("writing particle grid");
    FILE_UTILITIES::Write_To_File<RW>("particle_rle_grid"+f,removed_particle_processing.grid);
    LOG::Time("writing raw particle phi");
    FILE_UTILITIES::Write_To_File<RW>("raw_particle_rle_levelset"+f,removed_particle_processing.phi);
*/

    RLE_IMPLICIT_OBJECT<TV> water(water_grid,water_phi);
    T particle_contour=particle_grid.Minimum_Edge_Length()*removed_particle_processing.blending_parameter*removed_particle_processing.particle_power;

/*
    LOG::Time("writing merged levelset");
    {ARRAY<T> merged_phi(particle_grid.number_of_cells,false);ARRAY<T>::copy(particle_grid.positive_bandwidth,merged_phi);
    for(CELL_ITERATOR cell(particle_grid,particle_grid.number_of_ghost_cells);cell;cell++)if(cell.Short()){int c=cell.Cell();TV X=cell.X();
        merged_phi(c)=water(X)+particle_phi(c);}
    FILE_UTILITIES::Write_To_File<RW>("merged_rle_levelset"+f,merged_phi);}
*/

/*
    LOG::Time("writing particle phi");
    {ARRAY<T> phi(particle_phi);phi+=particle_contour;
    Fast_March(particle_grid,phi);
    FILE_UTILITIES::Write_To_File<RW>("particle_rle_levelset"+f,phi);}
*/

/*
    LOG::Time("writing refined water phi");
    {ARRAY<T> refined_water_phi(particle_grid.number_of_cells,false);ARRAY<T>::copy(particle_grid.positive_bandwidth,refined_water_phi);
    for(CELL_ITERATOR cell(particle_grid,particle_grid.number_of_ghost_cells);cell;cell++)if(cell.Short()){int c=cell.Cell();TV X=cell.X();
        refined_water_phi(c)=water(X);}
    Fast_March(particle_grid,refined_water_phi);
    FILE_UTILITIES::Write_To_File<RW>("refined_water_rle_levelset"+f,refined_water_phi);}
*/

    {LOG::SCOPE scope("UNION","unioning");
    particle_phi+=particle_contour;
    LOG::Time("fast marching");
    Fast_March(particle_grid,particle_phi); // fast march to prevent stair stepping artifacts
    if(sigma) Curvature_Motion(particle_grid,particle_phi);
    LOG::Time("writing grid");
    FILE_UTILITIES::Write_To_File<RW>("union_rle_grid"+f,particle_grid);
    //FILE_UTILITIES::Write_To_File<RW>("smoothed_particle_rle_levelset"+f,particle_phi);
    LOG::Time("unioning");
    for(CELL_ITERATOR cell(particle_grid,particle_grid.number_of_ghost_cells);cell;cell++)if(cell.Short()){int c=cell.Cell();TV X=cell.X();
        particle_phi(c)=min(particle_phi(c),water(X));}
    LOG::Time("fast marching");
    Fast_March(particle_grid,particle_phi);
    LOG::Time("writing levelset");
    FILE_UTILITIES::Write_To_File<RW>("union_rle_levelset"+f,particle_phi);}
}
//#####################################################################
// MAIN
//#####################################################################
int main(int argc, char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

    PARSE_ARGS parse_args;
    parse_args.Add_Double_Argument("-scale",.25);
    parse_args.Add_Double_Argument("-refine",1);
    parse_args.Add_Double_Argument("-bandwidth",3);
    parse_args.Add_Option_Argument("-use_velocity_scaling");
    parse_args.Add_Option_Argument("-preserve_volume");
    parse_args.Add_Double_Argument("-dt",(double)1/60,"dt","dt (time step)");
    parse_args.Add_Option_Argument("-write_particle_levelset");
    parse_args.Add_Option_Argument("-write_refined_levelset");
    parse_args.Add_Option_Argument("-write_union");
    parse_args.Add_String_Argument("-o","","directory","output directory");
    parse_args.Add_Double_Argument("-b",.8,"blending","blending parameter");
    parse_args.Add_Double_Argument("-power",1,"power","particle power");
    parse_args.Add_Double_Argument("-sigma",0);
    parse_args.Add_Double_Argument("-curvature_cfl",1);
    parse_args.Add_Option_Argument("-extrapolate_into_objects");
    parse_args.Set_Extra_Arguments(1,"<frame number>");
    parse_args.Parse(argc,argv);

    int frame=atoi(parse_args.Extra_Arg(1).c_str());
    PROCESSOR<float>(parse_args).Process<float>(frame);
    LOG::cout<<std::endl;
    return 0;
}
//#####################################################################
